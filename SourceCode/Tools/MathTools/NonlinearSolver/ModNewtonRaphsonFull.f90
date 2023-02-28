module modNewtonRaphsonFull
    use modNonLinearSystemOfEquations
    use NonlinearSolver
    implicit none

    type, extends(ClassNonlinearSolver) :: ClassNewtonRaphsonFull
        real(8) :: tol
        integer :: itmax
        integer :: NormType = 2 , MatrixType = 2
        logical :: ShowInfo = .true.

    contains
        procedure :: Solve => NewtonRaphsonFull_Solve
        !procedure :: Constructor => NewtonRaphsonFull_Constructor
        procedure :: ReadSolverParameters => NewtonRaphsonFull_ReadSolverParameters
    end type

    type ClassNewtonRaphsonFullErrors
        integer :: MaxNumberOfIteration = 1
        integer :: UserEvaluateSystemReportedError = 2
        integer :: UserEvaluateGradientReportedError = 3
        integer :: LinearSystemError = 4
    end type

    type ClassNewtonRaphsonFullNormTypes
        integer :: L2Norm = 1
        integer :: MaximumAbsoluteValue=2
    end type

    type ClassNewtonRaphsonFullMatrixTypes
        integer :: Full = 1
        integer :: Sparse = 2
    end type


    type (ClassNewtonRaphsonFullNormTypes)   , parameter :: NewtonRaphsonFull_NormTypes = ClassNewtonRaphsonFullNormTypes()
    type (ClassNewtonRaphsonFullErrors)      , parameter :: NewtonRaphsonFull_Errors = ClassNewtonRaphsonFullErrors()
    type (ClassNewtonRaphsonFullMatrixTypes) , parameter :: NewtonRaphsonFull_MatrixTypes = ClassNewtonRaphsonFullMatrixTypes()


contains
!-----------------------------------------------------------------

    subroutine NewtonRaphsonFull_Solve(this,SOE,Xguess,X)

        use GlobalSparseMatrix
        use MathRoutines
        class(ClassNewtonRaphsonFull) :: this
        class(ClassNonLinearSystemOfEquations)      :: SOE
        real(8),dimension(:)          :: Xguess , X

        integer :: it, i
        real(8),allocatable,dimension(:) :: R , DX, DXFull
        real(8) :: normR , normR0

        real(8),dimension(:,:),pointer :: GFull
        class(ClassGlobalSparseMatrix),pointer :: GSparse


        call SOE%Status%SetSuccess
        call this%Status%SetSuccess

        it = 0
        SOE%it = it
        X=Xguess

        if (SOE%isPeriodic) then
            allocate(R(SOE%nDOF),DX(SOE%nDOF),DXFull(size(X))) !SOE%nDOF = DOF reduced system
        else
            allocate(R(size(X)),DX(size(X)))
        endif

        LOOP: do while (.true.)
         

            !---------------------------------------------------------------------------------------------------------------
            ! Evaluating Residual - Nonlinear System of Equations
            !---------------------------------------------------------------------------------------------------------------
            call SOE%EvaluateSystem(X,R)

            if (SOE%Status%Error) then
                call this%Status%SetError(NewtonRaphsonFull_Errors%UserEvaluateSystemReportedError,'Error Evaluating system')
                return
            endif
            !---------------------------------------------------------------------------------------------------------------


            !---------------------------------------------------------------------------------------------------------------
            ! Evaluating the Residual Gradient
            !---------------------------------------------------------------------------------------------------------------
            select case (this%MatrixType)
                case (NewtonRaphsonFull_MatrixTypes%Full)
                    call SOE%EvaluateGradient(X,R,GFull)
                case (NewtonRaphsonFull_MatrixTypes%Sparse)
                    call SOE%EvaluateGradient(X,R,GSparse)
                case default
            end select

            if (SOE%Status%error) then
                call this%Status%SetError(NewtonRaphsonFull_Errors%UserEvaluateGradientReportedError,'Error Evaluating Gradient')
                return
            endif
            !---------------------------------------------------------------------------------------------------------------


            !---------------------------------------------------------------------------------------------------------------
            ! Computing the Residual Norm
            !---------------------------------------------------------------------------------------------------------------
            select case (this%normtype)
                case (NewtonRaphsonFull_NormTypes%L2Norm)
                    normR = norm(R)
                case (NewtonRaphsonFull_NormTypes%MaximumAbsoluteValue)
                    normR = maxval( dabs(R))
                case default
                    stop "NewtonRaphsonFull_Solve :: NormType not set"
            end select

            if (this%ShowInfo) write(*,'(12x,a,i3,a,e16.9)') 'IT: ',IT ,'  NORM: ',normR
            
            if (it==0) normR0=normR
            
            if ((it>0) .AND. (normR<normR0*this%tol)) then
                call this%Status%SetSuccess()
                if (this%ShowInfo) write(*,'(12x,a,i3,a)')'Converged in ',IT,' iterations'
                return
            elseif (it>= this%itmax) then
                call this%Status%SetError(NewtonRaphsonFull_Errors%MaxNumberOfIteration,'Maximum Number of Iterations reached!')
                return
            endif
            !---------------------------------------------------------------------------------------------------------------

            !---------------------------------------------------------------------------------------------------------------
            ! Update Iterations
            !---------------------------------------------------------------------------------------------------------------
            it=it+1
            SOE%it = it
            this%NumberOfIterations = it
            !---------------------------------------------------------------------------------------------------------------

            !---------------------------------------------------------------------------------------------------------------
            ! Solving the Linear System of Equations
            !---------------------------------------------------------------------------------------------------------------
            select case (this%MatrixType)
                case (NewtonRaphsonFull_MatrixTypes%Full)
                    call this%LinearSolver%Solve(GFull, -R, DX)
                case (NewtonRaphsonFull_MatrixTypes%Sparse)
                    call this%LinearSolver%Solve(GSparse, -R, DX)
                case default
                end select

            !---------------------------------------------------------------------------------------------------------------
            ! Update Unknown Variable and Additional Variables
            !---------------------------------------------------------------------------------------------------------------
                            
            call LineSearch(SOE, R, GSparse, DX, X, it)
                
            !if (SOE%isPeriodic) then
            !    call SOE%ExpandResult(DX,DXFull,it)
            !    X = X + DXFull
            !else
            !    X = X + DX
            !endif
            !call SOE%PostUpdate(X)

            !---------------------------------------------------------------------------------------------------------------

        end do LOOP

    end subroutine



    subroutine NewtonRaphsonFull_ReadSolverParameters(this,DataFile)
            use Parser
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNewtonRaphsonFull) :: this


            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser)::DataFile

		    !************************************************************************************

		    character(len=100),dimension(2)::ListOfOptions,ListOfValues

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
		    ListOfOptions=["tol","maxiter"]

		    call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)


            this%tol = ListOfValues(1)
            this%itmax = ListOfValues(2)

    end subroutine

    
    
        subroutine LineSearch(SOE, R, GSparse, DX, X, it)
            !************************************************************************************           
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNonLinearSystemOfEquations)      :: SOE
            ! Input variables
            ! ---------------------------------------------------------------------------------
            real(8),dimension(:)                        :: R , DX
            class(ClassGlobalSparseMatrix),pointer      :: GSparse
            integer                                     :: it

            ! Output variables
            real(8),dimension(:)                        :: X 
            
            !For Line Search
            real(8) :: R_scalar_0, R_scalar_eta, eta, eta_old, rho_LS, criteria_LS, alpha
            real(8) :: R_new_LS(size(DX)), X_new_LS(size(X)), DXFull(size(X))
            integer :: count_LS
            logical :: Divergence_LS
           
            !---------------------------------------------------------------------------------------------------------------
            ! Line-search / Bonet and Wood's Text Book  (2008)            
            R_scalar_0 = - dot_product(DX, R)
                
            eta = 1.0
            eta_old = eta
            rho_LS = 0.5
            Divergence_LS = .true.
            count_LS = 0
          
            if (SOE%isPeriodic) call SOE%ExpandResult(DX,DXFull,it)
            
            LS: do while (Divergence_LS)
                    
                count_LS = count_LS + 1
                write(*,'(12x,a,i3, a,e16.9)') 'Line Search Iteration: ',count_LS ,'  Step length : ',eta
                
                if (SOE%isPeriodic) then
                    X_new_LS = X + eta*DXFull
                else
                    X_new_LS = X + eta*DX
                endif
                                
                call SOE%PostUpdate(X_new_LS)
                call SOE%EvaluateSystem(X_new_LS,R_new_LS)
                call SOE%EvaluateGradient(X_new_LS,R_new_LS,GSparse)
                               
                R_scalar_eta = dot_product(DX, R_new_LS)
                    
                criteria_LS = abs(R_scalar_eta)/abs(R_scalar_0)
                    
                if (criteria_LS.lt.rho_LS) then
                    Divergence_LS = .false.
                    X = X_new_LS
                endif
                
                alpha = R_scalar_0/R_scalar_eta
                    
                if (alpha.lt.0.0) then
                    eta = (alpha/2) + sqrt((alpha/2)**2 - alpha)
                else
                    eta = alpha/2
                endif
                    
                if (eta.gt.1.0) then
                    eta = 0.99
                elseif(eta.lt.1.0e-3) then
                    eta = 1.0e-2                        
                end if
                    
                if ( (eta_old-eta).lt.1e-12) then
                    Divergence_LS = .false.
                    X = X_new_LS        
                end if
                
                eta_old = eta
                    
            enddo LS
            
        endsubroutine
        !==========================================================================================  
    
    
    
    
    
    end module

