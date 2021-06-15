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
        real(8) :: normR , norma, eta
        real(8),allocatable,dimension(:) :: R , DX, DXFull

        real(8),dimension(:,:),pointer :: GFull
        class(ClassGlobalSparseMatrix),pointer :: GSparse


        call SOE%Status%SetSuccess
        call this%Status%SetSuccess

        it = 0
        X=Xguess

        if (SOE%isPeriodic) then
            allocate(R(SOE%nDOF),DX(SOE%nDOF),DXFull(size(X))) !SOE%nDOF = DOF reduced system
        else
            allocate(R(size(X)),DX(size(X)))
        endif

        LOOP: do while (.true.)
            
            SOE%it = it

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

            if (normR<this%tol) then
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

            if (it>15) then
                eta = (sqrt(5.0D0)-1)/2
            else
                eta = 1.0D0
            endif
                
            !---------------------------------------------------------------------------------------------------------------
            ! Update Unknown Variable and Additional Variables
            !---------------------------------------------------------------------------------------------------------------
            if (SOE%isPeriodic) then
                call SOE%ExpandResult(DX,DXFull)
                X = X + eta*DXFull
            else
                X = X + eta*DX
            endif

            call SOE%PostUpdate(X)
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
        !==========================================================================================
!        subroutine NewtonRaphsonFull_Constructor (this)
!		    !************************************************************************************
!            ! DECLARATIONS OF VARIABLES
!		    !************************************************************************************
!            ! Object
!            ! ---------------------------------------------------------------------------------
!            class(ClassNewtonRaphsonFull) :: this
!		    !************************************************************************************
!
!        endsubroutine

    end module

