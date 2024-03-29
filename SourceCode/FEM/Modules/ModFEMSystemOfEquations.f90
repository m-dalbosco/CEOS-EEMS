module modFEMSystemOfEquations

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use BoundaryConditions
    use ElementLibrary
    use GlobalSparseMatrix

    implicit none

    type , extends(ClassNonLinearSystemOfEquations) :: ClassFEMSystemOfEquations

        real(8), dimension(:), allocatable                   :: Fint , Fext , UBar
        real(8)                                              :: Time
        integer, dimension(:), pointer                       :: DispDOF

        integer, dimension(:), allocatable                   :: PrescDispSparseMapZERO
        integer, dimension(:), allocatable                   :: PrescDispSparseMapONE
        integer, dimension(:), allocatable                   :: FixedSupportSparseMapZERO
        integer, dimension(:), allocatable                   :: FixedSupportSparseMapONE
        
        type  (ClassElementsWrapper), dimension(:), pointer  :: ElementList
        type  (ClassNodes), dimension(:), pointer            :: GlobalNodesList
        type  (ClassAnalysis)                                :: AnalysisSettings
        class (ClassBoundaryConditions), pointer             :: BC
        type  (ClassGlobalSparseMatrix), pointer             :: Kg


    contains

        procedure :: EvaluateSystem => EvaluateR
        procedure :: EvaluateGradientSparse => EvaluateKt
        procedure :: PostUpdate => FEMUpdateMesh

    end type

    contains

!--------------------------------------------------------------------------------------------------
    subroutine EvaluateR(this,X,R)

        use Interfaces
        class(ClassFEMSystemOfEquations) :: this
        real(8),dimension(:) :: X,R

            ! Update stress and internal variables
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, X, this%Status)

            ! Constitutive Model Failed. Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Internal Force
            call InternalForce(this%ElementList , this%AnalysisSettings , this%Fint, this%Status)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Residual
            R = this%Fint - this%Fext

    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine EvaluateKt(this,X,R,G,flagG)

        use Interfaces
        use MathRoutines
        class(ClassFEMSystemOfEquations)        :: this
        class (ClassGlobalSparseMatrix), pointer :: G
        real(8),dimension(:) :: X , R
        real(8) :: norma
        integer :: nDOF
        logical :: flagG

        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

        call TangentStiffnessMatrix(this%AnalysisSettings , this%ElementList , nDOF, this%Kg )

        this%BC%it = this%it
        ! As CC de deslocamento prescrito est�o sendo aplicadas no sistema Kx=R e n�o em Kx=-R!!!
        R = -R
        !call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X   )
        call this%BC%ApplyBoundaryConditionsNEW(  this%Kg , R , this%DispDOF, this%Ubar , X, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        R = -R

        G => this%Kg

    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine FEMUpdateMesh(this,X)
        use Interfaces
        class(ClassFEMSystemOfEquations) :: this
        real(8),dimension(:)::X

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,X)
        endif

    end subroutine





end module

