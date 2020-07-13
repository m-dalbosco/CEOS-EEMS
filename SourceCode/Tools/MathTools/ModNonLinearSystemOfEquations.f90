module modNonLinearSystemOfEquations

    use modStatus
    use GlobalSparseMatrix

    type ClassNonLinearSystemOfEquations

        type(ClassStatus) :: Status
        integer           :: it, nDOF
        logical           :: isPeriodic

    contains

        procedure :: EvaluateSystem => EvaluateSystemBase
        procedure :: EvaluateGradientFull => EvaluateGradientBase
        procedure :: EvaluateGradientSparse => EvaluateGradientSparseBase
        generic   :: EvaluateGradient => EvaluateGradientFull , EvaluateGradientSparse
        procedure :: PostUpdate => PostUpdateBase
        procedure :: ReduceSystem => ReduceSystemBase
        procedure :: ExpandResult => ExpandResultBase

    end type

    contains

!__________________________________________________________________________________________________
    subroutine EvaluateSystemBase(this,x,R)
        class(ClassNonLinearSystemOfEquations)::this
        real(8),dimension(:)::x,R
        stop "EvaluateSystem Not Implemented"
    end subroutine
!__________________________________________________________________________________________________
    subroutine EvaluateGradientBase(this,x,R,G)
        class(ClassNonLinearSystemOfEquations)::this
        real(8),dimension(:)::x,R
        real(8),dimension(:,:),pointer::G
        stop "EvaluateGradient Not Implemented"
    end subroutine
!__________________________________________________________________________________________________
    subroutine EvaluateGradientSparseBase(this,x,R,G)
        class(ClassNonLinearSystemOfEquations)::this
        class(ClassGlobalSparseMatrix) , pointer :: G
        real(8),dimension(:)::x,R
        stop "EvaluateGradient Not Implemented"
    end subroutine
!__________________________________________________________________________________________________
    subroutine PostUpdateBase(this,X)
        class(ClassNonLinearSystemOfEquations)::this
        real(8),dimension(:)::X
    end subroutine
!__________________________________________________________________________________________________
    subroutine ReduceSystemBase(this,R,G) !For systems with periodicity
        class(ClassNonLinearSystemOfEquations)::this
        class(ClassGlobalSparseMatrix) , pointer :: G
        real(8),dimension(:)::R
        stop "ReduceSystem Not Implemented"
    end subroutine
!__________________________________________________________________________________________________
    subroutine ExpandResultBase(this,Rred,Rfull) !For systems with periodicity
        class(ClassNonLinearSystemOfEquations)::this
        real(8),dimension(:)::Rred,Rfull
        stop "ExpandResult Not Implemented"
    end subroutine    
    
end module
