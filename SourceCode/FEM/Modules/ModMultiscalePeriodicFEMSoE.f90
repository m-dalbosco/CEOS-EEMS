module ModMultiscalePeriodicFEMSoE

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use BoundaryConditions
    use ElementLibrary
    use GlobalSparseMatrix
    use SparseMatrixRoutines
    use MathRoutines

    implicit none

    type , extends(ClassNonLinearSystemOfEquations) :: ClassMultiscalePeriodicFEMSoE

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
        type  (ClassGlobalSparseMatrix), pointer             :: KgRed
        type  (ClassGlobalSparseMatrix), pointer             :: TMat
        character(len=1),dimension(6)                        :: TMatDescr
        integer                                              :: nDOFRed


    contains

        procedure :: EvaluateSystem => EvaluateR
        procedure :: EvaluateGradientSparse => EvaluateKt
        procedure :: PostUpdate => FEMUpdateMesh
        procedure :: BuildT

    end type

    contains

!--------------------------------------------------------------------------------------------------
    subroutine EvaluateR(this,X,R)

        use Interfaces
        class(ClassMultiscalePeriodicFEMSoE) :: this
        real(8),dimension(:)                 :: X,R          !Reduced system
        real(8),allocatable,dimension(:)     :: XFull,RFull  !Full system
        integer                              :: nDOF, nDOFRed
        
            call this%AnalysisSettings%GetTotalNumberOfDOF(this%GlobalNodesList,nDOF)
            
            nDOFRed = this%nDOFRed
            
            allocate(XFull(nDOF),RFull(nDOF))
            
            XFull = 0.0d0
            call mkl_dcsrmv('N', nDOF, nDOFRed, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap(1:(size(this%TMat%RowMap)-1)), this%TMat%RowMap(2:size(this%TMat%RowMap)), X, 0.0d0, XFull) !Calculate XFull from X (red)
            !XFull = XFull + this%Ubar
            
            ! Update stress and internal variables
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, XFull, this%Status) !Solve constitutive model (full system)

            ! Constitutive Model Failed. Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Internal Force
            call InternalForce(this%ElementList , this%AnalysisSettings , this%Fint, this%Status) !Calculate internal force (full system)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Residual (full system)
            RFull = this%Fint - this%Fext
            
            R = 0.0d0
            ! Calculate R (red) from RFull
            call mkl_dcsrmv('T', nDOF, nDOFRed, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap(1:(size(this%TMat%RowMap)-1)), this%TMat%RowMap(2:size(this%TMat%RowMap)), RFull, 0.0d0, R)
            
    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine EvaluateKt(this,X,R,G)

        use Interfaces
        use MathRoutines
        class(ClassMultiscalePeriodicFEMSoE)       :: this
        class(ClassGlobalSparseMatrix),pointer     :: G
        type(ClassGlobalSparseMatrix)              :: KgAux
        real(8),dimension(:)                       :: X,R          !Reduced system
        real(8),allocatable,dimension(:)           :: XFull,RFull  !Full system
        real(8) :: norma
        integer :: nDOF, nDOFRed, info, nzmax, ValDum, ColDum, LengthKgAuxVal, LengthKgRedVal
                     
        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)
        
        allocate(XFull(nDOF),RFull(nDOF))
        XFull = 0.0d0
        RFull = 0.0d0
        
        nDOFRed = this%nDOFRed
        call mkl_dcsrmv('N', nDOF, nDOFRed, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap(1:(size(this%TMat%RowMap)-1)), this%TMat%RowMap(2:size(this%TMat%RowMap)), X, 0.0d0, XFull) !Calculate XFull from X (red)
        call mkl_dcsrmv('N', nDOF, nDOFRed, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap(1:(size(this%TMat%RowMap)-1)), this%TMat%RowMap(2:size(this%TMat%RowMap)), R, 0.0d0, RFull) !Calculate RFull from R (red)
        !XFull = XFull + this%Ubar
        
        call TangentStiffnessMatrix(this%AnalysisSettings , this%ElementList , nDOF, this%Kg) !Calculate full stiffness matrix
        !Print for checking
        call OutputSparseMatrix(this%Kg,'Kg.txt',nDOF,nDOF)        

        ! As CC de deslocamento prescrito est�o sendo aplicadas no sistema Kx=R e n�o em Kx=-R
        RFull = -RFull
        call this%BC%ApplyBoundaryConditionsNEW( this%Kg , RFull , this%DispDOF, this%Ubar , XFull, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        RFull = -RFull
        
        
        allocate(KgAux%RowMap(nDOFRed+1))
        if (associated(this%KgRed%RowMap)) then
            deallocate(this%KgRed%RowMap)
        endif
        if (associated(this%KgRed%Val)) then
            deallocate(this%KgRed%Val)
        endif
        if (associated(this%KgRed%Col)) then
            deallocate(this%KgRed%Col)
        endif
               
        nzmax = nDOF*nDOFRed
        !Calculate length of KgAux = Tmat'*Kg
        call mkl_dcsrmultcsr('T', 1, 0, nDOF, nDOFRed, nDOF, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%Kg%Val, this%Kg%Col, this%Kg%RowMap, ValDum, ColDum, KgAux%RowMap, nzmax, info)
        LengthKgAuxVal = KgAux%RowMap(nDOFRed+1)-1
        allocate(KgAux%Val(LengthKgAuxVal))
        allocate(KgAux%Col(LengthKgAuxVal))
        KgAux%Val = 0.0d0
        KgAux%Col = 0
        !Calculate KgAux = Tmat'*Kg
        call mkl_dcsrmultcsr('T', 2, 0, nDOF, nDOFRed, nDOF, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%Kg%Val, this%Kg%Col, this%Kg%RowMap, KgAux%Val, KgAux%Col, KgAux%RowMap, nzmax, info)
        !Print for checking
        call OutputSparseMatrix(KgAux,'KgAux.txt',nDOFRed,nDOF)
        
        allocate(this%KgRed%RowMap(nDOFRed+1))
        !Calculate length of KgRed = KgAux*TMat
        call mkl_dcsrmultcsr('N', 1, 0, nDOFRed, nDOF, nDOFRed, KgAux%Val, KgAux%Col, KgAux%RowMap, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, ValDum, ColDum, this%KgRed%RowMap, nzmax, info)
        LengthKgRedVal = this%KgRed%RowMap(nDOFRed+1)-1
        allocate(this%KgRed%Val(LengthKgRedVal))
        allocate(this%KgRed%Col(LengthKgRedVal))
        this%KgRed%Val = 0.0d0
        this%KgRed%Col = 0
        !Calculate KgRed = KgAux*TMat
        call mkl_dcsrmultcsr('N', 2, 0, nDOFRed, nDOF, nDOFRed, KgAux%Val, KgAux%Col, KgAux%RowMap, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%KgRed%Val, this%KgRed%Col, this%KgRed%RowMap, nzmax, info)
        !Print for checking
        call OutputSparseMatrix(this%KgRed,'KgRed.txt',nDOFRed,nDOFRed)
        
        G => this%KgRed
        
        deallocate(KgAux%RowMap,KgAux%Val,KgAux%Col)
        
        R = 0.0d0
        ! Calculate R (red) from RFull
        call mkl_dcsrmv('T', nDOF, nDOFRed, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap(1:(size(this%TMat%RowMap)-1)), this%TMat%RowMap(2:size(this%TMat%RowMap)), RFull, 0.0d0, R)

    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine FEMUpdateMesh(this,X)
        use Interfaces
        class(ClassMultiscalePeriodicFEMSoE) :: this
        real(8),dimension(:)             :: X
        real(8),allocatable,dimension(:) :: XFull
        integer                          :: nDOF
        
        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)
        
        allocate(XFull(nDOF))
        
        XFull = 0.0d0
        call mkl_dcsrmv('N', nDOF, nDOF, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%TMat%RowMap(2), X, 0.0d0, XFull) !Calculate XFull from X (red)

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,XFull)
        endif

    end subroutine
    
!--------------------------------------------------------------------------------------------------

    subroutine BuildT(this)
        class(ClassMultiscalePeriodicFEMSoE) :: this
        type(SparseMatrix)                   :: TMatSparse
        integer                              :: idx, i, j, k, m, n, col, nDOF, nDOFRed, nNod, nNodBound, FileID_TMatFull, FileID_TMatRed
        integer                              :: nNodBX, nNodBY, nNodBZ, nNodBXY, nNodBXZ, nNodBYZ
        integer,allocatable,dimension(:)     :: Xm, Xp, Ym, Yp, Zm, Zp, XmYm, XmZm, XmYp, XmZp, XpYm, XpZm, XpYp, XpZp, YmZm, YmZp, YpZm, YpZp
        integer                              :: XmYmZm, XpYmZm, XmYpZm, XmYmZp, XpYpZp, XmYpZp, XpYmZp, XpYpZm
        integer                              :: countXm, countXp, countYm, countYp, countZm, countZp
        integer                              :: countXmYm, countXmZm, countXmYp, countXmZp, countXpYm, countXpZm, countXpYp, countXpZp, countYmZm, countYmZp, countYpZm, countYpZp
        real(8)                              :: Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
        integer,allocatable,dimension(:,:)   :: TMatFull, TMatRed
        
        
        nNod = size(this%GlobalNodesList)
        nNodBound = size(this%BC%BoundaryNodes(1)%Nodes)
        
        !Finding RVE borders
        !------------------------------------------
        Xmin=0.0d0
        Xmax=0.0d0
        Ymin=0.0d0
        Ymax=0.0d0
        Zmin=0.0d0
        Zmax=0.0d0
        do k=1,nNodBound
            idx = this%BC%BoundaryNodes(1)%Nodes(k)
            if (this%GlobalNodesList(idx)%CoordX(1) < Xmin) then
                Xmin=this%GlobalNodesList(idx)%CoordX(1)
            endif
            if (this%GlobalNodesList(idx)%CoordX(1) > Xmax) then
                Xmax=this%GlobalNodesList(idx)%CoordX(1)
            endif
            if (this%GlobalNodesList(idx)%CoordX(2) < Ymin) then
                Ymin=this%GlobalNodesList(idx)%CoordX(2)
            endif
            if (this%GlobalNodesList(idx)%CoordX(2) > Ymax) then
                Ymax=this%GlobalNodesList(idx)%CoordX(2)
            endif
            if (this%GlobalNodesList(idx)%CoordX(3) < Zmin) then
                Zmin=this%GlobalNodesList(idx)%CoordX(3)
            endif
            if (this%GlobalNodesList(idx)%CoordX(3) > Zmax) then
                Zmax=this%GlobalNodesList(idx)%CoordX(3)
            endif
        enddo
        
        
        !Finding nodes in vertices and number of nodes in edges and inside faces
        !------------------------------------------
        nNodBX=0.0d0
        nNodBY=0.0d0
        nNodBZ=0.0d0
        nNodBXY=0.0d0
        nNodBXZ=0.0d0
        nNodBYZ=0.0d0
        do k=1,nNodBound
            idx = this%BC%BoundaryNodes(1)%Nodes(k)
            if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in Xm
                if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in XmYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node XmYmZm
                        XmYmZm = idx
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node XmYmZp
                        XmYmZp = idx
                    else !node inside XmYm
                        nNodBXY = nNodBXY+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node in XmYp
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node XmYpZm
                        XmYpZm = idx
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node XmYpZp
                        XmYpZp = idx
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in XmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node XmYpZm
                        !node XmYpZm already found
                    else !node inside XmZm
                        nNodBXZ = nNodBXZ+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node in XmZp
                    !nNodBXZ already counted
                else !node inside Xm
                    nNodBX=nNodBX+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in Ym
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in XmYm
                    !nNodBXY already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node in XpYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node XpYmZm
                        XpYmZm = idx
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node XpYmZp
                        XpYmZp = idx
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in YmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node XpYmZm
                        !node XpYmZm already found
                    else !node inside YmZm
                        nNodBYZ = nNodBYZ+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node in YmZp
                    !nNodBYZ already counted
                else !node inside Ym
                    nNodBY=nNodBY+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in Zm
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in XmZm
                    !nNodBXZ already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node in XpZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node XpYmZm
                        !node XpYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node XpYpZm
                        XpYpZm = idx
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in YmZm
                    !nNodBYZ already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node in YpZm
                    !nNodBYZ already counted
                else !node inside Zm
                    nNodBZ=nNodBZ+1
                endif
            endif

            if ((abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) .AND. (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) .AND. (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12)) then
                XpYpZp = idx
            endif
                    
        enddo
        
        allocate(Xm(nNodBX), Xp(nNodBX), Ym(nNodBY), Yp(nNodBY), Zm(nNodBZ), Zp(nNodBZ), XmYm(nNodBXY), XmZm(nNodBXZ), XmYp(nNodBXY), XmZp(nNodBXZ), XpYm(nNodBXY), XpZm(nNodBXZ), XpYp(nNodBXY), XpZp(nNodBXZ), YmZm(nNodBYZ), YmZp(nNodBYZ), YpZm(nNodBYZ), YpZp(nNodBYZ))
        
        !Finding nodes in edges and inside faces
        !------------------------------------------
        countXm = 1
        countXp = 1
        countYm = 1
        countYp = 1
        countZm = 1
        countZp = 1
        countXmYm = 1
        countXmZm = 1
        countXmYp = 1
        countXmZp = 1
        countXpYm = 1
        countXpZm = 1
        countXpYp = 1
        countXpZp = 1
        countYmZm = 1 
        countYmZp = 1
        countYpZm = 1
        countYpZp = 1
        do k=1,nNodBound
            idx = this%BC%BoundaryNodes(1)%Nodes(k)
            if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in Xm
                if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in XmYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node XmYmZm
                        !XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node XmYmZp
                        !XmYmZp already found
                    else !node inside XmYm
                        XmYm(countXmYm) = idx
                        countXmYm = countXmYm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node in XmYp
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node XmYpZm
                        !XmYpZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node XmYpZp
                        !XmYpZp already found
                    else !node inside XmYp
                        XmYp(countXmYp) = idx
                        countXmYp = countXmYp+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in XmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node XmYpZm
                        !node XmYpZm already found
                    else !node inside XmZm
                        XmZm(countXmZm) = idx
                        countXmZm = countXmZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node in XmZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node XmYmZp
                        !node XmYmZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node XmYpZp
                        !node XmYpZp already found
                    else !node inside XmZp
                        XmZp(countXmZp) = idx
                        countXmZp = countXmZp+1
                    endif
                else !node inside Xm
                    Xm(countXm) = idx
                    countXm = countXm+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node in Xp
                if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in XpYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node XpYmZm
                        !XpYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node XpYmZp
                        !XpYmZp already found
                    else !node inside XpYm
                        XpYm(countXpYm) = idx
                        countXpYm = countXpYm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node in XpYp
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node XpYpZm
                        !XpYpZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node XpYpZp
                        !XpYpZp already found
                    else !node inside XpYp
                        XpYp(countXpYp) = idx
                        countXpYp = countXpYp+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in XpZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node XpYmZm
                        !node XpYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node XpYpZm
                        !node XpYpZm already found
                    else !node inside XpZm
                        XpZm(countXpZm) = idx
                        countXpZm = countXpZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node in XpZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node XpYmZp
                        !node XpYmZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node XpYpZp
                        !node XpYpZp already found
                    else !node inside XpZp
                        XpZp(countXpZp) = idx
                        countXpZp = countXpZp+1
                    endif
                else !node inside Xp
                    Xp(countXp) = idx
                    countXp = countXp+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in Ym
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in XmYm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node in XpYm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in YmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node XpYmZm
                        !node XpYmZm already found
                    else !node inside YmZm
                        YmZm(countYmZm) = idx
                        countYmZm = countYmZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node in YmZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node XmYmZp
                        !node XmYmZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node XpYmZp
                        !node XpYmZp already found
                    else !node inside YmZp
                        YmZp(countYmZp) = idx
                        countYmZp = countYmZp+1
                    endif
                else !node inside Ym
                    Ym(countYm) = idx
                    countYm = countYm+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node in Yp
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in XmYp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node in XpYp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in YpZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node XmYpZm
                        !node XmYpZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node XpYpZm
                        !node XpYpZm already found
                    else !node inside YpZm
                        YpZm(countYpZm) = idx
                        countYpZm = countYpZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node in YpZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node XmYpZp
                        !node XmYpZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node XpYpZp
                        !node XpYpZp already found
                    else !node inside YpZp
                        YpZp(countYpZp) = idx
                        countYpZp = countYpZp+1
                    endif
                else !node inside Yp
                    Yp(countYp) = idx
                    countYp = countYp+1
                endif            
            elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<1.0D-12) then !node in Zm
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in XmZm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node in XpZm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in YmZm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node in YpZm
                    !already counted
                else !node inside Zm
                    Zm(countZm) = idx
                    countZm = countZm+1
                endif
             elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<1.0D-12) then !node in Zp
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<1.0D-12) then !node in XmZp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<1.0D-12) then !node in XpZp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<1.0D-12) then !node in YmZp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<1.0D-12) then !node in YpZp
                    !already counted
                else !node inside Zp
                    Zp(countZp) = idx
                    countZp = countZp+1
                endif
            endif
                    
        enddo
                
        call this%AnalysisSettings%GetTotalNumberOfDOF(this%GlobalNodesList, nDOF)
        
        nDOFRed = nDOF - 3*(nNodBX+nNodBY+nNodBZ+(3*nNodBXY)+(3*nNodBXZ)+(3*nNodBYZ)+7)

        allocate(TMatFull(nDOF,nDOF))
        
        TMatFull = 0.0d0
        
        do i=1,nDOF
            TMatFull(i,i)=1
        enddo
                
        !Impose free DOF of vertex XmYmZm
        TMatFull(3*XmYmZm-2,3*XmYmZm-2)=1
        TMatFull(3*XmYmZm-1,3*XmYmZm-1)=1
        TMatFull(3*XmYmZm,3*XmYmZm)=1
        !Impose periodicity on the remaining vertices
        TMatFull(3*XpYmZm-2,3*XmYmZm-2)=1
        TMatFull(3*XpYmZm-1,3*XmYmZm-1)=1
        TMatFull(3*XpYmZm,3*XmYmZm)=1
        TMatFull(3*XpYmZm-2,3*XpYmZm-2)=0
        TMatFull(3*XpYmZm-1,3*XpYmZm-1)=0
        TMatFull(3*XpYmZm,3*XpYmZm)=0
        TMatFull(3*XmYpZm-2,3*XmYmZm-2)=1
        TMatFull(3*XmYpZm-1,3*XmYmZm-1)=1
        TMatFull(3*XmYpZm,3*XmYmZm)=1
        TMatFull(3*XmYpZm-2,3*XmYpZm-2)=0
        TMatFull(3*XmYpZm-1,3*XmYpZm-1)=0
        TMatFull(3*XmYpZm,3*XmYpZm)=0
        TMatFull(3*XmYmZp-2,3*XmYmZm-2)=1
        TMatFull(3*XmYmZp-1,3*XmYmZm-1)=1
        TMatFull(3*XmYmZp,3*XmYmZm)=1
        TMatFull(3*XmYmZp-2,3*XmYmZp-2)=0
        TMatFull(3*XmYmZp-1,3*XmYmZp-1)=0
        TMatFull(3*XmYmZp,3*XmYmZp)=0
        TMatFull(3*XpYpZp-2,3*XmYmZm-2)=1
        TMatFull(3*XpYpZp-1,3*XmYmZm-1)=1
        TMatFull(3*XpYpZp,3*XmYmZm)=1
        TMatFull(3*XpYpZp-2,3*XpYpZp-2)=0
        TMatFull(3*XpYpZp-1,3*XpYpZp-1)=0
        TMatFull(3*XpYpZp,3*XpYpZp)=0
        TMatFull(3*XmYpZp-2,3*XmYmZm-2)=1
        TMatFull(3*XmYpZp-1,3*XmYmZm-1)=1
        TMatFull(3*XmYpZp,3*XmYmZm)=1
        TMatFull(3*XmYpZp-2,3*XmYpZp-2)=0
        TMatFull(3*XmYpZp-1,3*XmYpZp-1)=0
        TMatFull(3*XmYpZp,3*XmYpZp)=0
        TMatFull(3*XpYmZp-2,3*XmYmZm-2)=1
        TMatFull(3*XpYmZp-1,3*XmYmZm-1)=1
        TMatFull(3*XpYmZp,3*XmYmZm)=1
        TMatFull(3*XpYmZp-2,3*XpYmZp-2)=0
        TMatFull(3*XpYmZp-1,3*XpYmZp-1)=0
        TMatFull(3*XpYmZp,3*XpYmZp)=0
        TMatFull(3*XpYpZm-2,3*XmYmZm-2)=1
        TMatFull(3*XpYpZm-1,3*XmYmZm-1)=1
        TMatFull(3*XpYpZm,3*XmYmZm)=1
        TMatFull(3*XpYpZm-2,3*XpYpZm-2)=0
        TMatFull(3*XpYpZm-1,3*XpYpZm-1)=0
        TMatFull(3*XpYpZm,3*XpYpZm)=0
        
        !Periodicity of XY edges
        do m=1,nNodBXY
            !Find periodic nodes in XmYp, XpYp, XpYm and imposes periodicity
            do n=1,nNodBXY
                if (abs(this%GlobalNodesList(XmYp(n))%CoordX(3) - this%GlobalNodesList(XmYm(m))%CoordX(3))<1.0D-12) then
                    TMatFull(3*XmYp(n)-2,3*XmYm(m)-2)=1
                    TMatFull(3*XmYp(n)-1,3*XmYm(m)-1)=1
                    TMatFull(3*XmYp(n),3*XmYm(m))=1
                    TMatFull(3*XmYp(n)-2,3*XmYp(n)-2)=0
                    TMatFull(3*XmYp(n)-1,3*XmYp(n)-1)=0
                    TMatFull(3*XmYp(n),3*XmYp(n))=0
                endif
                if (abs(this%GlobalNodesList(XpYp(n))%CoordX(3) - this%GlobalNodesList(XmYm(m))%CoordX(3))<1.0D-12) then
                    TMatFull(3*XpYp(n)-2,3*XmYm(m)-2)=1
                    TMatFull(3*XpYp(n)-1,3*XmYm(m)-1)=1
                    TMatFull(3*XpYp(n),3*XmYm(m))=1
                    TMatFull(3*XpYp(n)-2,3*XpYp(n)-2)=0
                    TMatFull(3*XpYp(n)-1,3*XpYp(n)-1)=0
                    TMatFull(3*XpYp(n),3*XpYp(n))=0
                endif
                if (abs(this%GlobalNodesList(XpYm(n))%CoordX(3) - this%GlobalNodesList(XmYm(m))%CoordX(3))<1.0D-12) then
                    TMatFull(3*XpYm(n)-2,3*XmYm(m)-2)=1
                    TMatFull(3*XpYm(n)-1,3*XmYm(m)-1)=1
                    TMatFull(3*XpYm(n),3*XmYm(m))=1
                    TMatFull(3*XpYm(n)-2,3*XpYm(n)-2)=0
                    TMatFull(3*XpYm(n)-1,3*XpYm(n)-1)=0
                    TMatFull(3*XpYm(n),3*XpYm(n))=0
                endif
            enddo
        enddo
        
        !Periodicity of XZ edges
        do m=1,nNodBXZ
            !Find periodic nodes in XmZp, XpZp, XpZm and imposes periodicity
            do n=1,nNodBXZ
                if (abs(this%GlobalNodesList(XmZp(n))%CoordX(2) - this%GlobalNodesList(XmZm(m))%CoordX(2))<1.0D-12) then
                    TMatFull(3*XmZp(n)-2,3*XmZm(m)-2)=1
                    TMatFull(3*XmZp(n)-1,3*XmZm(m)-1)=1
                    TMatFull(3*XmZp(n),3*XmZm(m))=1
                    TMatFull(3*XmZp(n)-2,3*XmZp(n)-2)=0
                    TMatFull(3*XmZp(n)-1,3*XmZp(n)-1)=0
                    TMatFull(3*XmZp(n),3*XmZp(n))=0
                endif
                if (abs(this%GlobalNodesList(XpZp(n))%CoordX(2) - this%GlobalNodesList(XmZm(m))%CoordX(2))<1.0D-12) then
                    TMatFull(3*XpZp(n)-2,3*XmZm(m)-2)=1
                    TMatFull(3*XpZp(n)-1,3*XmZm(m)-1)=1
                    TMatFull(3*XpZp(n),3*XmZm(m))=1
                    TMatFull(3*XpZp(n)-2,3*XpZp(n)-2)=0
                    TMatFull(3*XpZp(n)-1,3*XpZp(n)-1)=0
                    TMatFull(3*XpZp(n),3*XpZp(n))=0
                endif
                if (abs(this%GlobalNodesList(XpZm(n))%CoordX(2) - this%GlobalNodesList(XmZm(m))%CoordX(2))<1.0D-12) then
                    TMatFull(3*XpZm(n)-2,3*XmZm(m)-2)=1
                    TMatFull(3*XpZm(n)-1,3*XmZm(m)-1)=1
                    TMatFull(3*XpZm(n),3*XmZm(m))=1
                    TMatFull(3*XpZm(n)-2,3*XpZm(n)-2)=0
                    TMatFull(3*XpZm(n)-1,3*XpZm(n)-1)=0
                    TMatFull(3*XpZm(n),3*XpZm(n))=0
                endif
            enddo
        enddo
        
        !Periodicity of YZ edges
        do m=1,nNodBYZ
            !Find periodic nodes in XmZp, XpZp, XpZm and imposes periodicity
            do n=1,nNodBYZ
                if (abs(this%GlobalNodesList(YmZp(n))%CoordX(1) - this%GlobalNodesList(YmZm(m))%CoordX(1))<1.0D-12) then
                    TMatFull(3*YmZp(n)-2,3*YmZm(m)-2)=1
                    TMatFull(3*YmZp(n)-1,3*YmZm(m)-1)=1
                    TMatFull(3*YmZp(n),3*YmZm(m))=1
                    TMatFull(3*YmZp(n)-2,3*YmZp(n)-2)=0
                    TMatFull(3*YmZp(n)-1,3*YmZp(n)-1)=0
                    TMatFull(3*YmZp(n),3*YmZp(n))=0
                endif
                if (abs(this%GlobalNodesList(YpZp(n))%CoordX(1) - this%GlobalNodesList(YmZm(m))%CoordX(1))<1.0D-12) then
                    TMatFull(3*YpZp(n)-2,3*YmZm(m)-2)=1
                    TMatFull(3*YpZp(n)-1,3*YmZm(m)-1)=1
                    TMatFull(3*YpZp(n),3*YmZm(m))=1
                    TMatFull(3*YpZp(n)-2,3*YpZp(n)-2)=0
                    TMatFull(3*YpZp(n)-1,3*YpZp(n)-1)=0
                    TMatFull(3*YpZp(n),3*YpZp(n))=0
                endif
                if (abs(this%GlobalNodesList(YpZm(n))%CoordX(1) - this%GlobalNodesList(YmZm(m))%CoordX(1))<1.0D-12) then
                    TMatFull(3*YpZm(n)-2,3*YmZm(m)-2)=1
                    TMatFull(3*YpZm(n)-1,3*YmZm(m)-1)=1
                    TMatFull(3*YpZm(n),3*YmZm(m))=1
                    TMatFull(3*YpZm(n)-2,3*YpZm(n)-2)=0
                    TMatFull(3*YpZm(n)-1,3*YpZm(n)-1)=0
                    TMatFull(3*YpZm(n),3*YpZm(n))=0
                endif
            enddo
        enddo
        
        !Periodicity of X face
        do m=1,nNodBX
            !Find periodic nodes in Xp and imposes periodicity
            do n=1,nNodBX
                if ((abs(this%GlobalNodesList(Xp(n))%CoordX(2) - this%GlobalNodesList(Xm(m))%CoordX(2))<1.0D-12) .AND. (abs(this%GlobalNodesList(Xp(n))%CoordX(3) - this%GlobalNodesList(Xm(m))%CoordX(3))<1.0D-12)) then
                    TMatFull(3*Xp(n)-2,3*Xm(m)-2)=1
                    TMatFull(3*Xp(n)-1,3*Xm(m)-1)=1
                    TMatFull(3*Xp(n),3*Xm(m))=1
                    TMatFull(3*Xp(n)-2,3*Xp(n)-2)=0
                    TMatFull(3*Xp(n)-1,3*Xp(n)-1)=0
                    TMatFull(3*Xp(n),3*Xp(n))=0
                endif
            enddo
        enddo
        
        !Periodicity of Y face
        do m=1,nNodBY
            !Find periodic nodes in Yp and imposes periodicity
            do n=1,nNodBY
                if ((abs(this%GlobalNodesList(Yp(n))%CoordX(1) - this%GlobalNodesList(Ym(m))%CoordX(1))<1.0D-12) .AND. (abs(this%GlobalNodesList(Yp(n))%CoordX(3) - this%GlobalNodesList(Ym(m))%CoordX(3))<1.0D-12)) then
                    TMatFull(3*Yp(n)-2,3*Ym(m)-2)=1
                    TMatFull(3*Yp(n)-1,3*Ym(m)-1)=1
                    TMatFull(3*Yp(n),3*Ym(m))=1
                    TMatFull(3*Yp(n)-2,3*Yp(n)-2)=0
                    TMatFull(3*Yp(n)-1,3*Yp(n)-1)=0
                    TMatFull(3*Yp(n),3*Yp(n))=0
                endif
            enddo
        enddo
        
        !Periodicity of Z face
        do m=1,nNodBZ
            !Finds periodic nodes in Zp and imposes periodicity
            do n=1,nNodBZ
                if ((abs(this%GlobalNodesList(Zp(n))%CoordX(1) - this%GlobalNodesList(Zm(m))%CoordX(1))<1.0D-12) .AND. (abs(this%GlobalNodesList(Zp(n))%CoordX(2) - this%GlobalNodesList(Zm(m))%CoordX(2))<1.0D-12)) then
                    TMatFull(3*Zp(n)-2,3*Zm(m)-2)=1
                    TMatFull(3*Zp(n)-1,3*Zm(m)-1)=1
                    TMatFull(3*Zp(n),3*Zm(m))=1
                    TMatFull(3*Zp(n)-2,3*Zp(n)-2)=0
                    TMatFull(3*Zp(n)-1,3*Zp(n)-1)=0
                    TMatFull(3*Zp(n),3*Zp(n))=0
                endif
            enddo
        enddo

        allocate(TMatRed(nDOF,nDOFRed))
        col = 1
        do j=1,nDOF
            if (maxval(TMatFull(:,j))>0) then
                TMatRed(:,col) = TMatFull(:,j)
                col = col+1
            endif
        enddo
        
        !Output matrix and release memory
        call OutputIntMatrix(TMatFull,'TMatFull.txt')
        deallocate(TMatFull)
        
        call SparseMatrixInit(TMatSparse , nDOF)
        
        !do i=1,nDOF
        !    do j=1,nDOFRed
        !        if (TMatRed(i,j)>0) then
        !            call SparseMatrixSetVal( i , j , 1.0d0 , TMatSparse )
        !        endif
        !    enddo
        !enddo
        
        !Identity
        do i=1,nDOF
            do j=1,nDOF
                if (i==j) then
                    call SparseMatrixSetVal( i , j , 1.0d0 , TMatSparse )
                endif
            enddo
        enddo
        nDOFRed = nDOF
        
        !Output matrix and release memory
        call OutputIntMatrix(TMatRed,'TMatRed.txt')
        deallocate(TMatRed)

        !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
        call ConvertToCoordinateFormat( TMatSparse , this%TMat%Row , this%TMat%Col , this%TMat%Val , this%TMat%RowMap)

        !Output matrix and release memory
        call OutputSparseMatrix(this%TMat,'TMatRedSparse.txt',nDOF,nDOFRed)
        call SparseMatrixKill(TMatSparse)
        
        this%nDOFRed = nDOFRed
        
    end subroutine


    !---------------------------------------------------------------------------------------------------------------------
    
    subroutine OutputIntMatrix(Mat,filename)
        integer,allocatable,dimension(:,:) :: Mat        
        character(len=*)                   :: filename
        character(len=30)                  :: Form
        integer                            :: nLin, nCol, FileID, i, j

        nLin = size(Mat,1)
        nCol = size(Mat,2)
        
        if (nCol == 24) then
            Form = '24(1X,I1)'
        elseif (nCol == 81) then
            Form = '81(1X,I1)'
        elseif (nCol == 192) then
            Form = '192(1X,I1)'
        elseif (nCol == 375) then
            Form = '375(1X,I1)'
        endif
                
        FileID = 73
        open(FileID,file=filename,status='unknown')
        do i=1,nLin
            write(FileID,'('//trim(Form)//')') (Mat(i,j), j=1,nCol)
        enddo
        close(FileID)
        
    end subroutine
    
    !---------------------------------------------------------------------------------------------------------------------
    
    subroutine OutputRealMatrix(Mat,filename)
        real(8),allocatable,dimension(:,:) :: Mat        
        character(len=*)                   :: filename
        character(len=30)                  :: Form
        integer                            :: nLin, nCol, FileID, i, j

        nLin = size(Mat,1)
        nCol = size(Mat,2)
        
        if (nCol == 24) then
            Form = '24(1X,E16.9)'
        elseif (nCol == 81) then
            Form = '81(1X,E16.9)'
        elseif (nCol == 192) then
            Form = '192(1X,E16.9)'
        elseif (nCol == 375) then
            Form = '375(1X,E16.9)'
        endif
                
        FileID = 37
        open(FileID,file=filename,status='unknown')
        do i=1,nLin
            write(FileID,'('//trim(Form)//')') (Mat(i,j), j=1,nCol)
        enddo
        close(FileID)
        
    end subroutine
    
    !---------------------------------------------------------------------------------------------------------------------
    
    subroutine OutputSparseMatrix(Mat,filename,nLin,nCol)
        type(ClassGlobalSparseMatrix)      :: Mat
        character(len=*)                   :: filename
        integer                            :: i, j, nLin, nCol, nVals, FileID
        real(8),allocatable,dimension(:,:) :: MatFull
        real(8),allocatable,dimension(:)   :: Vals
        integer,allocatable,dimension(:)   :: Cols
        
        allocate(MatFull(nLin,nCol))
        MatFull = 0.0d0
        do i=1,nLin
            nVals = size(Mat%Val(Mat%RowMap(i):(Mat%RowMap(i+1)-1)))
            allocate(Vals(nVals),Cols(nVals))
            Vals = Mat%Val(Mat%RowMap(i):(Mat%RowMap(i+1)-1))
            Cols = Mat%Col(Mat%RowMap(i):(Mat%RowMap(i+1)-1))
            do j=1,nVals
                MatFull(i,Cols(j)) = Vals(j)
            enddo
            deallocate(Vals,Cols)
        enddo
        
        call OutputRealMatrix(MatFull,filename)
        deallocate(MatFull)
        
    end subroutine
    
    

end module

