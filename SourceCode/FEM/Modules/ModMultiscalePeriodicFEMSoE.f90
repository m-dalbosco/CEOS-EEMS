module ModMultiscalePeriodicFEMSoE

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use BoundaryConditions
    use ElementLibrary
    use GlobalSparseMatrix
    use SparseMatrixRoutines

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
        integer                              :: nDOF
        
            call this%AnalysisSettings%GetTotalNumberOfDOF(this%GlobalNodesList,nDOF)
            
            allocate(XFull(nDOF),RFull(nDOF))
            
            XFull = 0.0d0
            call mkl_dcsrgemv('N', nDOF, this%TMat%Val, this%TMat%RowMap, this%TMat%Col, X, XFull) !Calculate XFull from X (red)
            
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
            call mkl_dcsrgemv('T', nDOF, this%TMat%Val, this%TMat%RowMap, this%TMat%Col, RFull, R)

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
        integer :: nDOF, info, nzmax
             
        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)
        
        allocate(XFull(nDOF),RFull(nDOF))

        allocate(KgAux%RowMap(size(this%Kg%RowMap)))
        allocate(KgAux%Val(size(this%Kg%Val)))
        allocate(KgAux%Col(size(this%Kg%Col)))
        KgAux%RowMap = 0.0d0
        KgAux%Val = 0.0d0
        KgAux%Col = 0.0d0
                    
        XFull = 0.0d0
        RFull = 0.0d0
        call mkl_dcsrgemv('N', nDOF, this%TMat%Val, this%TMat%RowMap, this%TMat%Col, X, XFull) !Calculate XFull from X (red)
        call mkl_dcsrgemv('N', nDOF, this%TMat%Val, this%TMat%RowMap, this%TMat%Col, R, RFull) !Calculate RFull from R (red)

        call TangentStiffnessMatrix(this%AnalysisSettings , this%ElementList , nDOF, this%Kg) !Calculate full stiffness matrix

        ! As CC de deslocamento prescrito estão sendo aplicadas no sistema Kx=R e não em Kx=-R
        RFull = -RFull
        call this%BC%ApplyBoundaryConditionsNEW( this%Kg , RFull , this%DispDOF, this%Ubar , XFull, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        RFull = -RFull
        
        nzmax = nDOF*nDOF
        call mkl_dcsrmultcsr('T', 0, 0, nDOF, nDOF, nDOF, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%Kg%Val, this%Kg%Col, this%Kg%RowMap, KgAux%Val, KgAux%Col, KgAux%RowMap, nzmax, info)
        call mkl_dcsrmultcsr('N', 0, 0, nDOF, nDOF, nDOF, KgAux%Val, KgAux%Col, KgAux%RowMap, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%KgRed%Val, this%KgRed%Col, this%KgRed%RowMap, nzmax, info)
        
        G => this%KgRed
        
        deallocate(KgAux%RowMap,KgAux%Val,KgAux%Col)
        
        R = 0.0d0
        ! Calculate R (red) from RFull
        call mkl_dcsrgemv('T', nDOF, this%TMat%Val, this%TMat%RowMap, this%TMat%Col, RFull, R)

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
        call mkl_dcsrgemv('N', nDOF, this%TMat%Val, this%TMat%RowMap, this%TMat%Col, X, XFull) !Calculate XFull from X (red)

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,XFull)
        endif

    end subroutine
    
!--------------------------------------------------------------------------------------------------

    subroutine BuildT(this,nDOFRed)
        class(ClassMultiscalePeriodicFEMSoE)   :: this
        type(SparseMatrix)                     :: TMatSparse
        integer                                :: idx, i, j, k, nDOF, nDOFRed, nNod, nNodBound
        integer                                :: nNodBX, nNodBY, nNodBZ, nNodBXY, nNodBXZ, nNodBYZ
        integer,allocatable,dimension(:)       :: Xm, Xp, Ym, Yp, Zm, Zp, XmYm, XmZm, XmYp, XmZp, XpYm, XpZm, XpYp, XpZp, YmZm, YmZp, YpZm, YpZp
        integer                                :: XmYmZm, XpYmZm, XmYpZm, XmYmZp, XpYpZp, XmYpZp, XpYmZp, XpYpZm
        integer                                :: countXm, countXp, countYm, countYp, countZm, countZp
        integer                                :: countXmYm, countXmZm, countXmYp, countXmZp, countXpYm, countXpZm, countXpYp, countXpZp, countYmZm, countYmZp, countYpZm, countYpZp
        real(8)                                :: Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
        
        
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
        
        call SparseMatrixInit(TMatSparse , nDOF)
            
        do i=1,nDOF
            do j=1,nDOF
                if (i==j) then
                    call SparseMatrixSetVal( i , j , 1.0d0 , TMatSparse )
                endif
            enddo
        enddo
        
        nDOFRed = nDOF
        
        !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
        call ConvertToCoordinateFormat( TMatSparse , this%TMat%Row , this%TMat%Col , this%TMat%Val , this%TMat%RowMap)

        !Releasing memory
        call SparseMatrixKill(TMatSparse)

    end subroutine





end module

