!##################################################################################################
! This module has a FEM Analysis
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module FEMAnalysis

	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ElementLibrary
    use Nodes
    use ModAnalysis
    use BoundaryConditions
    use GlobalSparseMatrix
    use NonLinearSolver

    implicit none




	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassFEMAnalysis: Definitions of FEM analysis
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassFEMAnalysis

		! Class Attributes
		!----------------------------------------------------------------------------------------
        type  (ClassElementsWrapper)    , pointer     , dimension(:) :: ElementList
        type  (ClassNodes)              , pointer     , dimension(:) :: GlobalNodesList
        type  (ClassAnalysis)           , pointer                    :: AnalysisSettings
        class (ClassBoundaryConditions) , pointer                    :: BC
        type  (ClassGlobalSparseMatrix) , pointer                    :: Kg


        class (ClassNonLinearSolver)    , pointer                    :: NLSolver

        ! Para usar no Probe...
        real(8), pointer, dimension(:) :: U => null()
        real (8) :: Time
        integer :: LoadCase

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
            procedure :: ReadInputData
            procedure :: AllocateKgSparse
            procedure :: AllocateKgSparseUpperTriangular
            procedure :: AllocateKgSparseMultiscaleMinimal
            procedure :: AllocateKgSparseMultiscaleMinimalUpperTriangular
            procedure :: AllocateKgSparseMultiscaleMinimalPhases
            procedure :: AllocateKgSparseMultiscaleMinimalPhasesUpperTriangular
            procedure :: AllocateKgSparseMultiscaleMinimalLinearD1
            procedure :: AllocateKgSparseMultiscaleMinimalLinearD1UpperTriangular
            procedure :: AllocateKgSparseMultiscaleMinimalLinearD3
            procedure :: AllocateKgSparseMultiscaleMinimalLinearD3UpperTriangular            
            procedure :: Solve => SolveFEMAnalysis
            procedure :: AdditionalMaterialModelRoutine

    end type

    contains


        !==========================================================================================
        ! Method ClassFEMAnalysis:
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  ReadInputData(this,FileName)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use modReadInputFile , only : readinputfile

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            character (len=*) :: FileName
		    !************************************************************************************

 		    !************************************************************************************
            ! SELECT PARAMETERS OF THE analysis type
		    !************************************************************************************

            allocate(this%BC)
            allocate(this%Kg)
            ! Reading the input files
            !************************************************************************************
            call ReadInputFile( FileName, this%AnalysisSettings , this%GlobalNodesList , this%ElementList , &
                                this%BC , this%NLSolver )

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparse (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Ke
            integer , pointer , dimension(:)    :: GM
            integer ::  e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element stiffness matrix used to allocate memory (module Analysis)
            Ke_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , nDOF )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)

                Ke => Ke_Memory( 1:nDOFel , 1:nDOFel )
                GM => GM_Memory( 1:nDOFel )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                call SparseMatrixSetArray( GM, GM, Ke, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################

              !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparseUpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Ke
            integer , pointer , dimension(:)    :: GM
            integer ::  e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element stiffness matrix used to allocate memory (module Analysis)
            Ke_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , nDOF )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)

                Ke => Ke_Memory( 1:nDOFel , 1:nDOFel )
                GM => GM_Memory( 1:nDOFel )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings, GM )

                call SparseMatrixSetArray( GM, GM, Ke, KgSparse, OPT_SET )
                
            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################

        


        !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimal (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,12
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################
        
                !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalUpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,12
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################
        
        

        !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalPhases (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF, NumberOfMaterials

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            NumberOfMaterials = this%AnalysisSettings%NumberOfMaterials

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+9*NumberOfMaterials+3) )


            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+9*NumberOfMaterials+3) , 1:(nDOFel+9*NumberOfMaterials+3) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+9*NumberOfMaterials+3) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,9*NumberOfMaterials+3
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################

                !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalPhasesUpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF, NumberOfMaterials

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            NumberOfMaterials = this%AnalysisSettings%NumberOfMaterials

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+9*NumberOfMaterials+3) )


            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+9*NumberOfMaterials+3) , 1:(nDOFel+9*NumberOfMaterials+3) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+9*NumberOfMaterials+3) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,9*NumberOfMaterials+3
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################
        
        !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD1 (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,12
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################
        
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD1UpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,12
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################

        !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD3 (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,12
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################
        
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD3UpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use SparseMatrixRoutines
            use ModAnalysis
            use Element
            use GlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                do i=1,12
                    GM( nDOFel + i ) = nDOF + i
                enddo
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################        
        
        
        !==========================================================================================
        ! Method ClassFEMAnalysis:
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  SolveFEMAnalysis( this )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer :: nDOF

 		    !************************************************************************************
            ! SELECT PARAMETERS OF THE analysis type
		    !************************************************************************************


            ! Calling the quasi-static analysis routine
            !************************************************************************************
            select case ( this%AnalysisSettings%AnalysisType )

                case ( AnalysisTypes%Quasi_Static )

                    if (this%AnalysisSettings%MultiscaleAnalysis) then

                        if ((this%AnalysisSettings%MultiscaleModel == MultiscaleModels%Taylor) .or. (this%AnalysisSettings%MultiscaleModel == MultiscaleModels%Linear) ) then

                            call this%AdditionalMaterialModelRoutine()
                            call QuasiStaticAnalysisFEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                         this%BC, this%Kg, this%NLSolver )

                        elseif (this%AnalysisSettings%MultiscaleModel == MultiscaleModels%Minimal) then

                            call this%AdditionalMaterialModelRoutine()
                            call QuasiStaticAnalysisMultiscaleMinimalFEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                                          this%BC, this%Kg, this%NLSolver )

                        elseif (this%AnalysisSettings%MultiscaleModel == MultiscaleModels%MinimalPhases) then

                            call this%AdditionalMaterialModelRoutine()
                            call QuasiStaticAnalysisMultiscaleMinimalPhasesFEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                                          this%BC, this%Kg, this%NLSolver )

                        elseif (this%AnalysisSettings%MultiscaleModel == MultiscaleModels%MinimalLinearD1) then

                            call this%AdditionalMaterialModelRoutine()
                            call QuasiStaticAnalysisMultiscaleMinimalLinearD1FEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                                          this%BC, this%Kg, this%NLSolver )

                        elseif (this%AnalysisSettings%MultiscaleModel == MultiscaleModels%MinimalLinearD3) then

                            call this%AdditionalMaterialModelRoutine()
                            call QuasiStaticAnalysisMultiscaleMinimalLinearD1FEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                                          this%BC, this%Kg, this%NLSolver )
                            
                        endif

                    else

                            call this%AdditionalMaterialModelRoutine()
                            call QuasiStaticAnalysisFEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                         this%BC, this%Kg, this%NLSolver )
                    endif

                case default
                    stop "Error in AnalysisType - ModFEMAnalysis"
            end select



		    !************************************************************************************

        end subroutine
        !==========================================================================================



        !==========================================================================================
        ! Method ClassFEMAnalysis:
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  WriteFEMResults( U, Time, LC, ST, CutBack, SubStep, Flag_EndStep, FileID, NumberOfIterations )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            !class(ClassFEMAnalysis) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: Time
            real(8), dimension(:) :: U
            integer :: i, FileID, LC, ST, CutBack, SubStep, Flag_EndStep, NumberOfIterations

 		    !************************************************************************************
            ! WRITING RESULTS
		    !************************************************************************************
            write(FileID,*) 'TIME =', Time
            write(FileID,*) 'LOAD CASE =', LC
            write(FileID,*) 'STEP =', ST
            write(FileID,*) 'CUT BACK =', CutBack
            write(FileID,*) 'SUBSTEP =', SubStep
            write(FileID,*) 'FLAG END STEP =', Flag_EndStep
            write(FileID,*) 'NUMBER OF ITERATIONS TO CONVERGE =', NumberOfIterations
            do i = 1,size(U)
                write(FileID,*) U(i)
            enddo

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !##################################################################################################
        ! This routine contains the procedures to solve a quasi-static analysis based in a incremental-
        ! iterative approach.
        !##################################################################################################
        subroutine QuasiStaticAnalysisFEM( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                           Kg , NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ElementLibrary
            use ModAnalysis
            use Nodes
            use BoundaryConditions
            use GlobalSparseMatrix
            use NonLinearSolver
            use Interfaces
            use MathRoutines
            use LoadHistoryData
            use modFEMSystemOfEquations

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type (ClassAnalysis)                                    :: AnalysisSettings
            type (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions),  pointer               :: BC
            type (ClassGlobalSparseMatrix),  pointer                :: Kg
            class(ClassNonLinearSolver),     pointer                :: NLSolver

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: U , R , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0, Uconverged
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: alpha, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases ,  CutBack, SubStep, e,gp, nDOF, FileID_FEMAnalysisResults, Flag_EndStep
            real(8), parameter :: GR= (1.0d0 + dsqrt(5.0d0))/2.0d0

            integer, allocatable, dimension(:) :: KgValZERO, KgValONE
            integer :: contZERO, contONE
            
            type(ClassFEMSystemOfEquations) :: FEMSoE

            FileID_FEMAnalysisResults = 42
            open (FileID_FEMAnalysisResults,file='FEMAnalysis.result',status='unknown')

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOF)

            write(FileID_FEMAnalysisResults,*) 'Total Number of DOF = ', nDOF

            FEMSoE % ElementList => ElementList
            FEMSoE % AnalysisSettings = AnalysisSettings
            FEMSoE % GlobalNodesList => GlobalNodesList
            FEMSoE % BC => BC
            FEMSoE % Kg => Kg
            allocate( FEMSoE % Fint(nDOF) , FEMSoE % Fext(nDOF) , FEMSoE % Ubar(nDOF) )


            ! Allocating arrays
            allocate( R(nDOF), DeltaFext(nDOF),   Fext_alpha0(nDOF) )
            allocate( U(nDOF), DeltaUPresc(nDOF), Ubar_alpha0(nDOF), Uconverged(nDOF)  )


            U = 0.0d0
            Ubar_alpha0 = 0.0d0

            nLoadCases = BC%GetNumberOfLoadCases()

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
            Flag_EndStep = 1
            call WriteFEMResults( U, 0.0d0, 1, 1, 0, 0, Flag_EndStep, FileID_FEMAnalysisResults, NumberOfIterations=0  )


            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    call BC%GetBoundaryConditions(AnalysisSettings, LC, ST, Fext_alpha0, DeltaFext,FEMSoE%DispDOF, U, DeltaUPresc)

                    
                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação
                    ! da CC de deslocamento prescrito
                    !-----------------------------------------------------------------------------------
                    if ( (LC == 1) .and. (ST == 1) ) then

                        allocate( KgValZERO(size(FEMSoE%Kg%Val)), KgValONE(size(FEMSoE%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoE%Kg, FEMSoE%DispDOF, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%PrescDispSparseMapZERO(contZERO), FEMSoE%PrescDispSparseMapONE(contONE) )

                        FEMSoE%PrescDispSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%PrescDispSparseMapONE(:) = KgValONE(1:contONE)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoE%Kg, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%FixedSupportSparseMapZERO(contZERO), FEMSoE%FixedSupportSparseMapONE(contONE) )

                        FEMSoE%FixedSupportSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%FixedSupportSparseMapONE(:) = KgValONE(1:contONE)

                        deallocate( KgValZERO, KgValONE )


                    end if
                    !-----------------------------------------------------------------------------------
                    
                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime)

                    ! Prescribed Incremental Displacement
                    Ubar_alpha0 = U
                    Uconverged = U

                    alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    alpha = alpha_max

                    CutBack = 0 ; SubStep = 0

                    SUBSTEPS: do while(.true.)


                        write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'


                        FEMSoE % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoE % Fext = Fext_alpha0 + alpha*DeltaFext
                        FEMSoE % Ubar = Ubar_alpha0 + alpha*DeltaUPresc


                        call NLSolver%Solve( FEMSoE , XGuess = Uconverged , X = U )

                        IF (NLSolver%Status%Error) then

                            write(*,'(12x,a)') 'Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoE%Status%ErrorDescription)
                            write(*,*)''

                            alpha = alpha_min + (1.0d0-1.0d0/GR)*( alpha - alpha_min )

                            U = Uconverged

                            ! Update Mesh Coordinates
                            if (AnalysisSettings%NLAnalysis == .true.) then
                                call UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,U)
                            endif

                            CutBack = CutBack + 1
                            SubStep = 1
                            if ( CutBack .gt. AnalysisSettings%MaxCutBack ) then
                                write(*,'(a,i3,a,i3,a,i3,a)') 'Load Case: ',LC,' Step: ', ST , ' did not converge with ', AnalysisSettings%MaxCutBack, ' cut backs.'
                                stop
                            endif

                            write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            !---------------------------------------------------------------------------
                        ELSEIF (alpha==1.0d0) then

                            SubStep = SubStep + 1

                            Flag_EndStep = 1
                            call WriteFEMResults( U, FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults, NLSolver%NumberOfIterations )

                            exit SUBSTEPS

                            !---------------------------------------------------------------------------
                        ELSE

                            SubStep = SubStep + 1

                            alpha_aux = alpha_min

                            alpha_min = alpha

                            alpha = min(alpha + GR*(alpha - alpha_aux),1.0d0)

                            Uconverged = U

                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            Flag_EndStep = 0
                            call WriteFEMResults( U, FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults,  NLSolver%NumberOfIterations  )

                        ENDIF


                    enddo SUBSTEPS



                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()
                        enddo
                    enddo
                    ! -----------------------------------------------------------------------------------




                    write(*,'(4x,a,i3)')'End Step: ',ST
                    write(*,*)''

                enddo STEPS

                write(*,'(a,i3)')'End Load Case: ',LC
                write(*,*)''
                write(*,*)''

            enddo LOAD_CASE

            close (FileID_FEMAnalysisResults)
            !************************************************************************************


        end subroutine
        !##################################################################################################


        !##################################################################################################
        ! This routine contains the procedures to solve a quasi-static analysis based in a incremental-
        ! iterative approach.
        !##################################################################################################
        subroutine QuasiStaticAnalysisMultiscaleMinimalFEM( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                                            Kg , NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ElementLibrary
            use ModAnalysis
            use Nodes
            use BoundaryConditions
            use GlobalSparseMatrix
            use NonLinearSolver
            use Interfaces
            use MathRoutines
            use LoadHistoryData
            use ModMultiscaleMinimalFEMSoE

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type  (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions),  pointer                :: BC
            type  (ClassGlobalSparseMatrix),  pointer                :: Kg
            class (ClassNonLinearSolver),     pointer                :: NLSolver

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: X , R , DeltaFext, DeltaUPresc, Fext_alpha0, Xbar_alpha0, Xconverged
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: alpha, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases ,  CutBack, SubStep, e,gp, nDOF, FileID_FEMAnalysisResults, Flag_EndStep
            real(8), parameter :: GR = (1.0d0 + dsqrt(5.0d0))/2.0d0

            integer, allocatable, dimension(:) :: KgValZERO, KgValONE
            integer :: contZERO, contONE            
            
            type(ClassMultiscaleMinimalFEMSoE) :: FEMSoE

            FileID_FEMAnalysisResults = 42
            open (FileID_FEMAnalysisResults,file='FEMAnalysis.result',status='unknown')

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOF)

            write(FileID_FEMAnalysisResults,*) 'Total Number of DOF = ', nDOF

            FEMSoE % ElementList => ElementList
            FEMSoE % AnalysisSettings = AnalysisSettings
            FEMSoE % GlobalNodesList => GlobalNodesList
            FEMSoE % BC => BC
            FEMSoE % Kg => Kg
            allocate( FEMSoE%Fint(nDOF) , FEMSoE%Fext(nDOF) , FEMSoE%Ubar(nDOF) , FEMSoE%Fmacro_current(9) )


            ! Allocating arrays
            allocate( R(nDOF), DeltaFext(nDOF),   Fext_alpha0(nDOF) )
            allocate( X(nDOF+12), DeltaUPresc(nDOF), Xbar_alpha0(nDOF), Xconverged(nDOF+12)  )

           ! Initial Guess
            X = 0.0d0

            nLoadCases = BC%GetNumberOfLoadCases()

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
            Flag_EndStep = 1
            call WriteFEMResults( X(1:nDOF), 0.0d0, 1, 1, 0, 0, Flag_EndStep, FileID_FEMAnalysisResults, NumberOfIterations=0  )


            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    ! Rotina usada somente para transportar o gradiente de deformação macro, no instante corrente (n+1), na variável DeltaFext.
                    !-------------------------------------------------------------
                    call BC%GetBoundaryConditions(AnalysisSettings, LC, ST, Fext_alpha0, DeltaFext,FEMSoE%DispDOF, X, DeltaUPresc)              
    
                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação
                    ! da CC de deslocamento prescrito
                    !-----------------------------------------------------------------------------------
                    if ( (LC == 1) .and. (ST == 1) ) then

                        allocate( KgValZERO(size(FEMSoE%Kg%Val)), KgValONE(size(FEMSoE%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoE%Kg, FEMSoE%DispDOF, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%PrescDispSparseMapZERO(contZERO), FEMSoE%PrescDispSparseMapONE(contONE) )

                        FEMSoE%PrescDispSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%PrescDispSparseMapONE(:) = KgValONE(1:contONE)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoE%Kg, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%FixedSupportSparseMapZERO(contZERO), FEMSoE%FixedSupportSparseMapONE(contONE) )

                        FEMSoE%FixedSupportSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%FixedSupportSparseMapONE(:) = KgValONE(1:contONE)

                        deallocate( KgValZERO, KgValONE )


                    end if
                    !-----------------------------------------------------------------------------------                    
                    FEMSoE%Fmacro_current(1:9) = DeltaFext(1:9)

                    !-------------------------------------------------------------


                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime)

                    ! Prescribed Incremental Displacement
                    Xconverged = X

                    alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    alpha = alpha_max

                    CutBack = 0 ; SubStep = 0

                    SUBSTEPS: do while(.true.)


                        write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'


                        FEMSoE % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoE%Fmacro_current(1:9) = alpha*DeltaFext(1:9)

                        call NLSolver%Solve( FEMSoE , XGuess = Xconverged , X=X )

                        IF (NLSolver%Status%Error) then

                            write(*,'(12x,a)') 'Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoE%Status%ErrorDescription)
                            write(*,*)''

                            alpha = alpha_min + (1.0d0-1.0d0/GR)*( alpha - alpha_min )

                            X = Xconverged

                            ! Update Mesh Coordinates
                            if (AnalysisSettings%NLAnalysis == .true.) then
                                call UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,X)
                            endif

                            CutBack = CutBack + 1
                            SubStep = 1
                            if ( CutBack .gt. AnalysisSettings%MaxCutBack ) then
                                write(*,'(a,i3,a,i3,a,i3,a)') 'Load Case: ',LC,' Step: ', ST , ' did not converge with ', AnalysisSettings%MaxCutBack, ' cut backs.'
                                stop
                            endif

                            write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            !---------------------------------------------------------------------------
                        ELSEIF (alpha==1.0d0) then

                            SubStep = SubStep + 1

                            Flag_EndStep = 1
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults, NLSolver%NumberOfIterations )

                            exit SUBSTEPS

                            !---------------------------------------------------------------------------
                        ELSE

                            SubStep = SubStep + 1

                            alpha_aux = alpha_min

                            alpha_min = alpha

                            alpha = min(alpha + GR*(alpha - alpha_aux),1.0d0)

                            Xconverged = X

                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            Flag_EndStep = 0
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults,  NLSolver%NumberOfIterations  )

                        ENDIF


                    enddo SUBSTEPS



                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()
                        enddo
                    enddo
                    ! -----------------------------------------------------------------------------------




                    write(*,'(4x,a,i3)')'End Step: ',ST
                    write(*,*)''

                enddo STEPS

                write(*,'(a,i3)')'End Load Case: ',LC
                write(*,*)''
                write(*,*)''

            enddo LOAD_CASE

            close (FileID_FEMAnalysisResults)
            !************************************************************************************


        end subroutine
        !##################################################################################################

        !##################################################################################################
        ! This routine contains the procedures to solve a quasi-static analysis based in a incremental-
        ! iterative approach.
        !##################################################################################################
        subroutine QuasiStaticAnalysisMultiscaleMinimalPhasesFEM( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                                            Kg , NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ElementLibrary
            use ModAnalysis
            use Nodes
            use BoundaryConditions
            use GlobalSparseMatrix
            use NonLinearSolver
            use Interfaces
            use MathRoutines
            use LoadHistoryData
            use ModMultiscaleMinimalPhasesFEMSoE

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type  (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions),  pointer                :: BC
            type  (ClassGlobalSparseMatrix),  pointer                :: Kg
            class (ClassNonLinearSolver),     pointer                :: NLSolver

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: X , R , DeltaFext, DeltaUPresc, Fext_alpha0, Xbar_alpha0, Xconverged
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: alpha, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases ,  CutBack, SubStep, e,gp, nDOF, FileID_FEMAnalysisResults, Flag_EndStep
            real(8), parameter :: GR = (1.0d0 + dsqrt(5.0d0))/2.0d0
            integer :: NumberOfMaterials

            integer, allocatable, dimension(:) :: KgValZERO, KgValONE
            integer :: contZERO, contONE            
            
            type(ClassMultiscaleMinimalPhasesFEMSoE) :: FEMSoE

            FileID_FEMAnalysisResults = 42
            open (FileID_FEMAnalysisResults,file='FEMAnalysis.result',status='unknown')

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOF)

            write(FileID_FEMAnalysisResults,*) 'Total Number of DOF = ', nDOF

            FEMSoE % ElementList => ElementList
            FEMSoE % AnalysisSettings = AnalysisSettings
            FEMSoE % GlobalNodesList => GlobalNodesList
            FEMSoE % BC => BC
            FEMSoE % Kg => Kg
            allocate( FEMSoE%Fint(nDOF) , FEMSoE%Fext(nDOF) , FEMSoE%Ubar(nDOF) , FEMSoE%Fmacro_current(9) )

            NumberOfMaterials = AnalysisSettings%NumberOfMaterials

            ! Allocating arrays
            allocate( R(nDOF), DeltaFext(nDOF),   Fext_alpha0(nDOF) )
            allocate( X(nDOF+9*NumberOfMaterials+3), DeltaUPresc(nDOF), Xbar_alpha0(nDOF), Xconverged(nDOF+9*NumberOfMaterials+3)  )

           ! Initial Guess
            X = 0.0d0

            nLoadCases = BC%GetNumberOfLoadCases()

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
            Flag_EndStep = 1
            call WriteFEMResults( X(1:nDOF), 0.0d0, 1, 1, 0, 0, Flag_EndStep, FileID_FEMAnalysisResults, NumberOfIterations=0  )


            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    ! Rotina usada somente para transportar o gradiente de deformação macro, no instante corrente (n+1), na variável DeltaFext.
                    !-------------------------------------------------------------
                    call BC%GetBoundaryConditions(AnalysisSettings, LC, ST, Fext_alpha0, DeltaFext,FEMSoE%DispDOF, X, DeltaUPresc)

                    
                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação
                    ! da CC de deslocamento prescrito
                    !-----------------------------------------------------------------------------------
                    if ( (LC == 1) .and. (ST == 1) ) then

                        allocate( KgValZERO(size(FEMSoE%Kg%Val)), KgValONE(size(FEMSoE%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoE%Kg, FEMSoE%DispDOF, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%PrescDispSparseMapZERO(contZERO), FEMSoE%PrescDispSparseMapONE(contONE) )

                        FEMSoE%PrescDispSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%PrescDispSparseMapONE(:) = KgValONE(1:contONE)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoE%Kg, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%FixedSupportSparseMapZERO(contZERO), FEMSoE%FixedSupportSparseMapONE(contONE) )

                        FEMSoE%FixedSupportSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%FixedSupportSparseMapONE(:) = KgValONE(1:contONE)

                        deallocate( KgValZERO, KgValONE )


                    end if
                    !-----------------------------------------------------------------------------------               
                    FEMSoE%Fmacro_current(1:9) = DeltaFext(1:9)

                    !-------------------------------------------------------------


                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime)

                    ! Prescribed Incremental Displacement
                    Xconverged = X

                    alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    alpha = alpha_max

                    CutBack = 0 ; SubStep = 0

                    SUBSTEPS: do while(.true.)


                        write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'


                        FEMSoE % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoE%Fmacro_current(1:9) = alpha*DeltaFext(1:9)

                        call NLSolver%Solve( FEMSoE , XGuess = Xconverged , X=X )

                        IF (NLSolver%Status%Error) then

                            write(*,'(12x,a)') 'Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoE%Status%ErrorDescription)
                            write(*,*)''

                            alpha = alpha_min + (1.0d0-1.0d0/GR)*( alpha - alpha_min )

                            X = Xconverged

                            ! Update Mesh Coordinates
                            if (AnalysisSettings%NLAnalysis == .true.) then
                                call UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,X)
                            endif

                            CutBack = CutBack + 1
                            SubStep = 1
                            if ( CutBack .gt. AnalysisSettings%MaxCutBack ) then
                                write(*,'(a,i3,a,i3,a,i3,a)') 'Load Case: ',LC,' Step: ', ST , ' did not converge with ', AnalysisSettings%MaxCutBack, ' cut backs.'
                                stop
                            endif

                            write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            !---------------------------------------------------------------------------
                        ELSEIF (alpha==1.0d0) then

                            SubStep = SubStep + 1

                            Flag_EndStep = 1
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults, NLSolver%NumberOfIterations )

                            exit SUBSTEPS

                            !---------------------------------------------------------------------------
                        ELSE

                            SubStep = SubStep + 1

                            alpha_aux = alpha_min

                            alpha_min = alpha

                            alpha = min(alpha + GR*(alpha - alpha_aux),1.0d0)

                            Xconverged = X

                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            Flag_EndStep = 0
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults,  NLSolver%NumberOfIterations  )

                        ENDIF


                    enddo SUBSTEPS



                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()
                        enddo
                    enddo
                    ! -----------------------------------------------------------------------------------




                    write(*,'(4x,a,i3)')'End Step: ',ST
                    write(*,*)''

                enddo STEPS

                write(*,'(a,i3)')'End Load Case: ',LC
                write(*,*)''
                write(*,*)''

            enddo LOAD_CASE

            close (FileID_FEMAnalysisResults)
            !************************************************************************************


        end subroutine
        !##################################################################################################

        !##################################################################################################
        ! This routine contains the procedures to solve a quasi-static analysis based in a incremental-
        ! iterative approach.
        !##################################################################################################
        subroutine QuasiStaticAnalysisMultiscaleMinimalLinearD1FEM( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                                            Kg , NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ElementLibrary
            use ModAnalysis
            use Nodes
            use BoundaryConditions
            use GlobalSparseMatrix
            use NonLinearSolver
            use Interfaces
            use MathRoutines
            use LoadHistoryData
            use ModMultiscaleMinimalLinearD1FEMSoE

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type  (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions),  pointer                :: BC
            type  (ClassGlobalSparseMatrix),  pointer                :: Kg
            class (ClassNonLinearSolver),     pointer                :: NLSolver

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: X , R , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0, Xconverged
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: alpha, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases ,  CutBack, SubStep, e,gp, nDOF, FileID_FEMAnalysisResults, Flag_EndStep
            real(8), parameter :: GR = (1.0d0 + dsqrt(5.0d0))/2.0d0

            integer, allocatable, dimension(:) :: KgValZERO, KgValONE
            integer :: contZERO, contONE            
            
            type(ClassMultiscaleMinimalLinearD1FEMSoE) :: FEMSoE

            FileID_FEMAnalysisResults = 42
            open (FileID_FEMAnalysisResults,file='FEMAnalysis.result',status='unknown')

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOF)

            write(FileID_FEMAnalysisResults,*) 'Total Number of DOF = ', nDOF

            FEMSoE % ElementList => ElementList
            FEMSoE % AnalysisSettings = AnalysisSettings
            FEMSoE % GlobalNodesList => GlobalNodesList
            FEMSoE % BC => BC
            FEMSoE % Kg => Kg
            allocate( FEMSoE%Fint(nDOF) , FEMSoE%Fext(nDOF) , FEMSoE%Ubar(nDOF) , FEMSoE%Fmacro_current(9) )


            ! Allocating arrays
            allocate( R(nDOF), DeltaFext(nDOF),   Fext_alpha0(nDOF) )
            allocate( X(nDOF+12), DeltaUPresc(nDOF), Ubar_alpha0(nDOF), Xconverged(nDOF+12)  )


           ! Initial Guess
            X = 0.0d0
            Ubar_alpha0 = 0.0d0

            nLoadCases = BC%GetNumberOfLoadCases()

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
            Flag_EndStep = 1
            call WriteFEMResults( X(1:nDOF), 0.0d0, 1, 1, 0, 0, Flag_EndStep, FileID_FEMAnalysisResults, NumberOfIterations=0  )


            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    ! Rotina usada somente para transportar o gradiente de deformação macro, no instante corrente (n+1), na variável DeltaFext.
                    !-------------------------------------------------------------
                    call BC%GetBoundaryConditions(AnalysisSettings, LC, ST, Fext_alpha0, DeltaFext,FEMSoE%DispDOF, X, DeltaUPresc)

                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação
                    ! da CC de deslocamento prescrito
                    !-----------------------------------------------------------------------------------
                    if ( (LC == 1) .and. (ST == 1) ) then

                        allocate( KgValZERO(size(FEMSoE%Kg%Val)), KgValONE(size(FEMSoE%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoE%Kg, FEMSoE%DispDOF, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%PrescDispSparseMapZERO(contZERO), FEMSoE%PrescDispSparseMapONE(contONE) )

                        FEMSoE%PrescDispSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%PrescDispSparseMapONE(:) = KgValONE(1:contONE)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoE%Kg, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%FixedSupportSparseMapZERO(contZERO), FEMSoE%FixedSupportSparseMapONE(contONE) )

                        FEMSoE%FixedSupportSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%FixedSupportSparseMapONE(:) = KgValONE(1:contONE)

                        deallocate( KgValZERO, KgValONE )


                    end if
                    !-----------------------------------------------------------------------------------                 
                    
                    FEMSoE%Fmacro_current(1:9) = DeltaFext(1:9)

                    !-------------------------------------------------------------
                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime)

                    ! Prescribed Incremental Displacement
                    Ubar_alpha0 = X(1:nDOF)
                    Xconverged = X

                    alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    alpha = alpha_max

                    CutBack = 0 ; SubStep = 0

                    SUBSTEPS: do while(.true.)


                        write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                        FEMSoE % Ubar = Ubar_alpha0 + alpha*DeltaUPresc
                        FEMSoE % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoE%Fmacro_current(1:9) = alpha*DeltaFext(1:9)

                        call NLSolver%Solve( FEMSoE , XGuess = Xconverged , X=X )

                        IF (NLSolver%Status%Error) then

                            write(*,'(12x,a)') 'Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoE%Status%ErrorDescription)
                            write(*,*)''

                            alpha = alpha_min + (1.0d0-1.0d0/GR)*( alpha - alpha_min )

                            X = Xconverged

                            ! Update Mesh Coordinates
                            if (AnalysisSettings%NLAnalysis == .true.) then
                                call UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,X)
                            endif

                            CutBack = CutBack + 1
                            SubStep = 1
                            if ( CutBack .gt. AnalysisSettings%MaxCutBack ) then
                                write(*,'(a,i3,a,i3,a,i3,a)') 'Load Case: ',LC,' Step: ', ST , ' did not converge with ', AnalysisSettings%MaxCutBack, ' cut backs.'
                                stop
                            endif

                            write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            !---------------------------------------------------------------------------
                        ELSEIF (alpha==1.0d0) then

                            SubStep = SubStep + 1

                            Flag_EndStep = 1
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults, NLSolver%NumberOfIterations )

                            exit SUBSTEPS

                            !---------------------------------------------------------------------------
                        ELSE

                            SubStep = SubStep + 1

                            alpha_aux = alpha_min

                            alpha_min = alpha

                            alpha = min(alpha + GR*(alpha - alpha_aux),1.0d0)

                            Xconverged = X

                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            Flag_EndStep = 0
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults,  NLSolver%NumberOfIterations  )

                        ENDIF


                    enddo SUBSTEPS



                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()
                        enddo
                    enddo
                    ! -----------------------------------------------------------------------------------




                    write(*,'(4x,a,i3)')'End Step: ',ST
                    write(*,*)''

                enddo STEPS

                write(*,'(a,i3)')'End Load Case: ',LC
                write(*,*)''
                write(*,*)''

            enddo LOAD_CASE

            close (FileID_FEMAnalysisResults)
            !************************************************************************************


                                                            end subroutine
        !##################################################################################################

        !##################################################################################################
        ! This routine contains the procedures to solve a quasi-static analysis based in a incremental-
        ! iterative approach.
        !##################################################################################################
        subroutine QuasiStaticAnalysisMultiscaleMinimalLinearD3FEM( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                                            Kg , NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ElementLibrary
            use ModAnalysis
            use Nodes
            use BoundaryConditions
            use GlobalSparseMatrix
            use NonLinearSolver
            use Interfaces
            use MathRoutines
            use LoadHistoryData
            use ModMultiscaleMinimalLinearD3FEMSoE

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type  (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions),  pointer                :: BC
            type  (ClassGlobalSparseMatrix),  pointer                :: Kg
            class (ClassNonLinearSolver),     pointer                :: NLSolver

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: X , R , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0, Xconverged
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: alpha, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases ,  CutBack, SubStep, e,gp, nDOF, FileID_FEMAnalysisResults, Flag_EndStep
            real(8), parameter :: GR = (1.0d0 + dsqrt(5.0d0))/2.0d0

            integer, allocatable, dimension(:) :: KgValZERO, KgValONE
            integer :: contZERO, contONE            
            
            type(ClassMultiscaleMinimalLinearD3FEMSoE) :: FEMSoE

            FileID_FEMAnalysisResults = 42
            open (FileID_FEMAnalysisResults,file='FEMAnalysis.result',status='unknown')

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOF)

            write(FileID_FEMAnalysisResults,*) 'Total Number of DOF = ', nDOF

            FEMSoE % ElementList => ElementList
            FEMSoE % AnalysisSettings = AnalysisSettings
            FEMSoE % GlobalNodesList => GlobalNodesList
            FEMSoE % BC => BC
            FEMSoE % Kg => Kg
            allocate( FEMSoE%Fint(nDOF) , FEMSoE%Fext(nDOF) , FEMSoE%Ubar(nDOF) , FEMSoE%Fmacro_current(9) )


            ! Allocating arrays
            allocate( R(nDOF), DeltaFext(nDOF),   Fext_alpha0(nDOF) )
            allocate( X(nDOF+12), DeltaUPresc(nDOF), Ubar_alpha0(nDOF), Xconverged(nDOF+12)  )


           ! Initial Guess
            X = 0.0d0
            Ubar_alpha0 = 0.0d0

            nLoadCases = BC%GetNumberOfLoadCases()

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
            Flag_EndStep = 1
            call WriteFEMResults( X(1:nDOF), 0.0d0, 1, 1, 0, 0, Flag_EndStep, FileID_FEMAnalysisResults, NumberOfIterations=0  )


            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    ! Rotina usada somente para transportar o gradiente de deformação macro, no instante corrente (n+1), na variável DeltaFext.
                    !-------------------------------------------------------------
                    call BC%GetBoundaryConditions(AnalysisSettings, LC, ST, Fext_alpha0, DeltaFext,FEMSoE%DispDOF, X, DeltaUPresc)

                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação
                    ! da CC de deslocamento prescrito
                    !-----------------------------------------------------------------------------------
                    if ( (LC == 1) .and. (ST == 1) ) then

                        allocate( KgValZERO(size(FEMSoE%Kg%Val)), KgValONE(size(FEMSoE%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoE%Kg, FEMSoE%DispDOF, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%PrescDispSparseMapZERO(contZERO), FEMSoE%PrescDispSparseMapONE(contONE) )

                        FEMSoE%PrescDispSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%PrescDispSparseMapONE(:) = KgValONE(1:contONE)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoE%Kg, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%FixedSupportSparseMapZERO(contZERO), FEMSoE%FixedSupportSparseMapONE(contONE) )

                        FEMSoE%FixedSupportSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%FixedSupportSparseMapONE(:) = KgValONE(1:contONE)

                        deallocate( KgValZERO, KgValONE )


                    end if
                    !-----------------------------------------------------------------------------------                 
                    
                    FEMSoE%Fmacro_current(1:9) = DeltaFext(1:9)

                    !-------------------------------------------------------------
                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime)

                    ! Prescribed Incremental Displacement
                    Ubar_alpha0 = X(1:nDOF)
                    Xconverged = X

                    alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    alpha = alpha_max

                    CutBack = 0 ; SubStep = 0

                    SUBSTEPS: do while(.true.)


                        write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                        FEMSoE % Ubar = Ubar_alpha0 + alpha*DeltaUPresc
                        FEMSoE % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoE%Fmacro_current(1:9) = alpha*DeltaFext(1:9)

                        call NLSolver%Solve( FEMSoE , XGuess = Xconverged , X=X )

                        IF (NLSolver%Status%Error) then

                            write(*,'(12x,a)') 'Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoE%Status%ErrorDescription)
                            write(*,*)''

                            alpha = alpha_min + (1.0d0-1.0d0/GR)*( alpha - alpha_min )

                            X = Xconverged

                            ! Update Mesh Coordinates
                            if (AnalysisSettings%NLAnalysis == .true.) then
                                call UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,X)
                            endif

                            CutBack = CutBack + 1
                            SubStep = 1
                            if ( CutBack .gt. AnalysisSettings%MaxCutBack ) then
                                write(*,'(a,i3,a,i3,a,i3,a)') 'Load Case: ',LC,' Step: ', ST , ' did not converge with ', AnalysisSettings%MaxCutBack, ' cut backs.'
                                stop
                            endif

                            write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            !---------------------------------------------------------------------------
                        ELSEIF (alpha==1.0d0) then

                            SubStep = SubStep + 1

                            Flag_EndStep = 1
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults, NLSolver%NumberOfIterations )

                            exit SUBSTEPS

                            !---------------------------------------------------------------------------
                        ELSE

                            SubStep = SubStep + 1

                            alpha_aux = alpha_min

                            alpha_min = alpha

                            alpha = min(alpha + GR*(alpha - alpha_aux),1.0d0)

                            Xconverged = X

                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            Flag_EndStep = 0
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults,  NLSolver%NumberOfIterations  )

                        ENDIF


                    enddo SUBSTEPS



                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()
                        enddo
                    enddo
                    ! -----------------------------------------------------------------------------------




                    write(*,'(4x,a,i3)')'End Step: ',ST
                    write(*,*)''

                enddo STEPS

                write(*,'(a,i3)')'End Load Case: ',LC
                write(*,*)''
                write(*,*)''

            enddo LOAD_CASE

            close (FileID_FEMAnalysisResults)
            !************************************************************************************


        end subroutine
        !##################################################################################################                                                            


        !==========================================================================================
        ! Method ClassFEMAnalysis:
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  AdditionalMaterialModelRoutine( this )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use MathRoutines
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class (ClassFEMAnalysis) :: this


            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: R, L, pitch, hand, theta, Xgp, X0ref, tXgp, norm_mX, A0, L0, w
            real(8) :: mX(3), IP(3), NodalValuesX(50)
            integer :: ElemRef, NodeRef, e, gp, n, NumberOfNodes, i, j, nGP, El_ID
            
            real(8),dimension(9) :: FiberData
            logical              :: file_exists


            real(8) , pointer , dimension(:,:) :: NaturalCoord
            real(8) , pointer , dimension(:)   :: Weight
            type(ClassElementProfile)          :: ElProfile

 		    !************************************************************************************
            ! ADDITIONAL COMPUTATIONS ON GAUSS POINTS
		    !************************************************************************************
! TODO (Thiago#1#): Passar os enumeradores dos modelos para realizar as contas somente nos elementos que possuem o devido modelo material.

            !####################################################################################
            ! Cálculo das direções dos reforços de fibra
            !####################################################################################
            
            inquire(file='Fiber_info.tab',exist=file_exists)
            
            if (.not.file_exists) then
            
                write(*,*) 'File Fiber_info.tab not found'
                STOP
            
            else
            
                do e = 1 , size(this%ElementList)
                
                    !Call profile to check if element has reinforcement capabilities
                    call this%ElementList(e)%El%GetProfile(ElProfile)
    
                    if (ElProfile%AcceptFiberReinforcement == .true.) then
                
                        call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)
            
                        do gp = 1,size(NaturalCoord,dim=1) !matrix Gauss points - mX=0

                            mX(1) = 0.0d0
                            mX(2) = 0.0d0
                            mX(3) = 0.0d0
                            A0 = 0.0d0
                            L0 = 0.0d0
            
                            this%ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%mX = mX
                            this%ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%A0 = A0
                            this%ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%L0 = L0
                    
                        enddo
                
                        do gp = 1,size(this%ElementList(e)%El%ExtraGaussPoints) !fibers Gauss points
                            
                            open(87,file='Fiber_info.tab',status='old')
                
                            do i=1,e                
                                read(87,*)
                                read(87,*)
                                read(87,'(i)') El_ID
                                read(87,*)
                                read(87,'(i)') nGP
                                read(87,*)
                    
                                if (i==e) then
                                                
                                    do j=1,nGP
                                        if (j==gp) then
                                            read(87,*) FiberData(:)
                                            
                                            IP(1) = FiberData(1)
                                            IP(2) = FiberData(2)
                                            IP(3) = FiberData(3)
                                            w = FiberData(4)
                                            mX(1) = FiberData(5)
                                            mX(2) = FiberData(6)
                                            mX(3) = FiberData(7)
                                            L0 = FiberData(8)
                                            A0 = FiberData(9)
                                            
                                            this%ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%NaturalCoord = IP
                                            this%ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%Weight = w

                                            this%ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%mX = mX
                                            this%ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%A0 = A0
                                            this%ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%L0 = L0
                                            
                                        elseif (nGP /= 0) then
                                            read(87,*)
                                        endif
                                    enddo
                        
                                elseif (nGP /= 0) then
                        
                                do j=1,nGP
                                    read(87,*)
                                enddo

                                endif
                                    
                            enddo
                
                            close(87)
                    
                        enddo
                    
                    endif
                                
                enddo
                
            endif
            
            
            !####################################################################################
            ! Cálculo das tangentes da hélice
            !####################################################################################

    if ( n .eq. 1234123312 ) then

            ! Parâmetros da Hélice
            R = 0.05d0          !1 fibra de 3 voltas: 2.30d-3
            L = 1.00d0           !1 fibra de 3 voltas: 297.90d-3
            pitch =  1.0d0
            hand  =  -1.0d0
            theta =  90.d0 !270.0d0   !CUIDAR A ORDEM DO DESENHO NO SOLIDWORKS!!!!!   !1 fibra de 3 voltas: 0.0d0

            ! Elemento e Nó de Referência
            ElemRef = 120   !61608 !3228  !93697 !1 !16308  !1          !211 !166 !45271 !15226 !4401     !2 Fibras: 1776  !12 Fibras: 4526       !1 fibra de 3 voltas: 1501
            NodeRef = 256   !70457 !7357  !1667  !75 !19670  !75         !1064 !110 !289 !5150 !2088     !2 Fibras: 2868  !12 Fibras: 5950       !1 fibra de 3 voltas: 31

            !Obtendo o ID do Nó de Referência
            NumberOfNodes =  this%ElementList(ElemRef)%El%GetNumberOfNodes()
            do n = 1,NumberOfNodes
                if (this%ElementList(ElemRef)%El%ElementNodes(n)%Node%ID .eq. NodeRef) then
                    NodeRef = n
                    exit
                endif
            enddo


            ! Cálculando a tangente nos pontos de Gauss (Para toda a malha!!!!)
            NodalValuesX = 0.0d0
            do e = 1 , size(this%ElementList)

                NumberOfNodes = this%ElementList(e)%El%GetNumberOfNodes()

                do gp = 1,size(this%ElementList(e)%El%GaussPoints)

                    ! Obtendo as coordenadas nodais X do elemento
                    do n = 1,NumberOfNodes
                        NodalValuesX(n) = this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(1)
                    enddo

                    call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                    call this%ElementList(e)%El%ElementInterpolation(NodalValuesX(1:NumberOfNodes),NaturalCoord(gp,:),Xgp)

                    ! Coordenada X do nó de referência (ponto onde o parâmetro t=0)
                    X0ref = this%ElementList(ElemRef)%El%ElementNodes(NodeRef)%Node%CoordX(1)

                    tXgp = (2.0d0*Pi*pitch/L)*( Xgp - X0ref )

                    mX(1) =  L/(2.0d0*Pi*pitch )
                    mX(2) = -R*dsin( tXgp + (theta*Pi/180.0d0) )
                    mX(3) =  hand*R*dcos( tXgp + (theta*Pi/180.0d0) )

                    norm_mX = ( mX(1)*mX(1)+mX(2)*mX(2)+mX(3)*mX(3) )**0.50d0
                    mX = mX/norm_mX

                    ! Fibras Retas
                    !----------------------
                    !mX(1) = 1.0d0
                    !mX(2) = 0.0d0
                    !mX(3) = 0.0d0
                    !-----------------------

                    this%ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%mX = mX

                enddo
            enddo
            !####################################################################################

    endif

 		    !************************************************************************************


        end subroutine
        !==========================================================================================


end module


































