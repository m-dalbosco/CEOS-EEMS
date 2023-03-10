!##################################################################################################
! This module contains the explicit interfaces used in the code.
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
module Interfaces

    !==============================================================================================
    interface

        !==========================================================================================
        subroutine ElementConstructor( this, ElementNodes, ElementType, GlobalNodesList, MaterialID )

	        !************************************************************************************
	        ! DECLARATIONS OF VARIABLES
	        !************************************************************************************
	        ! Modules and implicit declarations
	        ! -----------------------------------------------------------------------------------
	        use ElementLibrary
	        use Nodes

	        implicit none

	        ! Object
	        ! -----------------------------------------------------------------------------------
	        class(ClassElement) , pointer :: this

	        ! Input variables
	        ! -----------------------------------------------------------------------------------
	        type(ClassNodes) , dimension(:) , pointer , intent(in) :: GlobalNodesList
	        integer				 , intent(in) :: ElementType
	        integer,dimension(:) , intent(in) :: ElementNodes
            
            ! Output Variables
            integer :: MaterialID

        end subroutine
        !==========================================================================================


        !==========================================================================================
        subroutine AssembleGlobalMatrix( GM , Ke , Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use GlobalSparseMatrix
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(in)   :: GM
            real(8) , dimension(:,:) , intent(in) :: Ke

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            type(ClassGlobalSparseMatrix) :: Kg

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine AssembleGlobalMatrixUpperTriangular( GM , Ke , Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use GlobalSparseMatrix
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(in)   :: GM
            real(8) , dimension(:,:) , intent(in) :: Ke

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            type(ClassGlobalSparseMatrix) :: Kg

        end subroutine
        !==========================================================================================
        
        
        !==========================================================================================
        !subroutine ApplyBoundaryConditions( Kg , R , Presc_Disp_DOF , Fixed_Disp_DOF , Ubar , U )
        !
        !    !************************************************************************************
        !    ! DECLARATIONS OF VARIABLES
        !    !************************************************************************************
        !    ! Modules and implicit declarations
        !    ! -----------------------------------------------------------------------------------
        !    use GlobalSparseMatrix
        !    implicit none
        !
        !    ! Input variables
        !    ! -----------------------------------------------------------------------------------
        !    integer , dimension(:) , intent(in) ::  Presc_Disp_DOF , Fixed_Disp_DOF
        !
        !    ! Input/Output variables
        !    ! -----------------------------------------------------------------------------------
        !    real(8) , dimension(:) , intent(inout) ::  R  , Ubar , U
        !    type(ClassGlobalSparseMatrix) :: Kg
        !
        !end subroutine
        !==========================================================================================

        !==========================================================================================
        !subroutine ExportResultFile( Time , u, nDOFnode, TotalnNodes , ElementList , GiDFile)
        !    use Element
        !    use GiDResultFile
        !    implicit none
        !    type (ClassElementsWrapper) , dimension(:) :: ElementList
        !    type(ClassGiDResultFile)                                       :: GidFile
        !    integer :: nDOFnode, TotalnNodes
        !    real(8) , dimension(:) :: u
        !    Real(8)::Time
        !end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine SolveConstitutiveModel( ElementList , AnalysisSettings , Time, U, Status)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ElementLibrary
            use ModAnalysis
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:) :: ElementList
            type(ClassAnalysis)                       :: AnalysisSettings
            type(ClassStatus)                         :: Status
            real(8)                    , dimension(:) :: U

            real(8)                                   :: Time
            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine InternalForce( ElementList, AnalysisSettings, Fint, Status )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ElementLibrary
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:)  :: Fint

        end subroutine
        !==========================================================================================

        !==========================================================================================
!        subroutine ExternalForce( BC, LC, nSteps, Fext, DeltaFext )
!
!            !************************************************************************************
!            ! DECLARATIONS OF VARIABLES
!            !************************************************************************************
!            ! Modules and implicit declarations
!            ! -----------------------------------------------------------------------------------
!            use BoundaryConditions
!
!            implicit none
!
!            ! Input variables
!            ! -----------------------------------------------------------------------------------
!            type(ClassBoundaryConditions)  , intent(in)  :: BC
!            integer                        , intent(in)  :: nSteps , LC
!
!            ! Output variables
!            ! -----------------------------------------------------------------------------------
!            real(8) , dimension(:)  , intent(out)  :: Fext , DeltaFext
!
!        end subroutine
        !==========================================================================================

        !==========================================================================================
!        subroutine PrescribedDisplacement ( BC , LC , ST, NodalDispDOF, U, DeltaUPresc )
!
!            !************************************************************************************
!            ! DECLARATIONS OF VARIABLES
!            !************************************************************************************
!            ! Modules and implicit declarations
!            ! -----------------------------------------------------------------------------------
!            use BoundaryConditions
!
!            implicit none
!
!            ! Input variables
!            ! -----------------------------------------------------------------------------------
!            type(ClassBoundaryConditions)  :: BC
!            integer                        :: LC, ST
!
!            ! Output variables
!            ! -----------------------------------------------------------------------------------
!            real(8) , dimension(:)              :: U, DeltaUPresc
!            integer , pointer , dimension(:)    :: NodalDispDOF
!
!            !************************************************************************************
!        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine TangentStiffnessMatrix( AnalysisSettings , ElementList , nDOF, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ElementLibrary
            use GlobalSparseMatrix

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(in) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in) :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in) :: Kg
            integer :: nDOF

            !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        subroutine UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,DeltaU)
             use Nodes
             use ModAnalysis
             implicit none
             type (ClassNodes)    , pointer , dimension(:) :: GlobalNodesList
             real(8) , dimension(:)                        :: DeltaU
             type(ClassAnalysis)                           :: AnalysisSettings
        end subroutine
        !==============================================================================================


        !==============================================================================================
        subroutine MicroscopicDisplacement ( AnalysisSettings, GlobalNodesList, U, Umicro )

             use Nodes
             use ModAnalysis

             implicit none

             type(ClassAnalysis)                           :: AnalysisSettings
             type (ClassNodes)    , pointer , dimension(:) :: GlobalNodesList
             real(8) , dimension(:)                        :: U, Umicro

         endsubroutine
        !==============================================================================================


        !==============================================================================================
        subroutine MaterialConstructor( Element, GlobalNodesList, Material, AnalysisSettings, e )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ElementLibrary
            use Nodes
            use ConstitutiveModelLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassElement) , pointer                         :: Element
            type(ClassAnalysis)                                   :: AnalysisSettings
            type(ClassNodes) , dimension(:) , pointer             :: GlobalNodesList
            class(ClassConstitutiveModelWrapper)  , pointer       :: Material
            integer                                               :: e

        end subroutine
        !==============================================================================================


        !==============================================================================================
        subroutine ExternalForceMultiscaleMinimal( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ElementLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u


            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

        end subroutine
        !==============================================================================================

        !==============================================================================================
        subroutine ExternalForceMultiscaleMinimalPhases( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext, PhaseVolX, TotalVolX )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ElementLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u
            real(8)                    , dimension(:)  :: PhaseVolX
            real(8)                                    :: TotalVolX
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext
        
        end subroutine
        !==============================================================================================
        
        !==============================================================================================
        subroutine ExternalForceMultiscaleMinimalLinearD1( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ElementLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u


            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

        end subroutine
        !==============================================================================================

        !==============================================================================================
        subroutine ExternalForceMultiscaleMinimalLinearD3( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ElementLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u


            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

        end subroutine
        !==============================================================================================

        
    end interface



end module
