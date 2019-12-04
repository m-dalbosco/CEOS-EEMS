!##################################################################################################
! This routine assembles the nodal contributions of the global internal force.
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

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: e , nDOFel
    integer , pointer , dimension(:) :: GM
    real(8) , pointer , dimension(:) :: Fe
    real(8) , pointer , dimension(:,:) :: Ge
    real(8) , pointer , dimension(:,:) :: Ne

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE INTERNAL FORCE
    !************************************************************************************

    Fext=0.0d0

     do  e = 1, size( ElementList )

        call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings , nDOFel)

        Fe => Fe_Memory( 1:nDOFel )
        Ge => Ge_Memory( 1:9 , 1:nDOFel )
        Ne => Ne_Memory( 1:3 , 1:nDOFel )
        GM => GM_Memory( 1:nDOFel )

        call ElementList(e)%El%GetGlobalMapping(AnalysisSettings,GM)

        call ElementList(e)%El%Matrix_Ne_and_Ge(AnalysisSettings, Ne, Ge)

        Fe = matmul(transpose(Ge),Lambda_F) + matmul(transpose(Ne),Lambda_u)
        !Fe = matmul(transpose(Ge),Lambda_F)

        Fext(GM) = Fext(GM) + Fe

    enddo

    !************************************************************************************

end subroutine

