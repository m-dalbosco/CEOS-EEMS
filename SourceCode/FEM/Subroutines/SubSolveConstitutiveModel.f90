!##################################################################################################
! This routine solves the constitutive equations.
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
subroutine SolveConstitutiveModel( ElementList , AnalysisSettings, Time, U, Status)

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use ElementLibrary
    use ModAnalysis
    use ConstitutiveModel
    use ModStatus


    implicit none

    ! Input variables
    ! -----------------------------------------------------------------------------------
    type(ClassElementsWrapper) , dimension(:)  :: ElementList
    type(ClassAnalysis)                        :: AnalysisSettings
    type(ClassStatus)                          :: Status
    real(8)                    , dimension(:)  :: U
    real(8)                                    :: Time

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    real(8) :: F(3,3)
    real(8) :: Volume, VolumeX, T
    integer :: e , gp , nDOFel, i, j, nNodes
    integer , pointer , dimension(:)   :: GM
    real(8) , pointer , dimension(:,:) :: NaturalCoord
    real(8) , pointer , dimension(:)   :: Weight
    real(8)                            :: ExtraNaturalCoord(3)
    logical                            :: ElementError_exists

    !************************************************************************************

    !************************************************************************************
    ! SOLVING THE GLOBAL CONSTITUTIVE MODEL
    !************************************************************************************


    ! Loop over the elements
    do e = 1 , size(ElementList)

        call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings , nDOFel)
        GM => GM_Memory( 1:nDOFel )

        call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)
        
        call ElementList(e)%El%GetGlobalMapping(AnalysisSettings,GM)

        call ElementList(e)%El%ElementVolume(AnalysisSettings,Volume,VolumeX, Status)

        if (Status%Error) then

            inquire(file='ErrorElement.txt',exist=ElementError_exists) !Check if file with failed elements already exists
            
            if (ElementError_exists) then
                open(37,file='ErrorElement.txt',status='old',access='append')
            else
                open(37,file='ErrorElement.txt',status='new')
            endif
                        
            write(37,'(a)') Trim(Status%ErrorDescription)
            write(37,'(a,i6,a,f7.4)') 'Failed in element ',e,' at time ',Time
            write(37,'(a)') 'Nodal coordinates (X Y Z):'
            
            nNodes = ElementList(e)%El%GetNumberOfNodes()
            
            do i=1,nNodes
                write(37,'(3(1X,f12.6))') (ElementList(e)%El%ElementNodes(i)%Node%Coord(j), j=1,3)
            enddo
            
            write (37,*) ''   
            
            close(37)
            return
            
        endif

        ! Armazendo o volume de cada elemento
        ElementList(e)%El%Volume  = Volume
        ElementList(e)%El%VolumeX = VolumeX

        ! Loop over the Gauss Points
        do gp = 1 , size(ElementList(e)%El%GaussPoints)

            call ElementList(e)%El%DeformationGradient( NaturalCoord(gp,:) , U(GM) , &
                                                        AnalysisSettings , F, Status )

            ElementList(e)%El%GaussPoints(gp)%F = F

            ! AdditionalVariables
            !----------------------------------------------------------------------------
            ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%Jbar = Volume/VolumeX
            !----------------------------------------------------------------------------

            ElementList(e)%El%GaussPoints(gp)%Time = Time

            call ElementList(e)%El%GaussPoints(gp)%UpdateStressAndStateVariables(Status)

        enddo
        
    !************************************************************************************
    ! SOLVING THE GLOBAL CONSTITUTIVE MODEL - EXTRA GAUSS POINTS
    !************************************************************************************
                            
        if (AnalysisSettings%EmbeddedElements) then
        
            do gp = 1 , size(ElementList(e)%El%ExtraGaussPoints)
                
                ExtraNaturalCoord = ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%NaturalCoord

                call ElementList(e)%El%DeformationGradient( ExtraNaturalCoord , U(GM) , AnalysisSettings , F, Status )

                ElementList(e)%El%ExtraGaussPoints(gp)%F = F

                ! AdditionalVariables
                !----------------------------------------------------------------------------
                ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%Jbar = Volume/VolumeX
                !----------------------------------------------------------------------------

                ElementList(e)%El%ExtraGaussPoints(gp)%Time = Time

                call ElementList(e)%El%ExtraGaussPoints(gp)%UpdateStressAndStateVariables(Status)

            enddo
        
        endif

    enddo

end subroutine


