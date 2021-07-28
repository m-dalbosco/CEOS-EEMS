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
    real(8) , dimension(:) :: Fint

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: e , nDOFel, i, j, nNodes
    integer , pointer , dimension(:) :: GM
    real(8) , pointer , dimension(:) :: Fe
    logical :: ElementError_exists

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE INTERNAL FORCE
    !************************************************************************************

    Fint=0.0d0

     do  e = 1, size( ElementList )

        call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings , nDOFel)

        Fe => Fe_Memory( 1:nDOFel )
        GM => GM_Memory( 1:nDOFel )

        call ElementList(e)%El%GetGlobalMapping(AnalysisSettings,GM)

        call ElementList(e)%El%ElementInternalForce(AnalysisSettings, Fe, Status)

        if (Status%Error) then

            inquire(file='ErrorElement.txt',exist=ElementError_exists) !Check if file with failed elements already exists
            
            if (ElementError_exists) then
                open(37,file='ErrorElement.txt',status='old',access='append')
            else
                open(37,file='ErrorElement.txt',status='new')
            endif
                        
            write(37,'(a)') Trim(Status%ErrorDescription)
            write(37,'(a,i6)') 'Failed in element ',e
            write(37,'(a)') 'Nodal coordinates (X Y Z):'
            
            nNodes = ElementList(e)%El%GetNumberOfNodes()
            
            do i=1,nNodes
                write(37,'(3(1X,f12.6))') (ElementList(e)%El%ElementNodes(i)%Node%Coord(j), j=1,3)
            enddo
            
            write (37,*) ''   
            
            close(37)
            return
            
        endif

        Fint(GM) = Fint(GM) + Fe

    enddo

    !************************************************************************************

end subroutine

