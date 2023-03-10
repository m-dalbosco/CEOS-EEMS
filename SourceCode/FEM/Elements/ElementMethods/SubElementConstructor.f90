!##################################################################################################
! This routine constructs the element.
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
    
	! Internal variables
	! -----------------------------------------------------------------------------------
	integer :: i, nNodes

	!************************************************************************************


	!************************************************************************************
	! CONSTRUCT THE ELEMENT
	!************************************************************************************



	call AllocateNewElement( this , ElementType )

    this%ElementMaterialID = MaterialID

	nNodes = this%GetNumberOfNodes()

	allocate( this%ElementNodes(nNodes) )

	do i=1,nNodes
		this%ElementNodes(i)%Node => GlobalNodesList( ElementNodes(i) )
	enddo

	!************************************************************************************



end subroutine
