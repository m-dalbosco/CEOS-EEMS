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
    type(ClassElementProfile)                             :: ElProfile

	! Internal variables
	! -----------------------------------------------------------------------------------
	integer :: i, nNodes, nGP, gp, nGPe, e

	!************************************************************************************


	!************************************************************************************
	! CONSTRUCT OF MATERIALS
	!************************************************************************************

	call Element%AllocateGaussPoints(nGP)

    if ( ( nGP <= 0 ) .or. (nGP > 1000) ) then
        call Error ("Error: Number of the Gauss Points <=0 or >1000")
    endif
    
    ! Allocate the constitutive model for all element's Gauss point
    ! -----------------------------------------------------------------------------------
	call AllocateConstitutiveModel( Material%ModelEnumerator , AnalysisSettings , nGP ,  Element%GaussPoints )
        
    ! Copy material properties from reference material (read in the settings file) to
    ! Gauss points
    ! -----------------------------------------------------------------------------------
	do gp = 1,nGP
        call Element%GaussPoints(gp)%CopyProperties(Material%Mat(1))
    enddo

    ! Construct the Constitutive Model
    ! -----------------------------------------------------------------------------------
    do gp=1,nGP

        allocate( Element%GaussPoints(gp)%Stress( AnalysisSettings%StressSize ) )

        Element%GaussPoints(gp)%Stress = 0.0d0

        call Element%GaussPoints(gp)%ConstitutiveModelDestructor()

        call Element%GaussPoints(gp)%ConstitutiveModelConstructor(AnalysisSettings)

    enddo
    
	!************************************************************************************

    !************************************************************************************
	! CONSTRUCT OF MATERIALS - EXTRA GAUSS POINTS
	!************************************************************************************

    !Call profile to check if element has reinforcement capabilities
    call Element%GetProfile(ElProfile)
    
    if (ElProfile%AcceptFiberReinforcement == .true.) then
        
        call Element%GetNumberOfExtraGaussPoints(e,nGPe)
        
        ! Allocate the constitutive model for extra Gauss point
        ! -----------------------------------------------------------------------------------
        call AllocateConstitutiveModel( Material%ModelEnumerator , AnalysisSettings , nGPe ,  Element%ExtraGaussPoints )
        
        ! Copy material properties from reference material (read in the settings file) to extra Gauss points
        ! --------------------------------------------------------------------------------------------------
        do gp = 1,nGPe
            call Element%ExtraGaussPoints(gp)%CopyProperties(Material%Mat(1))
        enddo
        
        ! Construct the Constitutive Model
        ! -----------------------------------------------------------------------------------
        do gp=1,nGPe
            allocate( Element%ExtraGaussPoints(gp)%Stress( AnalysisSettings%StressSize ) )
            Element%ExtraGaussPoints(gp)%Stress = 0.0d0
            call Element%ExtraGaussPoints(gp)%ConstitutiveModelDestructor()
            call Element%ExtraGaussPoints(gp)%ConstitutiveModelConstructor(AnalysisSettings)
        enddo
    
    endif

end subroutine
