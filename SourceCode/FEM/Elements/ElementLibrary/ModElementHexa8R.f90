!##################################################################################################
! This module has the attributes and methods of the eight-node hexahedra linear element with
! full integration and fiber reinforcement
! ID -> Hexa8R = 410
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
module ElementHexa8R

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use Element

	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordHexa8R => null()
    real(8), pointer , dimension(:)   :: WeightHexa8R       => null()

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementHexa8R: Attributes and methods of the element Hexa8R
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElement) :: ClassElementHexa8R

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile          => GetProfile_Hexa8R
            procedure :: GetGaussPoints      => GetGaussPoints_Hexa8R
            procedure :: GetNumberOfNodes    => GetNumberOfNodes_Hexa8R
            procedure :: GetShapeFunctions   => GetShapeFunctions_Hexa8R
            procedure :: GetDifShapeFunctions=> GetDifShapeFunctions_Hexa8R
            procedure :: AllocateGaussPoints => AllocateGaussPointsParameters_Hexa8R

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        subroutine GetProfile_Hexa8R(this,Profile)
            class(ClassElementHexa8R)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile( 0                 , &
            NumberOfNodes = 8                              , &
            IsQuadratic = .false.                          , &
            GeometryType = GeometryTypes %  Hexahedra      , &
            FullIntegrationCapable = .true.                , &
            MeanDilatationCapable=.true. , &
            FiberReinforcementCapable=.true. , &
            ElementDimension = 3 )

        end subroutine

        !==========================================================================================
        ! Method GetGaussPoints_Hexa8R:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetGaussPoints_Hexa8R(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa8R) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO Hexa8R METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordHexa8R
            Weight       => WeightHexa8R

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetNumberOfNodes_Hexa8R:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Hexa8R(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa8R) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES - Hexa8R
		    !************************************************************************************

            nNodes=8

		    !************************************************************************************

        end function
        !==========================================================================================

        !==========================================================================================
        ! Method GetShapeFunctions_Hexa8R:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetShapeFunctions_Hexa8R(this , NaturalCoord , ShapeFunctions )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexa8R) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i
            real(8)             :: xi , eta , zeta , id(8,3)
            real(8) , parameter :: R1=1.0d0

		    !************************************************************************************--

            !************************************************************************************
            ! SHAPE FUNTIONS - Hexa8R
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)=[ -R1 , -R1 , -R1 ]
            id(2,:)=[  R1 , -R1 , -R1 ]
            id(3,:)=[  R1 ,  R1 , -R1 ]
            id(4,:)=[ -R1 ,  R1 , -R1 ]
            id(5,:)=[ -R1 , -R1 ,  R1 ]
            id(6,:)=[  R1 , -R1 ,  R1 ]
            id(7,:)=[  R1 ,  R1 ,  R1 ]
            id(8,:)=[ -R1 ,  R1 ,  R1 ]

            do i=1,8
                ShapeFunctions(i) = (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta) / 8.0d0
            enddo

      	    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetDifShapeFunctions_Hexa8R:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Hexa8R(this , NaturalCoord , DifShapeFunctions )

  			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexa8R) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i
            real(8)             :: xi , eta , zeta , id(8,3)
            real(8) , parameter :: R1=1.0d0 , R0=0.0d0

		    !************************************************************************************

      	    !************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE- Hexa8R
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)=[ -R1 , -R1 , -R1 ]
            id(2,:)=[  R1 , -R1 , -R1 ]
            id(3,:)=[  R1 ,  R1 , -R1 ]
            id(4,:)=[ -R1 ,  R1 , -R1 ]
            id(5,:)=[ -R1 , -R1 ,  R1 ]
            id(6,:)=[  R1 , -R1 ,  R1 ]
            id(7,:)=[  R1 ,  R1 ,  R1 ]
            id(8,:)=[ -R1 ,  R1 ,  R1 ]


            do i=1,8
                DifShapeFunctions(i,1) = id(i,1)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta) / 8.0d0
            enddo
            do i=1,8
                DifShapeFunctions(i,2) = id(i,2)*(R1 + id(i,1)*xi)*(R1 + id(i,3)*zeta) / 8.0d0
            enddo
            do i=1,8
                DifShapeFunctions(i,3) = id(i,3)*(R1 + id(i,1)*xi)*(R1 + id(i,2)*eta) / 8.0d0
            enddo
            
      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method AllocateGaussPointsParameters_Hexa8R: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Hexa8R(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa8R) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: x , id(10,3)
            real(8),parameter::R1=1.0d0

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS - Hexa8R
		    !************************************************************************************

            !Number of Gauss Points
            nGP=10

            if (associated(NaturalCoordHexa8R)) return
            allocate( NaturalCoordHexa8R(nGP,3) , WeightHexa8R(nGP) )

            x=1.0d0/dsqrt(3.0d0)

            id(1,:)=[ -R1 , -R1 , -R1 ]
            id(2,:)=[  R1 , -R1 , -R1 ]
            id(3,:)=[  R1 ,  R1 , -R1 ]
            id(4,:)=[ -R1 ,  R1 , -R1 ]
            id(5,:)=[ -R1 , -R1 ,  R1 ]
            id(6,:)=[  R1 , -R1 ,  R1 ]
            id(7,:)=[  R1 ,  R1 ,  R1 ]
            id(8,:)=[ -R1 ,  R1 ,  R1 ]
            id(9,:)=[  0.0d0 ,  0.0d0 ,  R1 ]
            id(10,:)=[ 0.0d0 ,  0.0d0 ,  -R1 ]

            NaturalCoordHexa8R=id*x

            WeightHexa8R=1.0d0

		    !************************************************************************************

            end subroutine
        !==========================================================================================



end module
