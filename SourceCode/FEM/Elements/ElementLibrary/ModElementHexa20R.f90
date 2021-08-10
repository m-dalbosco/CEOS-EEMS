!##################################################################################################
! This module has the attributes and methods of the 20-node hexahedra quadratic element with
! full integration
! ID -> Hexa20 = 420
!--------------------------------------------------------------------------------------------------
! Date: 2020/02
!
! Authors:  Jose Luis M. Thiesen
!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author: 
!##################################################################################################
module ElementHexa20R

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use Element

	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordHexa20R => null()
    real(8), pointer , dimension(:)   :: WeightHexa20R       => null()

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementHexa20: Attributes and methods of the element Hexa20
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElement) :: ClassElementHexa20R

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile          => GetProfile_Hexa20R
            procedure :: GetGaussPoints      => GetGaussPoints_Hexa20R
            procedure :: GetNumberOfNodes    => GetNumberOfNodes_Hexa20R
            procedure :: GetShapeFunctions   => GetShapeFunctions_Hexa20R
            procedure :: GetDifShapeFunctions=> GetDifShapeFunctions_Hexa20R
            procedure :: AllocateGaussPoints => AllocateGaussPointsParameters_Hexa20R

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains
             
        !==========================================================================================
        subroutine GetProfile_Hexa20R(this,Profile)
            class(ClassElementHexa20R)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile(421                 , &
            NumberOfNodes = 20                             , &
            IsQuadratic = .true.                           , &
            GeometryType = GeometryTypes %  Hexahedra      , &
            FullIntegrationCapable = .false.               , &
            ReducedIntegrationCapable = .true.            , &
            MeanDilatationCapable=.false. , &
            ElementDimension = 3 )

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method GetGaussPoints_Hexa20:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetGaussPoints_Hexa20R(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa20R) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO HEXA20R METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordHexa20R
            Weight       => WeightHexa20R

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetNumberOfNodes_Hexa20R:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Hexa20R(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa20R) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES - HEXA20R
		    !************************************************************************************

            nNodes=20

		    !************************************************************************************

        end function
        !==========================================================================================

        !==========================================================================================
        ! Method GetShapeFunctions_Hexa20R:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetShapeFunctions_Hexa20R(this , NaturalCoord , ShapeFunctions )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexa20R) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i, j
            real(8)             :: xi , eta , zeta , id(20,3)
            real(8) , parameter :: R2 = 2.0d0, R1=1.0d0, R0 = 0.0d0

		    !************************************************************************************--

            !************************************************************************************
            ! SHAPE FUNTIONS - HEXA20R
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)= [ -R1 , -R1 , -R1 ]
            id(2,:)= [  R1 , -R1 , -R1 ]
            id(3,:)= [  R1 ,  R1 , -R1 ]
            id(4,:)= [ -R1 ,  R1 , -R1 ]
            id(5,:)= [ -R1 , -R1 ,  R1 ]
            id(6,:)= [  R1 , -R1 ,  R1 ]
            id(7,:)= [  R1 ,  R1 ,  R1 ]
            id(8,:)= [ -R1 ,  R1 ,  R1 ]
            id(9,:)= [  R0 , -R1 , -R1 ]
            id(10,:)=[  R1 ,  R0 , -R1 ]
            id(11,:)=[  R0 ,  R1 , -R1 ]
            id(12,:)=[ -R1 ,  R0 , -R1 ]
            id(13,:)=[  R0 , -R1 ,  R1 ]
            id(14,:)=[  R1 ,  R0 ,  R1 ]
            id(15,:)=[  R0 ,  R1 ,  R1 ]
            id(16,:)=[ -R1 ,  R0 ,  R1 ]
            id(17,:)=[ -R1 , -R1 ,  R0 ]
            id(18,:)=[  R1 , -R1 ,  R0 ]
            id(19,:)=[  R1 ,  R1 ,  R0 ]
            id(20,:)=[ -R1 ,  R1 ,  R0 ]

            do i = 1,8
                ShapeFunctions(i) = (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)*(-R2 + id(i,1)*xi + id(i,2)*eta + id(i,3)*zeta) / 8.0d0
            enddo
            
            do i = 9,20
                if (id(i,1) == R0) then
                    ShapeFunctions(i) = (R1 - (xi**2))*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)
                elseif (id(i,2) == R0) then
                    ShapeFunctions(i) = (R1 - (eta**2))*(R1 + id(i,1)*xi)*(R1 + id(i,3)*zeta)
                elseif (id(i,3) == R0) then
                    ShapeFunctions(i) = (R1 - (zeta**2))*(R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)
                endif
                ShapeFunctions(i) = ShapeFunctions(i)/4.0d0    
            enddo
            

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetDifShapeFunctions_Hexa20R:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Hexa20R(this , NaturalCoord , DifShapeFunctions )

  			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexa20R) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i, j
            real(8)             :: xi , eta , zeta , id(20,3)
            real(8) , parameter :: R2 = 2.0d0, R1=1.0d0 , R0=0.0d0

		    !************************************************************************************

      	    !************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE- HEXA20
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)= [ -R1 , -R1 , -R1 ]
            id(2,:)= [  R1 , -R1 , -R1 ]
            id(3,:)= [  R1 ,  R1 , -R1 ]
            id(4,:)= [ -R1 ,  R1 , -R1 ]
            id(5,:)= [ -R1 , -R1 ,  R1 ]
            id(6,:)= [  R1 , -R1 ,  R1 ]
            id(7,:)= [  R1 ,  R1 ,  R1 ]
            id(8,:)= [ -R1 ,  R1 ,  R1 ]
            id(9,:)= [  R0 , -R1 , -R1 ]
            id(10,:)=[  R1 ,  R0 , -R1 ]
            id(11,:)=[  R0 ,  R1 , -R1 ]
            id(12,:)=[ -R1 ,  R0 , -R1 ]
            id(13,:)=[  R0 , -R1 ,  R1 ]
            id(14,:)=[  R1 ,  R0 ,  R1 ]
            id(15,:)=[  R0 ,  R1 ,  R1 ]
            id(16,:)=[ -R1 ,  R0 ,  R1 ]
            id(17,:)=[ -R1 , -R1 ,  R0 ]
            id(18,:)=[  R1 , -R1 ,  R0 ]
            id(19,:)=[  R1 ,  R1 ,  R0 ]
            id(20,:)=[ -R1 ,  R1 ,  R0 ]

            do i=1,8
                DifShapeFunctions(i,1) = id(i,1)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)*(id(i,1)*xi + &
                    id(i,2)*eta + id(i,3)*zeta - R2) + (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + &
                    id(i,3)*zeta)*id(i,1)
                DifShapeFunctions(i,2) = id(i,2)*(R1 + id(i,1)*xi)*(R1 + id(i,3)*zeta)*(id(i,1)*xi + &
                    id(i,2)*eta + id(i,3)*zeta - R2) + (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + &
                    id(i,3)*zeta)*id(i,2)
                DifShapeFunctions(i,3) = id(i,3)*(R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(id(i,1)*xi + &
                    id(i,2)*eta + id(i,3)*zeta - R2) + (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + &
                    id(i,3)*zeta)*id(i,3)
                DifShapeFunctions(i,1) = DifShapeFunctions(i,1)/8.0d0
                DifShapeFunctions(i,2) = DifShapeFunctions(i,2)/8.0d0
                DifShapeFunctions(i,3) = DifShapeFunctions(i,3)/8.0d0
            enddo
            
            do i = 9,20
                    if (id(i,1) == R0) then
                        DifShapeFunctions(i,1) = -R2*xi*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,2) = id(i,2)*(R1 - xi**2)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,3) = id(i,3)*(R1 - xi**2)*(R1 + id(i,2)*eta)
                    elseif (id(i,2) == R0) then
                        DifShapeFunctions(i,1) = id(i,1)*(R1 - eta**2)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,2) = -R2*eta*(R1 + id(i,1)*xi)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,3) = id(i,3)*(R1 - eta**2)*(R1 + id(i,1)*xi)
                    elseif (id(i,3) == R0) then
                        DifShapeFunctions(i,1) = id(i,1)*(R1 + id(i,2)*eta)*(R1 - zeta**2)
                        DifShapeFunctions(i,2) = id(i,2)*(R1 + id(i,1)*xi)*(R1 - zeta**2)
                        DifShapeFunctions(i,3) = -R2*zeta*(R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)
                    endif
                DifShapeFunctions(i,1) = DifShapeFunctions(i,1)/4.0d0    
                DifShapeFunctions(i,2) = DifShapeFunctions(i,2)/4.0d0    
                DifShapeFunctions(i,3) = DifShapeFunctions(i,3)/4.0d0    
            enddo
            
      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method AllocateGaussPointsParameters_Hexa20: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Hexa20R(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa20R) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: x , id(8,3)
            real(8),parameter::R1=1.0d0

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS - HEXA20R
		    !************************************************************************************

            !Number of Gauss Points
            nGP=8

            if (associated(NaturalCoordHexa20R)) return
            allocate( NaturalCoordHexa20R(nGP,3) , WeightHexa20R(nGP) )

            x=1.0d0/dsqrt(3.0d0)

            id(1,:)=[ -R1 , -R1 , -R1 ]
            id(2,:)=[  R1 , -R1 , -R1 ]
            id(3,:)=[  R1 ,  R1 , -R1 ]
            id(4,:)=[ -R1 ,  R1 , -R1 ]
            id(5,:)=[ -R1 , -R1 ,  R1 ]
            id(6,:)=[  R1 , -R1 ,  R1 ]
            id(7,:)=[  R1 ,  R1 ,  R1 ]
            id(8,:)=[ -R1 ,  R1 ,  R1 ]

            NaturalCoordHexa20R=id*x

            WeightHexa20R=1.0d0

		    !************************************************************************************

            end subroutine
        !==========================================================================================


end module
