!##################################################################################################
! This module has the attributes and methods for the Hyperelastic material model.
!--------------------------------------------------------------------------------------------------
! Date: 2015/03
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module NeoHookeanFiberExp

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ConstitutiveModel
    use ModStatus
    use MathRoutines    
    implicit none


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! "NameOfTheMaterialModel"Properties: Material Properties
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type NeoHookeanFiberExpProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: Gm, Km, k1, k2

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassNeoHookeanFiberExp

        !Obs.: The ClassConstitutiveModel already has the variables:
        ! - stress (Voigt notation)
        ! - strain (Voigt notation)
        ! - F (Deformation Gradient) (3x3 Tensor Components)
        ! - Jbar - Mean Dilation variable related to Simo-Taylor-Pister Variational Approach


		! Class Attributes : Usually the internal variables
		!----------------------------------------------------------------------------------------

		! Class Attributes : Material Properties
		!----------------------------------------------------------------------------------------
        type (NeoHookeanFiberExpProperties), pointer :: Properties => null()


        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_NeoHookeanFiberExp
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_NeoHookeanFiberExp
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_NeoHookeanFiberExp
             procedure :: GetResult                    => GetResult_NeoHookeanFiberExp
             procedure :: SwitchConvergedState         => SwitchConvergedState_NeoHookeanFiberExp
             procedure :: CopyProperties               => CopyProperties_NeoHookeanFiberExp

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_3D: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassNeoHookeanFiberExp) :: ClassNeoHookeanFiberExp_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_NeoHookeanFiberExp_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_NeoHookeanFiberExp_3D

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelConstructor_NeoHookeanFiberExp(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanFiberExp) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

		    !************************************************************************************

 		    !************************************************************************************
            ! ALLOCATE THE STATE VARIABLES
		    !************************************************************************************
            
		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method ConstitutiveModelDestructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelDestructor_NeoHookeanFiberExp(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanFiberExp) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATE THE STATE VARIABLES
		    !************************************************************************************

		    !************************************************************************************

        end subroutine
        !==========================================================================================



        !==========================================================================================
        ! Method ReadMaterialParameters_"NameOfTheMaterialModel": Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadMaterialParameters_NeoHookeanFiberExp(this,DataFile)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use Parser

            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanFiberExp) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser) :: DataFile

            ! Internal variables
            ! ---------------------------------------------------------------------------------
		    character(len=100), dimension(4) :: ListOfOptions, ListOfValues
		    logical, dimension(2)            :: FoundOption
		    integer                          :: i
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! All constitutive models must allocate its own properties!
		    allocate (this%Properties)
		    !____________________________________________________________________________________

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************

            ! Inform how the properties are shown in the "Settings" file.
            !------------------------------------------------------------------------------------
            ListOfOptions=["Gm","Km","k1","k2"]
            !------------------------------------------------------------------------------------

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)
		    !____________________________________________________________________________________

            ! Set the material properties: this%Properties%"NameOfTheProperty"
            ! Obs.: ListOfValues index must match with ListOfOptions index
            !------------------------------------------------------------------------------------
            this%Properties%Gm = ListOfValues(1)
            this%Properties%Km = ListOfValues(2)
            this%Properties%k1 = ListOfValues(3)
            this%Properties%k2 = ListOfValues(4)
            !------------------------------------------------------------------------------------


            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method CopyProperties_"NameOfTheMaterialModel": Routine that associates the material
        ! parameters in the Gauss Points
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_NeoHookeanFiberExp(this,Reference)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanFiberExp) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            class(ClassConstitutiveModel) :: Reference

		    !************************************************************************************

            ! Change field: "class is ( Class"NameOfTheMaterialModel"Q1P0 )"
            !-----------------------------------------------------------------------------------
             select type ( Reference )

                 class is ( ClassNeoHookeanFiberExp )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error: Subroutine CopyProperties Neo-Hookean FiberExp"
            end select
            !-----------------------------------------------------------------------------------

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_3D: Routine that
        ! contains the algorithm employed to update the state variables.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_NeoHookeanFiberExp_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use MathRoutines

            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanFiberExp_3D) :: this
            type(ClassStatus) :: Status

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: F(3,3), C(3,3), Cinv(3,3),Ciso(3,3), I(3,3), S(3,3), Sfric(3,3), devSfric(3,3)
            real(8) :: mX(3), A(3,3)
            real(8) :: J, p, Km, Gm, k1, k2, I4

		    !************************************************************************************

            !___________________________________________________________________________________
            !______________________________    REMARK    _______________________________________
            !
            !  DUE TO THE UPDATED LAGRANGIAN FORMULATION, THE OUTPUT STRESS MUST BE THE
            !  CAUCHY STRESS IN VOIGT NOTATION.
            !___________________________________________________________________________________


            !************************************************************************************
            ! ALGORITHM THAT UPDATES STRESS AND STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            Km = this%Properties%Km
            Gm = this%Properties%Gm
            k1 = this%Properties%k1
            k2 = this%Properties%k2
            F  = this%F
            mX = this%AdditionalVariables%mX
            
            ! -----------------------------------------------------------------------------------
            ! Kinematic Variables - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0
            
            ! Jacobian
            J = det(F)
            
            ! Right-Cauchy Green Strain
            C = matmul(transpose(F),F)
            
            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)
            
            ! Isochoric part of the Right-Cauchy Green Strain
            Ciso = (J**(-2.0d0/3.0d0))*C
            
            !Material Structural Tensor
            A = Tensor_Product(mX,mX)

            !Fourth Invariant
            I4 = Tensor_Inner_Product(C,A)

            ! -----------------------------------------------------------------------------------

            if (norm(mX)==0) then !calculate matrix contribution

                ! Modified Second Piola-Kirchhoff Matrix Stress - Calculated in 3D Tensorial Format
                ! -----------------------------------------------------------------------------------

                ! Second Piola-Kirchhoff Frictional
                Sfric = Gm*I

                ! Hydrostatic Pressure
                !p = Km*( 1.0d0 - (1.0d0/J) )
                p = Km*( J - 1.0d0  )
                !p = (BulkModulus/2)*( J - (1.0d0/J) )

                ! Deviatoric part of the Second Piola-Kirchhoff Frictional
                devSfric = Sfric - (1.0d0/3.0d0)*Tensor_Inner_Product(Sfric,C)*Cinv

                ! Modified Second Piola-Kirchhoff Stress
                S = (J**(-2.0d0/3.0d0))*devSfric + J*p*Cinv
                ! -----------------------------------------------------------------------------------

                ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
                ! notation.
                ! -----------------------------------------------------------------------------------
                S = matmul(matmul(F,S),transpose(F))/J

                this%Stress = Convert_to_Voigt_3D_Sym( S )
                ! -----------------------------------------------------------------------------------
                
            else !calculate fiber contribution
                
                if (I4 .gt. 1.0d0) then
                
                    !Fiber Second Piola stress
                    S = 2.0d0*k1*(I4-1.0d0)*exp((k2*(I4-1.0d0)**2.0d0)-1)*A
                
                    S = matmul(matmul(F,S),transpose(F))/J
                    this%Stress = Convert_to_Voigt_3D_Sym( S )
                
                else
                    
                    S = 0
                    this%Stress = Convert_to_Voigt_3D_Sym( S )
                
                endif
                
            endif
            

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_NeoHookeanFiberExp_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use MathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanFiberExp_3D) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:),intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, Km, Gm, k1, k2, d2PSIvol_dJ2, I4
            real(8) :: F(3,3), C(3,3),Cinv(3,3), Ciso(3,3), Sfric(3,3), I(3,3)
            
            real(8) :: mX(3), A(3,3), Avoigt(6)

            real(8) :: CV(6), CinvV(6), CisoV(6), SfricV(6), devSfricV(6), SisoV(6)
            real(8) :: PmV(6,6) , PV(6,6), Diso(6,6), Dvol(6,6), Is(6,6), Daux(6,6)

            real(8) :: devSfric(3,3), Siso(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            Km = this%Properties%Km
            Gm = this%Properties%Gm
            k1 = this%Properties%k1
            k2 = this%Properties%k2
            F  = this%F
            mX = this%AdditionalVariables%mX
            ! -----------------------------------------------------------------------------------

            ! Quantities calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Kinematic Variables - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0
            
            ! Jacobian
            J = det(F)
            
            ! Right-Cauchy Green Strain
            C = matmul(transpose(F),F)
            
            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)
            
            ! Isochoric part of the Right-Cauchy Green Strain
            Ciso = (J**(-2.0d0/3.0d0))*C
            
            !Material Structural Tensor
            A = Tensor_Product(mX,mX)
            Avoigt = Convert_to_Voigt_3D_Sym(A)

            !Fourth Invariant
            I4 = Tensor_Inner_Product(C,A)
            
            if (norm(mX)==0) then !calculate matrix contribution
            
                ! Matrix tangent modulus
                ! -----------------------------------------------------------------------------------

                ! Second Piola-Kirchhoff Frictional
                Sfric = Gm*I

                ! Hydrostatic Pressure
                !p = Km*( 1.0d0 - (1.0d0/J) )
                p = Km*( J - 1.0d0  )
                !p = (Km/2)*( J - (1.0d0/J) )

                ! Derivative of Hydrostatic Pressure
                !d2PSIvol_dJ2 = Km/(J**2.0d0)
                d2PSIvol_dJ2 = Km
                !d2PSIvol_dJ2 = (Km/2)*( 1 + (1/(J**2.0d0)) )

                ! -----------------------------------------------------------------------------------
                ! The subsequent computations are made in Voigt notation
                ! -----------------------------------------------------------------------------------


                ! Material tangent modulus in referential configuration
                ! -----------------------------------------------------------------------------------

                ! Right-Cauchy Green Strain
                CV = Convert_to_Voigt(C)

                ! Inverse of Right-Cauchy Green Strain
                CinvV = Convert_to_Voigt(Cinv)

                ! Isochoric part of the Right-Cauchy Green Strain
                CisoV = Convert_to_Voigt(Ciso)

                ! Second Piola-Kirchhoff Frictional
                SfricV = Convert_to_Voigt(Sfric)

                ! Deviatoric part of the Second Piola-Kirchhoff Frictional
                devSfricV = SfricV - (1.0d0/3.0d0)*Inner_Product_Voigt(SfricV,CV)*CinvV

                ! Isochoric part of the Second Piola-Kirchhoff
                SisoV = (J**(-2.0d0/3.0d0))*devSfricV

                ! Modified Projection Operator
                PmV = Square_Voigt(CinvV,CinvV) - (1.0d0/3.0d0)*Ball_Voigt(CinvV,CinvV)


                ! Isochoric part of the material tangent modulus in referential configuration
                Diso = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Inner_Product_Voigt(SfricV,CV)*PmV - &
                        (2.0d0/3.0d0)*( Ball_Voigt(SisoV,CinvV) + Ball_Voigt(CinvV,SisoV) )

                ! Pressure component of the material tangent modulus in referential configuration
                Dvol  = J*( p + J*d2PSIvol_dJ2 )*Ball_Voigt(CinvV,CinvV) - 2.0d0*J*p*Square_Voigt(CinvV,CinvV)

                ! Material tangent modulus in referential configuration
                Daux = Diso + Dvol

                ! -----------------------------------------------------------------------------------
                
            else !calculate fiber contribution
                
                if (I4 .gt. 1.0d0) then
                
                    ! Fiber material tangent modulus
                    Daux = 4*k1*(1.0d0+(2.0d0*k2*((I4-1.0d0)**2.0d0)))*exp((k2*(I4-1.0d0)**2.0d0)-1)*Ball_Voigt(Avoigt,Avoigt)
                
                else
                    
                    Daux = 0
                    
                endif
                
            endif

                ! Push-Forward:
                ! Computation of the spatial tangent modulus
                ! -----------------------------------------------------------------------------------
                D = Push_Forward_Voigt(Daux,F)
                ! -----------------------------------------------------------------------------------

		    !************************************************************************************

        end subroutine
        !==========================================================================================



        !==========================================================================================
        ! Method SwitchConvergedState_"NameOfTheMaterialModel": Routine that save de converged state.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine SwitchConvergedState_NeoHookeanFiberExp(this)
            class(ClassNeoHookeanFiberExp) :: this
        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetResult_NeoHookeanFiberExp(this, ID , Name , Length , Variable , VariableType  )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use MathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanFiberExp) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer :: ID

            ! Output variables
            ! -----------------------------------------------------------------------------------
            character(len=*)            :: Name
            integer                     :: Length, VariableType
            real(8) , dimension(:)      :: Variable

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer, parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), I(3,3), e(3,3), eV(6)
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! Initializing variable name.
		    Name = ''
		    !____________________________________________________________________________________


            ! Template to Export Result to GiD
            !------------------------------------------------------------------------------------

            !case(0)
                ! Inform the number of results
                !Length = 3
            !case(1)
                !Name = 'Name of the Variable'
                !VariableType = 'Type of the Variable (Scalar,Vector,Tensor(in Voigt Notation))'
                !Length = 'Size of the Variable'
                !Variable = Result to be informed. Inform a Gauss Point result or compute a new
                !           variable.

            !------------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            select case (ID)
            
                case(0)
            
                    Length = 2
            
                case(1)
            
                    Name='Cauchy Stress'
                    VariableType=Tensor
                    Length=size(this%Stress)
                    Variable(1:Length) = this%Stress
            
                case (2)
            
                    Name='Almansi Strain'
                    VariableType = Tensor
                    Length=size(this%Stress)
                    !-------------------------------------------------------------
                    !Almansi Strain
                    !-------------------------------------------------------------
                    e = 0.50d0*( I - matmul(this%F, transpose(this%F) ))
                    eV = Convert_to_Voigt(e)
                    Variable(1:Length) = eV(1:Length)
                    !-------------------------------------------------------------
            
                case default
                    call Error("Error retrieving result :: GetResult")
            end select

		    !************************************************************************************

        end subroutine
        !==========================================================================================




    end module

