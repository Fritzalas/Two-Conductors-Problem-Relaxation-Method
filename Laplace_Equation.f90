program Laplace_Equation
  use,intrinsic    :: iso_fortran_env ! iso_fortran_env: Ensures portability and consistency by providing standardized data type parameters and constants.With the use,intrinsic we prepare for the use of the module of intrincic that contains the ready function iso_fortran_env
  implicit none
  !--------------------------------------------------------------------------------
  !Declaration of variables:
  integer,parameter         :: dp=real64                   !We are going to use real64 data type for our data
  logical,allocatable       :: isConductor(:,:)            !A logical array to determine if we are in space or the limits of our given space
  real(dp),allocatable      :: V(:,:)                      !An array for the potencial V(x,y)
  real(dp)                  :: V1,V2                       !The values in the conductors
  real(dp)                  :: epsilon                     !The accuracy of our method
  integer                   :: L                           !Number of total places in our problem
  real(dp),parameter        :: ZERO=0.0_dp,QUARTER=0.25_dp !Useful constants for our calculations
  real(dp)                  :: C                           !The capacity of our system
  real(dp)                  :: Q1,Q2                       !The charges in each conductor
  real(dp),allocatable      :: sigma1(:,:),sigma2(:,:)     !Useful matrices for the conductors surface charge density 
  !---------------------------------------------------------------------------------
  !User Interface:
  print *, '#Enter V1,V2:'
  read  *, V1,V2
  print *, '#Enter the accuracy epsilon:'
  read  *, epsilon
  print *, '#Enter the number L of places in our lattice:'
  read  *, L
  print *, '#V1,V2= ',V1,V2
  print *, '#The relaxation method is having accuracy: ', epsilon
  print *, '#Number of places in lattice: ',L
  !---------------------------------------------------------------------------------
  !Array Allocation:
  ALLOCATE(isConductor(L,L))
  ALLOCATE(V(L,L))
  ALLOCATE(sigma1(L,L))
  ALLOCATE(sigma2(L,L))
  !---------------------------------------------------------------------------------
  !Calculations:
  call initialize_lattice(V,isConductor,L,V1,V2)               !Initializes the lattice and the test case for our V before the Gauss-Seidel method
  call laplace(V,isConductor,L,epsilon)                        !It runs our relaxation method
  call calc_charge_capacity(V,L,sigma1,sigma2,Q1,Q2,V1,V2,C)   !A subroutine that calculates the charge and also the capacity in each conductor
  call print_results(V,L,Q1,Q2,C,V1,V2)                        !Prints the desired results
  print *, '#Q1= ',Q1
  print *, '#Q2=', Q2
  print *, '#C= ', C
  !Internal Procedure:
contains
  !=================================================================================
  subroutine initialize_lattice(V,isConductor,L,V1,V2)
    implicit none
    !-------------------------------------------------------------------------------
    !Declaration of variables:
    integer                 :: L
    logical,dimension(L,L)  :: isConductor
    real(dp),dimension(L,L) :: V
    real(dp)                :: V1,V2
    integer                 :: i,j    !For our loops two integers
    !-------------------------------------------------------------------------------
    !Calculations:
    !We taking that every element in potencial V(x,y) has a zero initial value and also that every place in lattice is not occupied:
    V(:,:)=ZERO              !We taking that every element in the matrix has zero value. Also because of contains internal procedure we don't need to declare the variable ZERO again.
    isConductor(:,:)=.FALSE. !In most places in our lattice there is not a conductor
    !-------------------------------------------------------------------------------
    !First we have the outer conductor:
    do i=1,L
       V(1,i)=V1
       isConductor(1,i)=.TRUE.
       V(i,1)=V1
       isConductor(i,1)=.TRUE.
       V(i,L)=V1
       isConductor(i,L)=.TRUE.
       V(L,i)=V1
       isConductor(L,i)=.TRUE.
    end do
    !-------------------------------------------------------------------------------
    !Next we have tre inner conductor:
    do i=(2*L/5),(3*L/5)
       V(2*L/5,i)=V2
       isConductor(2*L/5,i)=.TRUE.
       V(3*L/5,i)=V2
       isConductor(3*L/5,i)=.TRUE.
       V(i,2*L/5)=V2
       isConductor(i,2*L/5)=.TRUE.
       V(i,3*L/5)=V2
       isConductor(i,3*L/5)=.TRUE.
    end do
    !-------------------------------------------------------------------------------
  end subroutine initialize_lattice
  !=================================================================================
  subroutine laplace(V,isConductor,L,epsilon)
    implicit none
    !-------------------------------------------------------------------------------
    !Declaration of variables:
    integer                 :: L
    logical,dimension(L,L)  :: isConductor
    real(dp),dimension(L,L) :: V
    real(dp)                :: epsilon
    integer                 :: i,j          !For do loops because of arrays
    integer                 :: icount       !It counts how many sweeps we have did in our lattice
    real(dp)                :: Vav,dV,error !Vav it determines the average value of V for four neighbouring points,dV counts the difference between new and past sweep and error calculates tha maximum dV
    !-------------------------------------------------------------------------------
    !Calculations:
    icount=0  !We give an initial value for our counter
    !We making an infinite do loop that stops when we have the neccessary accuracy
    do while (.TRUE.)
       error=ZERO    !We don't have an error in every begging of our loop
       !We starting from place 2 to L-1 because we know the values for 1 and L and also we want four neighbours
       do i=2,L-1
          do j=2,L-1
             !Change potencial when we have empty space
             if (.NOT. isConductor(i,j)) then !The not makes a True to False and also the False to True. So it runs only for empty False points
                Vav=(V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1))*QUARTER !Multiplications are better from divisions in fortran
                dV=ABS(Vav-V(i,j)) !Find the difference with the previous value
                if (dV>error) then
                   error=dV !We take by the end of the sweep the maximum difference. So it stops the infinite loop when the max error is smaller than the accuracy
                end if
                V(i,j)=Vav
                !------------------------------------------------------------------
             end if
             !---------------------------------------------------------------------
          end do
          !-------------------------------------------------------------------------
       end do
       icount=icount+1 !We had one sweep
       print *, '#Number of try: ',icount, '#Biggest error is: ',error
       if (error<epsilon) then
          return !The same with exit
       end if
       !----------------------------------------------------------------------------
    end do
    !-------------------------------------------------------------------------------
  end subroutine laplace
  !=================================================================================
  subroutine print_results(V,L,Q1,Q2,C,V1,V2)
    implicit none
    !-------------------------------------------------------------------------------
    !Declarations of variables:
    integer                 :: L
    real(dp),dimension(L,L) :: V
    integer                 :: i,j   !For the do loops
    real(dp)                :: Q1,Q2,C
    real(dp)                :: V1,V2
    !-------------------------------------------------------------------------------
    !Open a file that we are going to put our data for V
    open(unit=11,file="laplace.dat")
    do i=1,L
       do j=1,L
          write(11,*) i,j,V(i,j) !We saving our data inside our file
       end do
       write(11,*) ''            !Empty line that is very useful for our calculations
       !----------------------------------------------------------------------------
    end do
    !Have a file with the charges of the conductors and the capacity of the system:
    open(unit=12,file="charge_capacity.dat")
    write(12,*) Q1,Q2,C,ABS(V1-V2)
    !Have a file with the capacity of the system of conductors and the length of our system:
    open(unit=13,file="C_L.dat")
    write(13,*) L,C
    close(13)
    close(12)
    close(11)
    !-------------------------------------------------------------------------------
  end subroutine print_results
  !=================================================================================
  subroutine calc_charge_capacity(V,L,sigma1,sigma2,Q1,Q2,V1,V2,C)
    implicit none
    !-------------------------------------------------------------------------------
    !Declaration of variables:
    integer                  :: L
    real(dp),dimension(L,L)  :: V,sigma1,sigma2
    real(dp)                 :: C
    real(dp)                 :: V1,V2
    real(dp)                 :: Q1,Q2
    real(dp),parameter       :: ONE=1.0_dp, ZERO=0.0_dp
    real(dp),parameter       :: PI=atan2(ZERO,-ONE)
    real(dp),parameter       :: PI_4=4*PI
    integer                  :: i,j
    real(dp)                 :: dummy
    !-------------------------------------------------------------------------------
    !Calculations:
    !We set initial values for our calculations:
    sigma1(:,:)=ZERO !We set all values of the surface density equal to zero because we have zero charges to the free space
    sigma2(:,:)=ZERO
    !Calculate the surface density for the inner conductor and also the charge Q2 with the help of a dummy index:
    !-------------------------------------------------------------------------------
    Q2=0.0_dp
    do i=(2*L/5),(3*L/5)
       dummy=-(V(2*L/5,i)-V((2*L/5)-1,i))/(PI_4)
       Q2=Q2+dummy
       dummy=-(V(3*L/5,i)-V((3*L/5)+1,i))/(PI_4)
       Q2=Q2+dummy
       dummy=-(V(i,2*L/5)-V(i,(2*L/5)-1))/(PI_4)
       Q2=Q2+dummy
       dummy=-(V(i,3*L/5)-V(i,(3*L/5)+1))/(PI_4)
       Q2=Q2+dummy
    end do
    !-------------------------------------------------------------------------------
    !Calculate the surface density for the outter conductor and also the charge Q1 with the help of a dummy index:
    Q1=0.0_dp
    do i=1,L
       dummy=-(V(1,i)-V(2,i))/(PI_4)
       Q1=Q1+dummy
       dummy=-(V(i,1)-V(i,2))/(PI_4)
       Q1=Q1+dummy
       dummy=-(V(i,L)-V(i,L-1))/(PI_4)
       Q1=Q1+dummy
       dummy=-(V(L,i)-V(L-1,i))/(PI_4)
       Q1=Q1+dummy
    end do
    !-------------------------------------------------------------------------------
    !Calculation of the capacity:
    C=ABS(Q1/(V1-V2))
  end subroutine calc_charge_capacity
  !=================================================================================
end program Laplace_Equation
!===================================================================================
