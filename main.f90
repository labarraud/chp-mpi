program main
  use basicgradconj
  use mod_precision
  use boundarycondition

  implicit none
  ! real(wp),dimension(3,3)::A
  ! real(wp),dimension(3)::X,B
  real*8 ,dimension(:),allocatable::X,Y
  real*8,dimension(:),allocatable::F,G,H,U,var
  real*8::Lx,Ly,dt,D
  integer::Nx,Ny
  
!   A(1,:)=(/1,1,0/)  
!   A(2,:)=(/1,1,1/)  
!   A(3,:)=(/0,1,1/)
  
!   B=(/1,1,4/)

! print*,'A=',A
! print*,'B=',B
  
! call gradconj(A,X,B,0.0d01,100)

! print*,'X=',X
! print*,'AX=',matmul(A,X)

call cree_tableau_pas(X,0.0d0,1.0d0,0.1d0)
call cree_tableau_pas(Y,0.0d0,1.0d0,0.1d0)

allocate(F(size(X)*size(Y)))
allocate(G(size(X)*size(Y)))
allocate(H(size(X)*size(Y)))
allocate(U(size(X)*size(Y)))

D=1.
F=0.
! G=F
! H=F
G=0.
H=1.
U=0.
dt=12.
Lx=1.0
Ly=1.0
Nx=12
Ny=12

allocate(var(Nx*Nx))

var=createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt)

! print*,'X=',X
! print*,'Y=',Y
! print*,'F=',F
! print*,'G=',G
! print*,'H=',H
! print*,'U=',U


print*,'RHS=',var

end program main
