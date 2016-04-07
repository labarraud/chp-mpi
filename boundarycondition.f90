module boundarycondition
  implicit none


contains

  function createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt)
    implicit none
    real*8,intent(in)::D,Lx,Ly,dt
    integer,intent(in)::Nx,Ny
    real*8,intent(in),dimension(:)::G,U,F,H
    real*8,dimension(Nx*Nx)::createRHS
    real*8::dx,dy
    integer::i,j

    dx=Lx/(1+Nx)
    dy=Ly/(1+Ny)
    j=1

    do i=1,size(F)


       createRHS(i)=F(i)+U(i)/dt

       if ( i-j == 0 ) then
          createRHS(i)=createRHS(i)-D*H(i)/(dy*dy)
          j=j+Nx-1
       end if

       if ( i < Nx+1 ) then
          createRHS(i)=createRHS(i)-D*G(i)/(dx*dx)
       end if

       if ( i > Nx*(Ny-1) ) then
          createRHS(i)=createRHS(i)-D*G(i)/(dx*dx)
       end if


    end do
  end function createRHS

  subroutine cree_tableau_pas(T,a,b,h)
    implicit none
    real*8 ,intent(in)::a,b,h
    real*8 ,dimension(:) ,allocatable,intent(inout)::T
    integer::n,i
    real*8::z
    !n=floor((b-a)/h)+1
    n=1
    z=a
    do while (z<b)
       z=z+h
       n=n+1
    end do
    allocate(T(n))
    T(1)=a
    do i=2,n 
       T(i)=T(i-1)+h
    end do
  end subroutine cree_tableau_pas

  subroutine test_boundary_conditions()
    implicit none
    real*8 ,dimension(:),allocatable::X,Y
    real*8,dimension(:),allocatable::F,G,H,U,var
    real*8::Lx,Ly,dt,D
    integer::Nx,Ny

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
  end subroutine test_boundary_conditions
end module boundarycondition

