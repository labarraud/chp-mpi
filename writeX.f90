module writeX
  use boundarycondition


  implicit none

contains

  subroutine writevec(name,sta,Nx,X)
    implicit none
    character(len=*),intent(in)::name,sta
    integer,intent(in)::Nx
    real*8,dimension(Nx*Nx),intent(in)::X
    integer::i,j
    open(unit=10,file=name,status=sta)


    do i=1,Nx
       do j=1,Nx

          write(10,'(F18.6)',ADVANCE='NO') X( i+(j-1)*Nx )

       end do
       write(10,*)''
    end do
    close(10)

    print*,'fichier créé: ',name

  end subroutine writevec

  subroutine test_writing()
    implicit none
    real*8 ,dimension(:),allocatable::X,Y
    real*8,dimension(:),allocatable::F,G,H,U,var
    real*8::Lx,Ly,dt,D
    integer::Nx,Ny,i

    allocate(F(144))
    allocate(G(144))
    allocate(H(144))
    allocate(U(144))

    D=1.
    F=0.
    ! G=F
    ! H=F
    G=(1./12)*(1./12)
    H=(1./12)*(1./12)
    U=0.
    dt=12.
    Lx=1.0
    Ly=1.0
    Nx=12
    Ny=12

    allocate(var(Nx*Nx))

    var=createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt)

    call writevec('test.dat','unknown',Nx,var)

  end subroutine test_writing




end module writeX
