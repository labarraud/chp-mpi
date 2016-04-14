program main
  use basicgradconj
  use mod_precision
  use boundarycondition
  use gradconjadapt
  use writeX
  
  implicit none
  real*8,dimension(:),allocatable::F,G,H,U,RHS
  real*8::Lx,Ly,dt,dx,dy,D,x,y,alpha,beta,gamma,epsilon
  integer::Nx,Ny,i,j



  
  call test_boundary_conditions()
  
  call test_gradconjA()
  
  call test_writing()


  !------- PREMIER CAS --------
  epsilon=0.00001
  Lx=1.0
  Ly=1.0
  Nx=10
  Ny=Nx

  allocate(F(Nx*Ny))
  allocate(G(Nx*Ny))
  allocate(H(Nx*Ny))
  allocate(U(Nx*Ny))
  allocate(RHS(Nx*Ny))
  
  D=1.
  dx=Lx/(Nx+1)
  dy=dx

  !Create F
  do i=1,Nx
     do j=1,Ny    
        x=dx*i
        y=dy*j
        F((i-1)*Nx+j)=2*(x-x*x+y-y*y)
     end do
  end do

!!$  do i=1,Nx*Nx
!!$     print*,F(i)
!!$  end do
!!$  
  G=0.0d0
  H=0.0d0
  U=0.0d0
  dt=1
        
  alpha=(1./dt)+2*D*(1./(dx*dx)+1./(dy*dy))
  beta=-D/(dx*dx)
  gamma=-D/(dy*dy)

 
  do i=1,10
     RHS = createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt)
     call gradconjA(alpha,beta,gamma,Nx,U,RHS,epsilon)
  end do
  
  call writevec('cas1.dat','unknown',Nx,U,dx,dy)
  
  deallocate(F)
  deallocate(G)
  deallocate(H)
  deallocate(U)
  deallocate(RHS)

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  !------- DEUXIEME CAS --------
  !------------------------------------------------------------------
  !------------------------------------------------------------------
  epsilon=0.00001
  Lx=1.0
  Ly=1.0
  Nx=20
  Ny=Nx

  allocate(F(Nx*Ny))
  allocate(G(Nx*Ny))
  allocate(H(Nx*Ny))
  allocate(U(Nx*Ny))
  allocate(RHS(Nx*Ny))
  
  D=1.
  dx=Lx/(Nx+1)
  dy=dx

  !Create F
  do i=1,Nx
     do j=1,Ny    
        x=dx*i
        y=dy*j
        F((i-1)*Nx+j)=sin(x)+cos(y)
     end do
  end do


  G=0.0d0
  H=0.0d0
  !Create G
  do i=1,Nx
     do j=1,Ny    
        x=dx*i
        y=dy*j
        G((i-1)*Nx+j)=sin(x)+cos(y)
        H((i-1)*Nx+j)=sin(x)+cos(y)
     end do
  end do
  
  
!!$  do i=1,Nx*Nx
!!$     print*,F(i)
!!$  end do
!!$  
 
  U=0.0d0
  dt=1
        
  alpha=(1./dt)+2*D*(1./(dx*dx)+1./(dy*dy))
  beta=-D/(dx*dx)
  gamma=-D/(dy*dy)

 
  do i=1,10
     RHS = createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt)
     call gradconjA(alpha,beta,gamma,Nx,U,RHS,epsilon)
  end do
  
  call writevec('cas2.dat','unknown',Nx,U,dx,dy)
  
  deallocate(F)
  deallocate(G)
  deallocate(H)
  deallocate(U)
  deallocate(RHS)

end program main
