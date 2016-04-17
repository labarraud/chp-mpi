program main
  use basicgradconj
  use mod_precision
  use boundarycondition
  use gradconjadapt
  use writeX
  
  implicit none
  real*8,dimension(:),allocatable::F,G,H,U,RHS
  real*8::Lx,Ly,dt,dx,dy,D,x,y,alpha,beta,gamma,epsilon,t
  integer::Nx,Ny,i,j,n,nt

  !------------------------------------------------------------------!
  !-                                                                -!
  !-                        TEST UNITAIRE                           -!
  !-                                                                -!
  !------------------------------------------------------------------!

  !call test_matvecA()

  !call test_gradconjA()
  
  !call test_boundary_conditions()

  !call test_writing()



  !------------------------------------------------------------------!
  !-                                                                -!
  !-                        PREMIER CAS TEST                        -!
  !-                                                                -!
  !------------------------------------------------------------------!
  
  epsilon=0.00001
  Lx=1.0
  Ly=1.0
  Nx=10
  Ny=Nx

  allocate(F(Nx*Ny))
  allocate(G(2*Nx))
  allocate(H(2*Nx))
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


  
  !------------------------------------------------------------------!
  !-                                                                -!
  !-                       DEUXIEME CAS TEST                        -!
  !-                                                                -!
  !------------------------------------------------------------------!

  epsilon=0.000001
  Lx=1.0
  Ly=1.0
  Nx=10
  Ny=Nx

  allocate(F(Nx*Ny))
  allocate(G(2*Nx))
  allocate(H(2*Ny))
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
        F(i+(j-1)*Ny)=sin(x)+cos(y)
     end do
  end do

  G=0.0d0
  H=0.0d0

  !G bas
  do i=1,Nx
     x=dx*i
     y=dy*0
     G(i)=sin(x)+cos(y)
  end do


  !G haut
  do i=1,Nx
     x=dx*i
     y=dy*(Ny+1)
     G(i+Nx)=sin(x)+cos(y)
  end do

  !H gauche
  do i=0,Ny-1
     x=dx*0
     y=dy*(i+1)
     H(2*i+1)=sin(x)+cos(y)
  end do  

  !H droit
  do i=1,Ny
     x=dx*(Nx+1)
     y=dy*i
     H(2*i)=sin(x)+cos(y)
  end do  

  U=0.0d0
  dt=1
        
  alpha=+(1./dt)+2*D*(1./(dx*dx)+1./(dy*dy))
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


  !------------------------------------------------------------------!
  !-                                                                -!
  !-                        TROISIEME CAS TEST                      -!
  !-                                                                -!
  !------------------------------------------------------------------!
  
  epsilon=0.00001
  Lx=1.0
  Ly=1.0
  Nx=40
  Ny=Nx
  D=1.
  dx=Lx/(Nx+1)
  dy=dx

  allocate(F(Nx*Ny))
  allocate(G(2*Nx))
  allocate(H(2*Nx))
  allocate(U(Nx*Ny))
  allocate(RHS(Nx*Ny))

  G=0.0d0
  H=1.0d0
  U=0.0d0
  dt=0.1
        
  alpha=(1./dt)+2*D*(1./(dx*dx)+1./(dy*dy))
  beta=-D/(dx*dx)
  gamma=-D/(dy*dy)
  
  n=1
  do nt=1,5
     t=n*dt
     !create F
     do i=1,Nx
        do j=1,Ny    
           x=dx*i
           y=dy*j
           F((i-1)*Nx+j)=exp(-(x-(Lx/2))**2)*exp(-(y-(Ly/2))**2)*cos((pi/2)*t)
        end do
     end do
     
     RHS = createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt)
     call gradconjA(alpha,beta,gamma,Nx,U,RHS,epsilon)
     n=n+1
  end do
  
  call writevec('cas3_5.dat','unknown',Nx,U,dx,dy)
  
  print*,n*dt

  do nt=1,60
     t=n*dt
     !create F
     do i=1,Nx
        do j=1,Ny    
           x=dx*i
           y=dy*j
           F((i-1)*Nx+j)=exp(-(x-(Lx/2))**2)*exp(-(y-(Ly/2))**2)*cos((pi/2)*t)
        end do
     end do
     
     RHS = createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt)
     call gradconjA(alpha,beta,gamma,Nx,U,RHS,epsilon)
     n=n+1
  end do
  
  call writevec('cas3_10.dat','unknown',Nx,U,dx,dy)
  
  deallocate(F)
  deallocate(G)
  deallocate(H)
  deallocate(U)
  deallocate(RHS)


end program main
