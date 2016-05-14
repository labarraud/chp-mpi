program main
  use mod_boundarycondition
  use mod_gradconj
  use mod_read_write
  use modcharge
  use mod_matconst
  
  implicit none
  real*8,dimension(:),allocatable::F,G,H,U,RHS
  real*8::Lx,Ly,dt,dx,dy,D,x,y,alpha,beta,gamma,epsilon,t,starttime,endtime
  character*20::nameresult,namestat
  character(len=:), allocatable::folderresult,folderstat
  integer::Nx,Ny,i,j,n,nt,statinfo,Me,Np,i1,i2,ivar,cas

  
  call MPI_INIT(statinfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Me,statinfo)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
 
  

  call read_param(epsilon,Lx,Ly,Nx,Ny,D,dt,Nt,folderresult,cas)

  
  dx=Lx/(Nx+1)
  dy=Ly/(Ny+1)

  
  call charge(Nx*Ny,Np,me,i1,i2)
  
  print*,Me,i1,i2



  allocate(F(i1:i2))
  allocate(G(2*Nx))
  allocate(H(2*Ny))
  allocate(U(max(i1-Nx,1):min(Nx*Ny,i2+Nx)))
  allocate(RHS(i1:i2))
    
  select case (cas)
     
  case (1)
     
     !------------------------------------------------------------------!
     !-                                                                -!
     !-                        PREMIER CAS TEST                        -!
     !-                                                                -!
     !------------------------------------------------------------------!
     
     !Creation de F
     do ivar=i1,i2
        j=1+(ivar-1)/Nx
        i=1+mod(ivar-1,Nx)
        x=dx*i
        y=dy*j
        F((j-1)*Nx+i)=2*(x-x*x+y-y*y)
     end do
     
     
     G=0.0d0
     H=0.0d0
     U=0.0d0
     
     alpha=(1./dt)+2*D*(1./(dx*dx)+1./(dy*dy))
     beta=-D/(dx*dx)
     gamma=-D/(dy*dy)
     
     
     call MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
     starttime = MPI_Wtime();
     do i=1,Nt
        RHS = createRHS(G,U(i1:i2),F(i1:i2),H,D,Lx,Ly,Nx,Ny,dt,i1,i2,Me)
        call gradconjA(alpha,beta,gamma,Nx,Ny,U(max(i1-Nx,1):min(Nx*Ny,i2+Nx)),RHS(i1:i2),epsilon,Me,i1,i2,Np)
     end do
       endtime= MPI_Wtime();

  case(2)
  
     !------------------------------------------------------------------!
     !-                                                                -!
     !-                       DEUXIEME CAS TEST                        -!
     !-                                                                -!
     !------------------------------------------------------------------!
     
     
     !Creation de F
     do ivar=i1,i2
        j=1+(ivar-1)/Nx
        i=1+mod(ivar-1,Nx)   
        x=dx*i
        y=dy*j
        F(i+(j-1)*Nx)=sin(x)+cos(y)
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
     
     alpha=+(1./dt)+2*D*(1./(dx*dx)+1./(dy*dy))
     beta=-D/(dx*dx)
     gamma=-D/(dy*dy)

     call MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
     starttime = MPI_Wtime();
     do i=1,Nt
        RHS = createRHS(G,U(i1:i2),F(i1:i2),H,D,Lx,Ly,Nx,Ny,dt,i1,i2,Me)
        call gradconjA(alpha,beta,gamma,Nx,Ny,U(max(i1-Nx,1):min(Nx*Ny,i2+Nx)),RHS(i1:i2),epsilon,Me,i1,i2,Np)
     end do
     
     endtime= MPI_Wtime();
      
  case(3)
     
     !------------------------------------------------------------------!
     !-                                                                -!
     !-                        TROISIEME CAS TEST                      -!
     !-                                                                -!
     !------------------------------------------------------------------!
     
     G=0.0d0
     H=1.0d0
     U=0.0d0
     
     alpha=(1./dt)+2*D*(1./(dx*dx)+1./(dy*dy))
     beta=-D/(dx*dx)
     gamma=-D/(dy*dy)
     
     n=1
     
     call MPI_BARRIER(MPI_COMM_WORLD,Statinfo)
     starttime = MPI_Wtime();

     do nt=1,Nt
        t=n*dt
        !creation de F
        do ivar=i1,i2
           j=1+(ivar-1)/Nx
           i=1+mod(ivar-1,Nx)   
           x=dx*i
           y=dy*j
           F((j-1)*Nx+i)=exp(-(x-(Lx/2))**2)*exp(-(y-(Ly/2))**2)*cos((pi/2)*t)
        end do
       
        
        RHS = createRHS(G,U(i1:i2),F(i1:i2),H,D,Lx,Ly,Nx,Ny,dt,i1,i2,Me)
        call gradconjA(alpha,beta,gamma,Nx,Ny,U(max(i1-Nx,1):min(Nx*Ny,i2+Nx)),RHS(i1:i2),epsilon,Me,i1,i2,Np)
        n=n+1
     end do
     endtime= MPI_Wtime();
     
  case default
     print*,'ERROR --- CAS NON VALIDE --- '
  end select


  
  
  !print du vecteur solution pour chaque proc
  call Rename(Me,nameresult,'sol')
  call writevec(folderresult//'/'//nameresult,'unknown',Nx,U(i1:i2),dx,dy,i1,i2)



  !print des statistiques d execution pour chaque proc 
  call Rename(Me,namestat,'stat')
  ALLOCATE(character(len=LEN(TRIM('Stat'))) :: folderstat)
  folderstat='Stat'
  folderstat=folderstat//'/'//namestat
  open(unit=5,file=folderstat,status='unknown')
  write(5,'(F20.15)') endtime-starttime
  close(5)

  !creation du script gnuplot par le processus 0
  call creatscriptgnuplot(folderresult,Np,cas)




  deallocate(F)
  deallocate(G)
  deallocate(H)
  deallocate(U)
  deallocate(RHS)
  deallocate(folderresult)
  deallocate(folderstat)

  call MPI_FINALIZE(statinfo)


end program main
