module gradconjadapt
  use mod_precision
  use basicgradconj
  implicit none

  

contains
  subroutine gradconjA(alpha,beta,gamma,Nx,X,B,epsilon)
    implicit none
    real(wp),intent(inout),dimension(:)::X
    real(wp),intent(in),dimension(:)::B
    real(wp),intent(in)::epsilon,alpha,beta,gamma
    integer,intent(in)::Nx
    real(wp),dimension(:),allocatable::R,R1,w,d
    real(wp)::varalpha,varbeta
    integer::n
    
    allocate(R(size(B)))
    allocate(R1(size(B)))
    allocate(w(size(B)))
    allocate(d(size(B)))
    
    X=0
    R=matvecA(alpha,beta,gamma,Nx,X,.false.)-B
    d=R
    n=0
    
    do while (n<size(B) .and. norme2(R)>epsilon)
       w=matvecA(alpha,beta,gamma,Nx,d,.false.)
       varalpha=prodscal(d,R)/prodscal(d,w)
       X=X-alpha*d
       R1=R-varalpha*w
       varbeta=(norme2(R1)**2)/(norme2(R)**2)
       R=R1
       d=R+varbeta*d
       n=n+1
       print*,X
    end do
    
    deallocate(R)
    deallocate(R1)
    deallocate(w)
    deallocate(d)
    
  end subroutine gradconjA


  function matvecA(alpha,beta,gamma,Nx,X,test)
    !compute A(matrix of problem)X
    implicit none
    !table of parameter(dt,Lx,Ly,Nx,Ny,D)
    real(wp),intent(in)::alpha,beta,gamma 
    real(wp),intent(in),dimension(:)::X
    logical,intent(in)::test
    real(wp),dimension(size(X))::matvecA 
    integer,intent(in)::Nx
    !test
    integer,dimension(:,:),allocatable::A
    !test
    integer::i,j,intbeta,intalpha,intgamma
    
    allocate(A(Nx*Nx,Nx*Nx))
    A=0
    intalpha=1
    intbeta=2
    intgamma=3




    do i=1,Nx
       if (mod(i,Nx)==1) then
          matvecA(i)=alpha*X(i)+beta*X(i+1)+gamma*X(i+Nx)
          !TEST
          if (test) then
             A(i,i)=intalpha
             A(i,i+1)=intbeta
             A(i,i+Nx)=intgamma
          end if
       else if (mod(i,Nx)==0)  then
          matvecA(i)=beta*X(i-1)+alpha*X(i)+gamma*X(i+Nx)
          !TEST
          if (test) then          
             A(i,i-1)=intbeta
             A(i,i)=intalpha
             A(i,i+Nx)=intgamma
          end if
       else
          matvecA(i)=beta*X(i-1)+alpha*X(i)+beta*X(i+1)+gamma*X(i+Nx)
          !TEST
          if (test) then
             A(i,i-1)=intbeta
             A(i,i)=intalpha
             A(i,i+1)=intbeta
             A(i,i+Nx)=intgamma
          end if
       end if
    end do
  
    do i=Nx+1,Nx*Nx-Nx
       if (mod(i,Nx)==1) then
          matvecA(i)=gamma*X(i-Nx)+alpha*X(i)+beta*X(i+1)+gamma*X(Nx+i)
          !TEST
          if (test) then
             A(i,i-Nx)=intgamma
             A(i,i)=intalpha
             A(i,i+1)=intbeta
             A(i,Nx+i)=intgamma
          end if
       else if (mod(i,Nx)==0)  then
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+gamma*X(Nx+i)
          !TEST
          A(i,i-Nx)=intgamma
          A(i,i-1)=intbeta
          A(i,i)=intalpha
          A(i,i+Nx)=intgamma
       else
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+beta*X(i+1)+gamma*X(Nx+i)
          !TEST
          if (test) then
             A(i,i-Nx)=intgamma
             A(i,i-1)=intbeta
             A(i,i)=intalpha
             A(i,i+1)=intbeta
             A(i,Nx+i)=intgamma
          end if
       end if
    end do

    
    do i=Nx*Nx-Nx+1,Nx*Nx
       if (mod(i,Nx)==1) then
          matvecA(i)=gamma*X(i-Nx)+alpha*X(i)+beta*X(i-1)
          !TEST
          if (test) then
             A(i,i-Nx)=intgamma
             A(i,i)=intalpha
             A(i,i+1)=intbeta
          end if
       else if (mod(i,Nx)==0) then
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)
          !TEST
          if (test) then
             A(i,i-Nx)=intgamma
             A(i,i-1)=intbeta
             A(i,i)=intalpha
          end if
       else
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+beta*X(i+1)
          !TEST
          if (test) then
             A(i,i-Nx)=intgamma
             A(i,i-1)=intbeta
             A(i,i)=intalpha
             A(i,i+1)=intbeta
          end if
       end if
    end do
    
    !TEST
    if (test) then
       do i=1,Nx*Nx
          do j=1,Nx*Nx
             write(*,'(I2)', ADVANCE='NO') A(i,j) 
          end do
          print*,' ',X(i)
       end do
    end if
    
  end function matvecA



  subroutine test_matvecA()
    implicit none
    real(wp),dimension(:,:),allocatable::A
    real(wp),dimension(:),allocatable::X,Y
    
    allocate(A(3*3,3*3))
    allocate(X(3*3))
    
    X=1.0d0
    
    Y=matvecA(1.0d0,2.0d0,3.0d0,3,X,.true.)
    
    print*,Y
  end subroutine test_matvecA


end module gradconjadapt
