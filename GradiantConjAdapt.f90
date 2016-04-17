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

    X=0.0d0
    R=matvecA(alpha,beta,gamma,Nx,X)-B
    d=R
    n=0
    
    do while (n<size(B) .and. norme2(R)>epsilon)
       w=matvecA(alpha,beta,gamma,Nx,d)
       varalpha=prodscal(d,R)/prodscal(d,w)
       X=X-varalpha*d
       R1=R-varalpha*w
       varbeta=(norme2(R1)**2)/(norme2(R)**2)
       R=R1
       d=R+varbeta*d
       n=n+1
    end do
    
    deallocate(R)
    deallocate(R1)
    deallocate(w)
    deallocate(d)
    
  end subroutine gradconjA


  function matvecA(alpha,beta,gamma,Nx,X)
    !compute A(matrix of problem)X
    implicit none
    !table of parameter(dt,Lx,Ly,Nx,Ny,D)
    real(wp),intent(in)::alpha,beta,gamma 
    real(wp),intent(in),dimension(:)::X
    real(wp),dimension(size(X))::matvecA 
    integer,intent(in)::Nx
    integer::i,j
    real(wp),dimension(:,:),allocatable::A

    matvecA=0.0d0


    do i=1,Nx
       if (mod(i,Nx)==1) then
          matvecA(i)=alpha*X(i)+beta*X(i+1)+gamma*X(i+Nx)
       else if (mod(i,Nx)==0)  then
          matvecA(i)=beta*X(i-1)+alpha*X(i)+gamma*X(i+Nx)
       else
          matvecA(i)=alpha*X(i)+beta*X(i-1)+beta*X(i+1)+gamma*X(i+Nx)
       end if
    end do
  
    do i=Nx+1,Nx*Nx-Nx
       if (mod(i,Nx)==1) then
          matvecA(i)=gamma*X(i-Nx)+alpha*X(i)+beta*X(i+1)+gamma*X(Nx+i)
       else if (mod(i,Nx)==0)  then
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+gamma*X(Nx+i)
       else
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+beta*X(i+1)+gamma*X(Nx+i)
       end if
    end do

    
    do i=Nx*Nx-Nx+1,Nx*Nx
       if (mod(i,Nx)==1) then
          matvecA(i)=gamma*X(i-Nx)+alpha*X(i)+beta*X(i+1)
       else if (mod(i,Nx)==0) then
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)
       else
          matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+beta*X(i+1)
       end if
    end do

  end function matvecA


  subroutine test_gradconjA()
    implicit none
    real(wp) ,dimension(:,:),allocatable::A
    real(wp),dimension(:),allocatable::X1,X2,Y1,Y2,B
    real(wp)::eps
    integer::n,i

    n=5

    eps=0.000001

    allocate(X1(n*n))
    allocate(Y1(n*n))
    allocate(X2(n*n))
    allocate(Y2(n*n))
    allocate(B(n*n))

    B=1.0d0

    print*,'---------Calcul classique avec  gradconj ------'
    call test_A(100,-2,-1,A,n)

    call  gradconj(A,X1,B,eps)

    Y1=matmul(A,X1)

    !print Bfinal
    print*,'B1final='
    do i=1,n*n
        print*,Y1(i)
    end do
    
    print*,'---------Calcul classique de  gradconjA------'

    call gradconjA(100.0d0,-2.0d0,-1.0d0,n,X2,B,eps)
    
    Y2=matmul(A,X2)
    
    !print Bfinal
    print*,'B2final='
    do i=1,n*n
       print*,Y1(i),Y2(i),Y1(i)-Y2(i)
    end do
  end subroutine test_gradconjA







  subroutine test_A(intalpha,intbeta,intgamma,M,N)
    !test
    real(wp),dimension(:,:),allocatable,intent(out)::M
    integer,intent(in)::n,intbeta,intalpha,intgamma
    integer,dimension(:,:),allocatable::A
    integer::i,j
    allocate(A(N*N,N*N))
    allocate(M(N*N,N*N))
    A=0.0d0

    
    do i=1,N
       if (mod(i,N)==1) then
          A(i,i)=intalpha
          A(i,i+1)=intbeta
          A(i,i+N)=intgamma
       else if (mod(i,N)==0)  then   
          A(i,i-1)=intbeta
          A(i,i)=intalpha
          A(i,i+N)=intgamma
       else
          A(i,i-1)=intbeta
          A(i,i)=intalpha
          A(i,i+1)=intbeta
          A(i,i+N)=intgamma
       end if
    end do
    
    do i=N+1,N*N-N
       if (mod(i,N)==1) then
          A(i,i-N)=intgamma
          A(i,i)=intalpha
          A(i,i+1)=intbeta
          A(i,N+i)=intgamma
       else if (mod(i,N)==0)  then
          A(i,i-N)=intgamma
          A(i,i-1)=intbeta
          A(i,i)=intalpha
          A(i,i+N)=intgamma
       else
          A(i,i-N)=intgamma
          A(i,i-1)=intbeta
          A(i,i)=intalpha
          A(i,i+1)=intbeta
          A(i,N+i)=intgamma
       end if
    end do
  
    
    do i=N*N-N+1,N*N
     if (mod(i,N)==1) then
        A(i,i-N)=intgamma
        A(i,i)=intalpha
        A(i,i+1)=intbeta
     else if (mod(i,N)==0) then
        A(i,i-N)=intgamma
        A(i,i-1)=intbeta
        A(i,i)=intalpha
     else
        A(i,i-N)=intgamma
        A(i,i-1)=intbeta
        A(i,i)=intalpha
        A(i,i+1)=intbeta
     end if
  end do
  
  do i=1,N*N
     do j=1,N*N
        write(*,'(I4)', ADVANCE='NO') A(i,j) 
     end do
     print*,' '
  end do

  M=A
end subroutine test_A

subroutine test_matvecA()
  implicit none
  real*8,dimension(:),allocatable::F1,G1,F2,G2
  real*8,dimension(:,:),allocatable::A
  integer::i,n
  
  n=5
  
  allocate (F1(n*n))
  allocate (G1(n*n))
  allocate (F2(n*n))
  allocate (G2(n*n))

  print*,"------- calcul matvecA -------"
  
  F1=0.000

  do i=1,n*n
     F1(i)=(i*0.1)**2
  end do

  G1=matvecA(100.0d0,-2.0d0,-1.0d0,n,F1)
  print*,"X="
  do i=1,n*n
     print*,F1(i)
  end do
    print*,"Y="
  do i=1,n*n
     print*,G1(i)
  end do
  
  print*,"------- calcul classique -------"
  print*,"A="
  call test_a(100,-2,-1,A,n)
  
  F2=F1

  G2=matmul(A,F2)
  
  print*,"X="
  do i=1,n*n
     print*,F1(i),F2(i)
  end do
  
  print*,"Y="
  do i=1,n*n
     print*,G1(i),G2(i),G1(i)-G2(i)
  end do

  deallocate (F1)
  deallocate (G1)
  deallocate (F2)
  deallocate (G2)

end subroutine test_matvecA



end module gradconjadapt
