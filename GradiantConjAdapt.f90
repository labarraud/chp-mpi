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
       print*,X
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
    integer::i,j,intbeta,intalpha,intgamma
    integer,dimension(:,:),allocatable::A
    logical::test

    test=.true.
    allocate(A(Nx*Nx,Nx*Nx))
    A=0
    intalpha=floor(alpha)
    intbeta=floor(beta)
    intgamma=floor(gamma)




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

    print*,matvecA
   
  end function matvecA


  subroutine test_gradconjA()
    implicit none
    real(wp) ,dimension(:,:),allocatable::A
    real(wp),dimension(:),allocatable::X1,X2,Y1,Y2,B
    real(wp)::eps
    integer::n,i

    n=3
    eps=0.000001

    allocate(X1(n*n))
    allocate(Y1(n*n))
    allocate(X2(n*n))
    allocate(Y2(n*n))
    allocate(B(n*n))

    B=2.0d0

    print*,'---------Calcul classique de A ------'
    call test_A(A,n)

    print*,A

    call  gradconj(A,X1,B,eps,n)

    Y1=matmul(A,X1)

    !print Bfinal
    print*,'B1final='
    do i=1,n*n
        print*,Y1(i)
    end do
    
    print*,'---------Calcul classique de Afonc ------'

    call gradconjA(10.0d0,-2.0d0,-1.0d0,3,X2,B,eps)
    
    Y2=matmul(A,X2)
    
    !print Bfinal
    print*,'B2final='
    do i=1,n*n
       print*,Y2(i)
    end do
    
    
  end subroutine test_gradconjA

  subroutine test_A(M,N)
    !test
    real(wp),dimension(:,:),allocatable,intent(out)::M
    integer,intent(in)::n
    integer,dimension(:,:),allocatable::A
    integer::i,j,intbeta,intalpha,intgamma
    allocate(A(N*N,N*N))
    allocate(M(N*N,N*N))
    A=0
    intalpha=10
    intbeta=-2
    intgamma=-1 
    
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
        write(*,'(I2)', ADVANCE='NO') A(i,j) 
     end do
     print*,' '
  end do

  M=A
end subroutine test_A

end module gradconjadapt
