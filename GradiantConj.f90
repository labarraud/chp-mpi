module basicgradconj
  use mod_precision
  implicit none

  

contains
  subroutine gradconj(A,X,B,epsilon,Nl)
    implicit none
    real(wp),intent(inout),dimension(:)::X
    real(wp),intent(in),dimension(:)::B
    real(wp),intent(in),dimension(:,:)::A
    real(wp),intent(in)::epsilon
    real(wp)::alpha,beta
    integer,intent(in)::Nl
    real(wp),dimension(:),allocatable::R,R1,w,d
    integer::n
    
    allocate(R(size(B)))
    allocate(R1(size(B)))
    allocate(w(size(B)))
    allocate(d(size(B)))

    X=0
    R=matvec(A,X)-B
    print*,matvec(A,X)
    d=R
    n=0
    
    do while (n<Nl .and. norme2(R)>epsilon)
       w=matvec(A,d)
       print*,matvec(A,d)
       alpha=prodscal(d,R)/prodscal(d,w)
       X=X-alpha*d
       R1=R-alpha*w
       beta=(norme2(R1)**2)/(norme2(R)**2)
       R=R1
       d=R+beta*d
       n=n+1
       print*,X
    end do
    
    deallocate(R)
    deallocate(R1)
    deallocate(w)
    deallocate(d)
    
  end subroutine gradconj

  function norme2(V)
    !compute norme of V
    implicit none
    real*8,intent(in),dimension(:)::V
    real*8::norme2
    integer::i
    norme2=0;
    do i=1,size(V)
       norme2=norme2+V(i)*V(i)
    end do
    norme2=sqrt(norme2)
  end function norme2

  function matvec(A,X)
    !compute AX
    implicit none
    real*8,intent(in),dimension(:,:)::A
    real*8,intent(in),dimension(:)::X
    real*8,dimension(size(X))::matvec
    integer::i,j
    
    do i=1,size(X)
       matvec(i)=0
       do j=1,size(X)
          matvec(i)=matvec(i)+A(i,j)*X(j)
       end do
    end do
  end function matvec
  
  function prodscal(X,Y)
    !compute XY
    implicit none
    real*8,dimension(:)::X,Y
    real*8::prodscal
    integer::i
    prodscal=0
    do i=1,size(X)
       prodscal=prodscal+X(i)*Y(i)
    end do
  end function prodscal


end module basicgradconj
