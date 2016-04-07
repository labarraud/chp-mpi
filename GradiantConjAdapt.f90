module gradconjadapt
  use mod_precision
  use basicgradconj
  implicit none

  

contains
  subroutine gradconjA(P,X,B,epsilon)
    implicit none
    real(wp),intent(inout),dimension(:)::X
    real(wp),intent(in),dimension(:)::B
    real(wp),intent(in),dimension(:,:)::A
    real(wp),intent(in)::epsilon
    real(wp)::alpha,beta
    !integer,intent(in)::Nl
    real(wp),dimension(:),allocatable::R,R1,w,d
    integer::n
    
    allocate(R(size(B)))
    allocate(R1(size(B)))
    allocate(w(size(B)))
    allocate(d(size(B)))

    X=0
    R=matvect(A,X)-B
    d=R
    n=0
    
    do while (n<size(B) .and. norme2(R)>epsilon)
       w=matmul(A,d)
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


  function matvecA(P,X)
    !compute A(matrix of problem)X
    implicit none
    !table of parameter(dt,Lx,Ly,Nx,Ny,D)
    real*8,intent(in),dimension(4)::P 
    real*8,intent(in),dimension(:)::X
    real*8,dimension(size(X))::matvec
    integer::i,j
    real(wp)::alpha,beta,gamma,dx,dy

    dx=P(2)/(1+P(4))
    dy=P(3)/(1+P(5))

    alpha=1./P(1)+(2*P(4))/dx**2+(2*P(4))/dy**2
    beta=-P(4)/dx**2
    gamma=-P(4)/dy**2
    
    matvecA(mod(i,Nx))=alpha*X(mod(i,Nx))+beta*X(mod(i+2),Nx))+gamma*X(mod(i,Nx))
    !do from 2 to Nx*Ny
    
    do i=Nx+1,size(X)
       matvecA(i)=gamma*X(mod(i,Nx))+alpha*X(mod(i/Nx,Nx)*Nx+1  )+beta*X(mod(i+2),Nx))+gamma*X(mod(i,Nx))
       do j=mod(i+2,Nx),
          matvec(i)=matvec(i)+A(i,j)*X(j)
       end do
    end do

end function matvecA
  
end module gradconjadapt
