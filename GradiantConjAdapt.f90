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
    
  end subroutine gradconjA


  function matvecA(alpha,beta,gamma,Nx,X)
    !compute A(matrix of problem)X
    implicit none
    !table of parameter(dt,Lx,Ly,Nx,Ny,D)
    real(wp),intent(in)::alpha,beta,gamma 
    real(wp),intent(in),dimension(:)::X
    real(wp),dimension(size(X))::matvec
    integer::i,j
    
    
    matvecA(1)=alpha*X(1)+beta*X(2)+gamma*X(Nx+1)
    do i=2,Nx-1
       matvecA(i)=beta*X(i-1)+alpha*X(i)+beta*X(i+1)+gamma*X(i+Nx+1)
    end do
    matvecA(Nx)=beta*X(Nx-1)+alpha*X(i)+beta*X(i+1)+gamma*X(i+Nx+1)
    

  
    do i=Nx+1,Nx+Nx
       if (mod(i,Nx)==1) then
          matvecA(i)=gamma*X(i-Nx)+alpha*X(Nx+1)+beta*X(Nx+2)+gamma*X(Nx+Nx+1)        
       elseif (mod(i,Nx)==0))
          matvecA(i)=gamma*X(i-Nx)
       else
          matvecA(i)=gamma*X(i-Nx)+beta*X(Nx+i-1)+alpha*X(i)+beta*X(i+1)+gamma*X(i+Nx+1)
       end if
    end do


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
