program main
  use basicgradconj
  use mod_precision
  
  implicit none
  real(wp),dimension(3,3)::A
  real(wp),dimension(3)::X,B
  
  A(1,:)=(/1,1,0/)  
  A(2,:)=(/1,1,1/)  
  A(3,:)=(/0,1,1/)
  
  B=(/1,1,4/)

print*,'A=',A
print*,'B=',B
  
call gradconj(A,X,B,0.0d01,100)

print*,'X=',X
print*,'AX=',matmul(A,X)
end program main
