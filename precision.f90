module mod_precision
  
  implicit none
  
  ! working precision
  integer, parameter :: wp = kind(1.0d0)
  
  ! pi
  real(wp), parameter :: pi = acos(-1.0_wp)
  
  type R2
     real(wp) :: x(8)
  end type R2
  
end module mod_precision
