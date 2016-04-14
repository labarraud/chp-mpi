program main
  use basicgradconj
  use mod_precision
  use boundarycondition
  use gradconjadapt
  use writeX


  implicit none

  call test_boundary_conditions()

  call test_matvecA()

  call test_writing()
  



end program main
