program main
  use basicgradconj
  use mod_precision
  use boundarycondition
  use gradconjadapt

  implicit none

  call test_boundary_conditions()

  call test_matvecA()



end program main
