!------------------------------------------------------------------------------
! Projet de Calcul Hautes Performance - ENSEIRB-MATMECA - Mai 2016
!------------------------------------------------------------------------------
!
! MODULE: Module modcharge
!
!> @author
!> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
!
! DESCRIPTION: 
!> Ce module contient la subroutine de calcul de répartition de la charge pour chaque processus.
!
! REVISION HISTORY:
! 14 Mai 2016 - Initial Version
!------------------------------------------------------------------------------


module modcharge
  implicit none
contains
!---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Pour le nombre d'elements total d'un tableau, un nombre maximun de processus, un numero de processus renvoie les indices de debut et de fin pour le processeur considéré.
  ! 
  !
  !> @param[in] n nombre d'elements total d'un tableau
  !> @param[in] Np nombre total de processus
  !> @param[in] me numero du processeur
  !> @param[out] i1 indice de debut 
  !> @param[out] i2 indice de fin
  !---------------------------------------------------------------------------  
  subroutine charge(n,Np,me,i1,i2)
    implicit none
    integer,intent(in)::n,Np,me
    integer,intent(out)::i1,i2
    integer::var,var2,var3

    var=n/Np
    var2=n-var*Np
    if ((me)<var2)then
       i1=var*me+me+1
       i2=var*me+me+var+1
    else
       i1=var*me+var2+1
       i2=var*me+var2+var
    end if
  end subroutine charge
  
end module modcharge
