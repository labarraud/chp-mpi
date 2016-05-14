!------------------------------------------------------------------------------
! Projet de Calcul Hautes Performance - ENSEIRB-MATMECA - Mai 2016
!------------------------------------------------------------------------------
!
! MODULE: Module mod_read_write
!
!> @author
!> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
!
! DESCRIPTION: 
!> Ce module contient l'ensemble des subroutines pour lire et ecrire dans des fichiers texte externes.
!
! REVISION HISTORY:
! 14 Mai 2016 - Initial Version
!------------------------------------------------------------------------------


module mod_read_write
  implicit none

contains
  
  !---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Pour un nom de fichier, un vecteur, les deltas, les nombres de mailles, les indices de début et de fin du vecteur ecrit le vecteur dans le fichier  
  ! 
  !
  !> @param[in] name nom du fichier 
  !> @param[in] sta statut du fichier
  !> @param[in] Nx nombre de mailles selon x
  !> @param[in] U vecteur a ecrire
  !> @param[in] dx taille d'une maille selon x
  !> @param[in] dy taille d'une maille selon y
  !> @param[in] i0 indice de debut du vecteur
  !> @param[in] iN indice de fin de vecteur
  !---------------------------------------------------------------------------  
  subroutine writevec(name,sta,Nx,U,dx,dy,i0,iN)
    implicit none
    character(len=*),intent(in)::name,sta
    integer,intent(in)::Nx,i0,iN
    real*8,intent(in)::dx,dy
    real*8,dimension(i0:iN),intent(in)::U
    integer::ivar,i,j,imin,imax,jmin,jmax
    
    open(unit=10,file=name,status=sta)
    
    do ivar=i0,iN
       j=1+(ivar-1)/Nx
       i=1+mod(ivar-1,Nx)
       write(10,*) dx*i,dy*j, U(ivar)
    end do
    
    close(10)
    
    print*,'fichier créé: ',name
    
  end subroutine writevec






  !---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Lit les parametres d'entree de configuration des variables globales dans la racine de l'executable nommé param.txt 
  ! 
  !
  !> @param[out] epsilon pas d'arret du gradiant conjugué 
  !> @param[out] Lx taille du domaine selon x
  !> @param[out] Ly taille du domaine selon y 
  !> @param[out] Nx nombre de mailles selon x 
  !> @param[out] Ny nombre de mailles selon y
  !> @param[out] D constante de diffusion 
  !> @param[out] dt pas de temps
  !> @param[out] Nt nombre de pas de temps maximun
  !> @param[out] folder chemin de sortie des résultats
  !> @param[out] cas cas test choisi
  !---------------------------------------------------------------------------
  subroutine read_param(epsilon,Lx,Ly,Nx,Ny,D,dt,Nt,folder,cas)
    implicit none
    integer,intent(out)::Nx,Ny,cas,Nt
    real*8,intent(out)::epsilon,Lx,Ly,D,dt
    CHARACTER(len=1000) input
    character(len=:), allocatable::folder
    
    open(unit=10,file='param.txt',status='old')
    read(10,*) epsilon
    read(10,*) Lx
    read(10,*) Ly
    read(10,*) Nx
    read(10,*) Ny
    read(10,*) D
    read(10,*) dt
    read(10,*) Nt
    read(10,*) input
    ALLOCATE(character(len=LEN(TRIM(input))) :: folder)
    folder=trim(input)
    read(10,*) cas
    
    close(10)
  end subroutine read_param


  
  !---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Pour un numero de processus et un type de fichier de sortie renvoie un nom de fichier correspondant
  ! 
  !
  !> @param[in] Me numero du processus
  !> @param[in] type type de fichier de sortie (sol,stat)
  !> @param[out] name nom du fichier correspondant
  !---------------------------------------------------------------------------
  subroutine Rename(Me,name,type)
    implicit none
    integer,intent(in)::Me
    character*20,intent(out)::name
    character(len=*),intent(in)::type
    character*3::tn
    integer::i1,i2,i3,i4
    i1=Me/100
    i2=(Me-100*i1)/10
    i3=Me - 100*i1 - 10*i2
    tn=char(i1+48)//char(i2+48)//char(i3+48)
    name=type//tn//'.dat'
  end subroutine Rename
  



!---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Pour un cas, nombre total de processeurs et le chemin de sortie des fichier solution génère le script gnuplot dans la racine de l'executable pour ploter directement l'ensemble des solutions de chaque processeur
  ! 
  !> @param[in] folder chemin de sortie des résultat
  !> @param[in] Np nombre total de processeurs
  !> @param[in] cas cas test choisi
  !---------------------------------------------------------------------------
  subroutine creatscriptgnuplot(folder,Np,cas)
    implicit none
    character(len=*),intent(in)::folder
    integer,intent(in)::Np,cas
    character*20::name
    character(len=3)::type
    integer::ivar,i,j,imin,imax,jmin,jmax
    
    type='sol'

    open(unit=5,file='plotresult.gnu',status='unknown')
    write(5,*) "set terminal postscript eps enhanced color" 
    write(5,*) "set output 'plotresult.eps' "
    select case (cas)
    case(1)
       write(5,*) "set view 60,75"
    case(2)
       write(5,*) "set view 40,30"
    case(3)
       write(5,*) "set view 60,15"
    case default
    end select

    write(5,'(A10)',advance='no') "splot "
    do i=0,Np-1
       call Rename(i,name,type)
       write(5,'(A35)',advance='no') adjustr ("'"//folder//'/'//name)//"' ,"
    end do
    select case (cas)
    case(1)
       write(5,'(A25)',advance='no') "x*(1-x)*y*(1-y) "
    case(2)
       write(5,'(A15)',advance='no') "sin(x)+cos(y) "
    case(3)
       write(5,'(A15)',advance='no') " "
    case default
    end select

  end subroutine creatscriptgnuplot
  
end module mod_read_write
