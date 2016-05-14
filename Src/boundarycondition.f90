!------------------------------------------------------------------------------
! Projet de Calcul Hautes Performance - ENSEIRB-MATMECA - Mai 2016
!------------------------------------------------------------------------------
!
! MODULE: Module mod_boudarycondition
!
!> @author
!> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
!
! DESCRIPTION: 
!> Ce module contient la fonction de calcul du vecteur RHS, qui contient l'ensemble des conditions aux bords 
!
! REVISION HISTORY:
! 14 Mai 2016 - Initial Version
!------------------------------------------------------------------------------
module mod_boundarycondition
  implicit none


contains
  !---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Fonction qui calcule le vecteur RHS sur la bande i0:in
  ! 
  !
  !> @param[in] G vecteur contenant les conditions aux bords droit et gauche 
  !> @param[in] U vecteur U de l'iteration precedente
  !> @param[in] F vecteur F contenant les sources
  !> @param[in] H vecteur contenant les conditions aux bords haut et bas
  !> @param[in] Lx taille du domaine selon x
  !> @param[in] Ly taille du domaine selon y 
  !> @param[in] Nx nombre de mailles selon x 
  !> @param[in] Ny nombre de mailles selon y
  !> @param[in] D constante de diffusion 
  !> @param[in] dt pas de temps
  !> @param[in] Me numero du proceseur
  !> @param[in] i0 indice de debut du vecteur X à calculer 
  !> @param[in] iN indice de fin de vecteur X à calculer
  !
  !> @return vecteur RHS(i0:iN)
  !---------------------------------------------------------------------------  
  function createRHS(G,U,F,H,D,Lx,Ly,Nx,Ny,dt,i0,iN,Me)
    implicit none
    real*8,intent(in)::D,Lx,Ly,dt
    integer,intent(in)::Nx,Ny,i0,iN,Me
    real*8,intent(in),dimension(i0:iN)::F,U
    real*8,intent(in),dimension(1:2*Nx)::G
    real*8,intent(in),dimension(1:2*Ny)::H
    real*8,dimension(i0:iN)::createRHS
    real*8::dx,dy
    integer::i,j,k

    dx=Lx/(1+Nx)
    dy=Ly/(1+Ny)
    
    do i=i0,iN

       createRHS(i)=F(i)+U(i)/dt

       if (mod(i,Nx) == 0) then
          createRHS(i)=createRHS(i)+D*H((i/Nx)*2)/(dx*dx)
       else if (mod(i,Nx) == 1) then
          createRHS(i)=createRHS(i)+D*H(2*(i/Nx)+1)/(dx*dx)
       end if

       if ( i < Nx+1 ) then
          createRHS(i)=createRHS(i)+D*G(i)/(dy*dy)
       else if (i<Nx*Ny .and. i> Nx*(Ny-1)) then
          createRHS(i)=createRHS(i)+D*G(mod(i,Nx)+Nx)/(dy*dy)
       else if (i==Nx*Ny) then
          createRHS(i)=createRHS(i)+D*G(2*Nx)/(dy*dy)
       end if
    end do
  end function createRHS

end module mod_boundarycondition

