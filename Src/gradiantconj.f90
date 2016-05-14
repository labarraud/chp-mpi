!------------------------------------------------------------------------------
! Projet de Calcul Hautes Performance - ENSEIRB-MATMECA - Mai 2016
!------------------------------------------------------------------------------
!
! MODULE: Module mod_gradconj
!
!> @author
!> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
!
! DESCRIPTION: 
!> Ce module contient la subroutine du gradiant conjugue parallelise ainsi que toutes les subroutines du produit matrice vecteur (Ax) et du calcul du produit scalaire
!
! REVISION HISTORY:
! 14 Mai 2016 - Initial Version
!------------------------------------------------------------------------------
module mod_gradconj
  use mpi
  implicit none

  
contains
  
  !---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> A partir de la fonction de calcul matrice vecteur(AX),du vecteur B, du numero du processeur, des indices de debut et de fin et du nombre total de processeurs, calcul le vecteur X du systeme (AX=B) par la methode du gradiant conjugue 
  ! 
  !
  !> @param[in] alpha variable alpha du probleme numerique
  !> @param[in] beta variable beta du probleme numerique
  !> @param[in] gamma variable gamma du probleme numerique
  !> @param[in] Nx nombre de mailles dans la direction x
  !> @param[out] X vecteur X que l'on cherche à determiner
  !> @param[in] B vecteur second membre B du systeme AX=B
  !> @param[in] Me numero du proceseur
  !> @param[in] i0 indice de debut du vecteur X à calculer 
  !> @param[in] iN indice de fin de vecteur X à calculer
  !> @param[in] Np nombre total de processeurs
  !---------------------------------------------------------------------------  
  subroutine gradconjA(alpha,beta,gamma,Nx,Ny,X,B,epsilon,Me,i0,iN,Np)
    implicit none
    integer,intent(in)::Nx,Ny,Me,Np,i0,iN
    real*8,intent(inout),dimension(max(i0-Nx,1):min(Nx*Ny,iN+Nx))::X
    real*8,intent(in),dimension(i0:iN)::B
    real*8,intent(in)::epsilon,alpha,beta,gamma
    real*8,dimension(max(i0-Nx,1):min(Nx*Ny,iN+Nx))::d
    real*8,dimension(i0:iN)::R,R1,w
    real*8::varalpha,varbeta,dr,dw,sumdr,sumdw,normR,normR1,sumnormR,sumnormR1
    integer::n,statinfo
    integer,dimension(MPI_STATUS_SIZE)::status
    
    X=0.0d0
    R(i0:iN)=matvecA(alpha,beta,gamma,Nx,Ny,X(max(i0-Nx,1):min(Nx*Ny,iN+Nx)),i0,iN)-B(i0:iN)

    !---------------------------------------------------------------------------
    ! calcul du carre de la norme de partie de vecteur R calculee par le processeur puis reduction avec sommation de toutes les parties calculees par les autres processeurs
    !---------------------------------------------------------------------------
    normR=prodscal(R(i0:iN),R(i0:iN),i0,iN)
    call MPI_ALLREDUCE(normR,sumnormR,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,statinfo)

    !---------------------------------------------------------------------------
    ! un fois le total du carre de la norme recu chaque processeur calcul sa racine 
    !--------------------------------------------------------------------------  
    sumnormR=sqrt(sumnormR)

    d(i0:iN)=R(i0:iN)
    
    n=0
    
    do while (n<Nx*Nx .and. sumnormR>epsilon)
       !---------------------------------------------------------------------------
       ! Si le processus n'est pas le premier, envoie de la partie haute du vecteur d (Nx premiers elements) au processeur precedent
       !---------------------------------------------------------------------------
       if (Me>0) then
          call MPI_SEND(d(i0:i0+Nx-1),Nx,MPI_REAL8,Me-1,n,MPI_COMM_WORLD,statinfo)
       end if

       !---------------------------------------------------------------------------
       ! Si le processeur n'est pas le dernier, envoie de la partie basse du vecteur d (Nx derniers elements) au processeur suivant 
       !---------------------------------------------------------------------------
       if (Me<Np-1) then
          call MPI_SEND(d(iN-Nx+1:iN),Nx,MPI_REAL8,Me+1,n,MPI_COMM_WORLD,statinfo)
       end if
       
       !---------------------------------------------------------------------------
       ! Si le processeur n'est pas le premier, reception de la partie basse du vecteur d du processeur precedent dans la partie extra haute du vecteur d du processeur
       !---------------------------------------------------------------------------
       if (Me>0) then
          call MPI_RECV(d(i0-Nx:i0-1),Nx,MPI_REAL8,Me-1,n,MPI_COMM_WORLD,status,statinfo)
       end if

       !---------------------------------------------------------------------------
       ! Si le processeur n'est pas le dernier, reception de la partie haute du vecteur d du processeur suivant dans la partie extra basse du vecteur d du processeur
       !---------------------------------------------------------------------------
       if (Me<Np-1) then
          call MPI_RECV(d(iN+1:iN+Nx),Nx,MPI_REAL8,Me+1,n,MPI_COMM_WORLD,status,statinfo)
       end if

       w(i0:iN)=matvecA(alpha,beta,gamma,Nx,Ny,d(max(i0-Nx,1):min(Nx*Nx,iN+Nx)),i0,iN)
       
       !---------------------------------------------------------------------------
       ! calcul du produit scalaire entre de partie de vecteur d et R calculee par le processeur puis reduction avec sommation de toutes les parties calculees par les autres processeurs
       !---------------------------------------------------------------------------
       dr=prodscal(d(i0:iN),R(i0:iN),i0,iN)
       call MPI_ALLREDUCE(dr,sumdr,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,statinfo)

       !---------------------------------------------------------------------------
       ! calcul du produit scalaire entre de partie de vecteur d et W calculee par le processeur puis reduction avec sommation de toutes les parties calculees par les autres processeurs
       !---------------------------------------------------------------------------
       dw=prodscal(d(i0:iN),w(i0:iN),i0,iN)
       call MPI_ALLREDUCE(dw,sumdw,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,statinfo)


       !---------------------------------------------------------------------------
       ! une fois le total de chaque produit scalaire recu calcul de varalpha par chaque processeur 
       !--------------------------------------------------------------------------         
       varalpha=sumdr/sumdw;
       X(i0:iN)=X(i0:iN)-varalpha*d(i0:iN)
       
       R1(i0:iN)=R(i0:iN)-varalpha*w(i0:iN)


       
       !---------------------------------------------------------------------------
       ! calcul du carre de la norme de partie de vecteur R calculee par le processeur puis reduction avec sommation de toutes les parties calculees par les autres processeurs
       !---------------------------------------------------------------------------
       normR=prodscal(R(i0:iN),R(i0:iN),i0,iN)
       call MPI_ALLREDUCE(normR,sumnormR,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,statinfo)

       !---------------------------------------------------------------------------
       ! calcul du carre de la norme de partie de vecteur R1 calculee par le processeur puis reduction avec sommation de toutes les parties calculees par les autres processeurs
       !---------------------------------------------------------------------------
       normR1=prodscal(R1(i0:iN),R1(i0:iN),i0,iN)
       call MPI_ALLREDUCE(normR1,sumnormR1,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,statinfo)   

       !---------------------------------------------------------------------------
       ! un fois le total de chaque carre de la norme recu chaque processeur calcul de varbeta
       !--------------------------------------------------------------------------  
       varbeta=sumnormR1/sumnormR

       !--------------------------------------------------------------------------
       ! un fois le total du carre de la norme recu chaque processeur calcul de la norme de R
       !--------------------------------------------------------------------------  
       sumnormR=sqrt(sumnormR1)
       R(i0:iN)=R1(i0:iN)
       
       d(i0:iN)=R(i0:iN)+varbeta*d(i0:iN)

       n=n+1
    end do
  end subroutine gradconjA


  !---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Fonction du produit matrice vecteur Ax en considérant que la Matrice est creuse et connue de manière itérative
  ! 
  !
  !> @param[in] alpha variable alpha du probleme numerique
  !> @param[in] beta variable beta du probleme numerique
  !> @param[in] gamma variable gamma du probleme numerique
  !> @param[in] Nx nombre de mailles dans la direction x
  !> @param[in] X vecteur X dont on veut faire le produit avec A, qui contient des extras de longeur Nx en partie haute et basse
  !> @param[in] i0 indice de debut du vecteur X à calculer 
  !> @param[in] iN indice de fin de vecteur X à calculer
  !
  !> @return resutat du produit AX sur la bande i0:iN
  !---------------------------------------------------------------------------  
  function matvecA(alpha,beta,gamma,Nx,Ny,X,i0,iN)
    !produit de la matrice A (du probleme) avec un vecteur X
    implicit none
    integer,intent(in)::Nx,Ny,i0,iN
    real*8,intent(in)::alpha,beta,gamma
    real*8,intent(in),dimension(max(i0-Nx,1):min(Nx*Ny,iN+Nx))::X
    real*8,dimension(i0:iN)::matvecA 
    integer::i,j,Statinfo
    real*8,dimension(:,:),allocatable::A

    matvecA=0.0d0
    

    if (i0<Nx) then 
    
       do i=max(i0,1),min(Nx,iN)
          if (mod(i,Nx)==1) then
             matvecA(i)=alpha*X(i)+beta*X(i+1)+gamma*X(i+Nx)
          else if (mod(i,Nx)==0)  then
             matvecA(i)=beta*X(i-1)+alpha*X(i)+gamma*X(i+Nx)
          else
             matvecA(i)=alpha*X(i)+beta*X(i-1)+beta*X(i+1)+gamma*X(i+Nx)
          end if
       end do
       
    end if

    if (iN>Nx+1) then   

       do i=max(Nx+1,i0),min(Nx*Ny-Nx,iN)
          if (mod(i,Nx)==1) then
             matvecA(i)=gamma*X(i-Nx)+alpha*X(i)+beta*X(i+1)+gamma*X(Nx+i)
          else if (mod(i,Nx)==0)  then
             matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+gamma*X(Nx+i)
          else
             matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+beta*X(i+1)+gamma*X(Nx+i)
          end if
       end do
       
    end if

    if (iN>Nx*Ny-Nx+1) then
    
       do i=max(Nx*Ny-Nx+1,i0),min(Nx*Ny,iN)
          if (mod(i,Nx)==1) then
             matvecA(i)=gamma*X(i-Nx)+alpha*X(i)+beta*X(i+1)
          else if (mod(i,Nx)==0) then
             matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)
          else
             matvecA(i)=gamma*X(i-Nx)+beta*X(i-1)+alpha*X(i)+beta*X(i+1)
          end if
       end do

    end if

    
  end function matvecA


!---------------------------------------------------------------------------  
  !> @author 
  !> Laurent Barraud et Clement Delibes (ENSEIRB-MATMECA)
  !
  ! DESCRIPTION: 
  !> Calcule le produit scalaire des vecteur X et Y pour des indices de debut i0 et fin iN
  ! 
  !
  !> @param[in] X vecteur X
  !> @param[in] Y vecteur Y
  !> @param[in] i0 indice de debut du vecteur X et Y 
  !> @param[in] iN indice de fin de vecteur X et Y
  !
  !> @return retourne le resultat du produit scalaire 
  !---------------------------------------------------------------------------    
  function prodscal(X,Y,i0,iN)
    !produit scalaire XY
    implicit none
    integer,intent(in)::i0,iN
    real*8,dimension(i0:iN),intent(in)::X,Y
    real*8::prodscal
    integer::i
    prodscal=0
    do i=i0,iN
       prodscal=prodscal+X(i)*Y(i)
    end do
  end function prodscal



  
end module mod_gradconj
