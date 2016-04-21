! Copyright Actimar 2016

! This software called SORT2RADIAL is a computer program whose purpose is to 
! derive radial current velocities from input SORT or USORT Wera radar 
! data using the MUSIC algorithm. It is a fortran version of 
! the Matlab code developped by "Laboratoire de Sondages Electromagnétiques
! de l'Environnement Terrestre, UMR 6017" now integrated into the MIO
! (Mediterranean Institute of Oceanography).

! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 

! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.

! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security.

! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.


subroutine vp_noise(norm_vp_noise,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
use init_params_module

implicit none

! Parameters structure
type(parameters)::params
! Eigen values for noise
real(kind=8),dimension(params%NAT)::norm_vp_noise
! Intermediate variables
real(kind=8)::niveau_eff,facteur_cadrage
integer::n_chirp,n_ant,MT,n_voie,n_range,k_ant
integer::ierr,t_OBS,n_OBS
integer,dimension(params%NAT)::list_ant
real(kind=8),allocatable,dimension(:,:,:,:)::TABSIG,TAB1
complex(kind=8),allocatable,dimension(:,:,:)::TRAW,SORTshort,SORTout,TOBS_OUT
complex(kind=8),allocatable,dimension(:,:)::obs
complex(kind=8),allocatable,dimension(:,:)::gamma
real(kind=8),allocatable,dimension(:,:)::mgammar,mgammai
complex(kind=8),allocatable,dimension(:,:,:)::mat_inter
real(kind=8),allocatable,dimension(:,:,:)::TRAWr,TRAWi
real(kind=8),allocatable,dimension(:)::TourChirp,wb,TourOBS,wbtime
integer::i,i_chirp,i_ant,i_range,i_COUR,i_debut,i_OBS,i_SUB,i_fdopp
real(kind=8)::pi,chouilla_seuil,enersig
complex::j
integer::lensav,lenwrk,NPTS
real(kind=8),allocatable,dimension(:)::work
real(kind=8),allocatable,dimension(:)::wsave
real(kind=8),allocatable,dimension(:)::wr,wi,valpro
real(kind=8),allocatable,dimension(:,:)::xr,xi,vecpro,toute_valpro

! Initialization
pi = acos(0.0)*2.
j=cmplx(0,1)

! Parametres propres a cette fonction purement illustratif (surement lie au
! filtrage). Inutile pour l'estimation des vp normalisees.
niveau_eff=1.
facteur_cadrage=(niveau_eff/7.071)*32556

! Parametres illustratifs d'une acquisition radar fictive
! Nombre de chirps
n_chirp=1024
! Nombre d'antennes
n_ant  =24
! Nombre d'echantillons par chirp
MT     =64
! Nombre de voies (I et Q)
n_voie =2
! Nombre de cellules distance considerees
n_range=4
! Nombre d'antennes selectionnees
k_ant  =params%NAT
! Tableau sur le nombres d'antennes selectionnees
list_ant=params%listant

! Parametres bases sur le traitement actuel
t_OBS=params%reg_radial*params%NSUB
n_OBS=size(params%list_OBS_ok)

lensav=2*MT+int(log(real(MT,kind=8)))+4
lenwrk=2*MT
allocate(wsave(lensav),stat=ierr)
allocate(work(lenwrk),stat=ierr)
call zfft1i(MT,wsave,lensav,ierr)


! Initialisation de la matrice des signaux
allocate(TABSIG(n_voie,MT,n_ant,n_chirp),stat=ierr)
allocate(TAB1(n_voie,MT,n_ant,n_chirp),stat=ierr)
allocate(TRAW(MT,n_ant,n_chirp),stat=ierr)
allocate(TRAWr(MT,n_ant,n_chirp),stat=ierr)
allocate(TRAWi(MT,n_ant,n_chirp),stat=ierr)
allocate(SORTshort(n_range,n_ant,n_chirp),stat=ierr)
allocate(SORTout(n_chirp,n_ant,n_range),stat=ierr)

TAB1=0.

! Generates random number
call random_number(TAB1)
TABSIG=facteur_cadrage*(2.*TAB1-1.)

! Fenetre de ponderation de Blackman (pour eviter l'effet de Gibbs)
allocate(TourChirp(MT),stat=ierr)
allocate(wb(MT),stat=ierr)

TourChirp=2.*pi*(/(i, i=0,MT-1, 1)/)/MT
wb=0.35875-(0.48829*cos(TourChirp))+(0.14128*cos(2.*TourChirp))-(0.01168*cos(3.*TourChirp));

! Construction du vecteur complexe
TRAW = TABSIG(1, :, :, :)+j*TABSIG(2, :, :, :)

! Application de la fenetre
! Boucles sur les tailles du tableau
do i_chirp=1,n_chirp
     do i_ant=1,n_ant          
         TRAW(:,i_ant,i_chirp)=wb*TRAW(:,i_ant,i_chirp)
     enddo
enddo

! Calcul de la fft
!! Transformation de Fourier
! Boucles sur les tailles du tableau
do i_chirp=1,n_chirp
     do i_ant=1,n_ant 
         call zfft1f(MT,1,TRAW(:,i_ant,i_chirp),MT,wsave,lensav,work,lenwrk,ierr)         
     enddo
enddo

! Multiplication due to FFT program
TRAW=MT*TRAW

! Fabrication du tableau par distance selon sens du chirp (normalement)
SORTshort=TRAW(1:n_range,:, :)

!! On rebrasse les dimensions pour avoir le numero de chirp en tete
do i_chirp=1,n_chirp
     do i_ant=1,n_ant
         do i_range=1,n_range
             SORTout(i_chirp,i_ant,i_range)=SORTshort(i_range,i_ant,i_chirp)
         enddo
     enddo
enddo

! Fabrication des observations
NPTS=nint(n_chirp/real(params%vaca_div,kind=8))

! Initialisation
allocate(TOBS_OUT(NPTS,n_ant,t_OBS),stat=ierr)
allocate(TourOBS(NPTS),stat=ierr)
allocate(wbtime(NPTS),stat=ierr)

! Overlapping (chevauchement) sub aquisitions
do i_OBS=1,t_OBS
     i_COUR=floor((i_OBS-1)/real(params%NSUB,kind=8))+1
     i_COUR=modulo(i_COUR-1,params%reg_radial)+1
     i_SUB=modulo(i_OBS-1,params%NSUB)+1   
     i_debut=floor(((i_SUB-1)*NPTS)/real(params%Kchevauch,kind=8))
     TOBS_OUT(:,:,i_OBS)=SORTout(i_debut+(/(i, i=1,NPTS, 1)/),:,int(i_COUR))
enddo

! Fenetre de ponderation de Blackman (pour eviter l'effet de Gibbs)
TourOBS=2.*pi*(/(i, i=0,NPTS-1, 1)/)/NPTS
wbtime=0.35875-(0.48829*cos(TourOBS))+ (0.14128*cos(2*TourOBS))-(0.01168*cos(3*TourOBS))

! Application de la fenetre
! Boucles sur les tailles du tableau
do i_obs=1,t_OBS
     do i_ant=1,n_ant            
         TOBS_OUT(:,i_ant,i_obs)=wbtime*TOBS_OUT(:,i_ant,i_obs)
     enddo
enddo

deallocate(wsave)
deallocate(work)
lensav=2*NPTS+int(log(real(NPTS,kind=8)))+4
lenwrk=2*NPTS
allocate(wsave(lensav),stat=ierr)
allocate(work(lenwrk),stat=ierr)
call zfft1i(NPTS,wsave,lensav,ierr)

! Calcul de la fft
do i_obs=1,t_OBS
     do i_ant=1,n_ant 
         call zfft1f(NPTS,1,TOBS_OUT(:,i_ant,i_obs),NPTS,wsave,lensav,work,lenwrk,ierr)         
     enddo
enddo

! Multiplication due to FFT program
TOBS_OUT=NPTS*TOBS_OUT

! Initialisation de la matrice interspectrale
allocate(mat_inter(k_ant,k_ant,NPTS),stat=ierr)
allocate(obs(k_ant,size(params%list_OBS_ok)),stat=ierr)
allocate(gamma(k_ant,k_ant),stat=ierr)
allocate(mgammar(k_ant,k_ant),stat=ierr)
allocate(mgammai(k_ant,k_ant),stat=ierr)

mat_inter=0.
! Fabrication de gamma pour toutes les antennes
! Boucle sur le nombre de raies
do i_fdopp=1,NPTS
     ! Extraction des donnees qui nous interessent
     obs=TOBS_OUT(i_fdopp,list_ant,params%list_OBS_ok)

     ! Calcul effectif de gamma, la matrice interspectrale pour la frequence
     ! Initialisation pour eviter les problemes de convergence 
     ! on va mettre un petit qqchose sur la diagonale
     chouilla_seuil=0.001
     
     gamma=0.
     do i_ant=1,k_ant
         gamma(k_ant,k_ant)=chouilla_seuil
     enddo

     ! Construction de gamma
     gamma=gamma+matmul(obs,transpose(conjg(obs)))

     ! Normalisation pour garder des amplitudes comparables
     gamma=gamma/(k_ant*n_OBS)

     ! Construction de la matrice interspectrale
     mat_inter(:,:,i_fdopp)=gamma
enddo

! Construction de la matrice de toutes les valeurs propres
allocate(wr(k_ant),stat=ierr)
allocate(wi(k_ant),stat=ierr)
allocate(valpro(k_ant),stat=ierr)
allocate(xr(k_ant,k_ant),stat=ierr)
allocate(xi(k_ant,k_ant),stat=ierr)
allocate(vecpro(k_ant,k_ant),stat=ierr)
allocate(toute_valpro(k_ant,NPTS),stat=ierr)

! Boucle sur le nombre de raies
do  i_fdopp=1,NPTS
     ! Diagonalisation de mgamma
     mgammar=real(mat_inter(:,:,i_fdopp))
     mgammai=aimag(mat_inter(:,:,i_fdopp))
     
     call cg(k_ant,k_ant, mgammar, mgammai, wr,  wi, 1, xr, xi, ierr )
     
     valpro=wr+j*wi
     vecpro=xr+j*xi

     ! Mise en format colonne des valeurs propres
     valpro=abs(valpro)
 
     ! Energie totale sur la raie normalisee
     enersig=sum(valpro)

!     ! Colonne de toutes les valeurs propres de la plus grande a la plus
     ! petite puis normalisation
     valpro=valpro/enersig

     ! Mise en forme de la matrice
     toute_valpro(:,i_fdopp)=valpro
enddo

!! Moyenne des valeurs propres normalisees
norm_vp_noise=sum(toute_valpro,2)/NPTS

return
end subroutine vp_noise
