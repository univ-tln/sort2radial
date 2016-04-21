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


module azim_module

! The aim of this module is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! -----------------
use init_params_module

implicit none

! Global variables
! ----------------
! Indices for files opening


contains


subroutine azim_proc(data_chevfft,var_ant,courant_dir,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
! Data after Doppler processing
complex(kind=8),dimension(params%NPTS,params%n_ant,params%n_OBS,params%n_range)::data_chevfft
! Antenna manifold
complex(kind=8),dimension(params%NAT,params%nazim)::var_ant
! Intermediate variables to load data
real(kind=8),dimension(params%prof,params%n_dir)::courant_dir,nrj_dir
real(kind=8),allocatable,dimension(:,:,:,:)::dopplertp_m
real(kind=8),allocatable,dimension(:,:,:)::mIcourant
real(kind=8),dimension(params%NPTS)::int_seuil
real(kind=8),dimension(params%NAT)::norm_vp_noise
! Intermediate variables
integer::i_couronne,n_porte,deltamax,ierr,ind,ii,jj,ind_pos,ind_neg,ind_all,i,ind_dir,ind_prof,n_nan,n_tot
logical::test_var
real(kind=8)::vmax,fmax,taux_hit
integer,allocatable,dimension(:)::liste_raies_vague_pos,liste_raies_vague_neg,liste_raies_vague  ! Intermediary table


! Estimation statistique de valeurs propres de bruit
call vp_noise(norm_vp_noise,params)

! Memory allocation
allocate(dopplertp_m(params%prof,params%n_dir,params%multiv2f,2),stat=ierr)
allocate(mIcourant(params%n_dir,params%multiv2f,3),stat=ierr)

courant_dir=999999999.

! Traitement DF par couronne
! !!!!!!!!!!!!!!!!!!!!!!!!!!
! Boucle sur les couronnes groupees
do i_couronne=1,params%prof

     ! Numero de porte
     n_porte=params%porte(i_couronne)
    
     ! Selection des raies candidates de la couronne selon la gamme de vitesse
     vmax=maxval(params%TVMAX(i_couronne,:))

     ! Definition selon vmax des bonnes raies pour recherche de sources 
     ! autour des raies de Bragg
     ! fmax en pts spectraux pour la vitesse max recherchee
     fmax=vmax*(2.*params%F/params%c_light)
    
     ! Indice entier correspondant
     deltamax=nint(fmax/params%DeltaF)
    
     ! On limite deltamax a Dibragg pour eviter l'ambiguite de vitesse
     if(deltamax>=params%DIbragg)then
         deltamax=params%DIbragg-1
     endif
    
     ! Liste des positions de raies utiles pour le signal
!     liste_raies_vague_pos=np.asarray(np.nonzero(abs(int_seuil-params%Ind_bragg[0])<=deltamax))
!     liste_raies_vague_neg=np.asarray(np.nonzero(abs(int_seuil-params%Ind_bragg[1])<=deltamax))
!     params%liste_raies_vague=np.union1d(liste_raies_vague_pos,liste_raies_vague_neg)
!
!     ! Il convient aussi de ne pas selectionner la zone raie zero
!     ! et la largeur peut dependre du site selon kzero 

     allocate(liste_raies_vague_pos(params%NPTS),stat=ierr)
     allocate(liste_raies_vague_neg(params%NPTS),stat=ierr)
     allocate(liste_raies_vague(params%NPTS),stat=ierr)
     
     int_seuil=(/(i, i=1,params%NPTS, 1)/)
     liste_raies_vague_pos=0
     liste_raies_vague_neg=0
     liste_raies_vague=0

     ind_pos=0
     ind_neg=0
     ind_all=0

     do ind=1,params%NPTS
         if(abs(int_seuil(ind)-params%Ind_bragg(1))<=deltamax)then
             ind_pos=ind_pos+1
             liste_raies_vague_pos(ind_pos)=ind
         endif
     enddo
     
     do ind=1,params%NPTS
         if(abs(int_seuil(ind)-params%Ind_bragg(2))<=deltamax)then
             ind_neg=ind_neg+1
             liste_raies_vague_neg(ind_neg)=ind
         endif
     enddo

     do ii=1,ind_pos
         test_var=.TRUE.
         do jj=1,ind_neg
             if(liste_raies_vague_pos(ii) .EQ. liste_raies_vague_neg(jj))then
                 test_var=.FALSE.
                 exit
             endif
         enddo
         if(test_var)then
             ind_all=ind_all+1
             liste_raies_vague(ind_all)=liste_raies_vague_pos(ii)
         endif
     enddo

     do jj=1,ind_neg
         do ii=1,ind_pos
         test_var=.TRUE.
             if(liste_raies_vague_pos(ii) .EQ. liste_raies_vague_neg(jj))then
                 test_var=.FALSE.
                 exit
             endif
         enddo
         if(test_var)then
             ind_all=ind_all+1
             liste_raies_vague(ind_all)=liste_raies_vague_neg(jj)
         endif
     enddo
          
     ind_pos=0
     ! The zero area of the largest width is not taken into account
     do ii=1,ind_all
         test_var=.TRUE.
         do jj=1,2*params%demi_zero+1
             if(liste_raies_vague(ii) .EQ. params%liste_raies_zero(jj))then
                 test_var=.FALSE.
                 exit
             endif
         enddo
         if(test_var)then
             ! Intermediate variable used for test purpose only
             ind_pos=ind_pos+1
             liste_raies_vague_pos(ind_pos)=liste_raies_vague(ii)
         endif
     enddo

     allocate(params%liste_raies_vague(ind_pos),stat=ierr)
     do ind =1,ind_pos
         params%liste_raies_vague(ind)=liste_raies_vague_pos(ind)
     enddo

     ! Generation du masque de support des vagues
     allocate( params%support_vague( params%NPTS),stat=ierr)
     params%support_vague=0
     params%support_vague(params%liste_raies_vague)=1
!
!     ! Selection des raies pour evaluation du seuil de bruit selon methode
!     ! hors du spectre de la mer: les bosses de Bragg et
!     ! on va prendre une garde pour ne pas travailler sur les bosses ordre 2
!     ! Liste des positions de raies utiles pour le bruit
!     liste_raies_bruit_pos=np.asarray(np.nonzero(abs(int_seuil-params%Ind_bragg[0])>delta2))
!     liste_raies_bruit_neg=np.asarray(np.nonzero(abs(int_seuil-params%Ind_bragg[1])>delta2))
!     liste_raies_bruit= np.intersect1d(np.squeeze(liste_raies_bruit_pos),np.squeeze(liste_raies_bruit_neg))
!
!     ! Il convient aussi de ne pas tenir compte non plus d un voisinage de zero
!     ! de plus grande largeur 

     ! Fonction de recherche des sources pour la couronne consideree
     call sourc_cour(data_chevfft,var_ant,norm_vp_noise,n_porte,dopplertp_m,params)

     ! Attribution de vitesse a toutes les sources
     ! Transformation doppler vers vitesse sur toutes les solutions,
     ! Sans recentrage du spectre, correction z, avec option egalisation DG
     ! Avec affaiblissement des solutions trop rapides par masque
     mIcourant=999999999.
     call multi_courant(dopplertp_m,i_couronne,mIcourant,params)

     ! Determination de la vitesse du courant, selon methode et params    
     ! puis elimination des directions trop faibles selon dB_dynadir
     ! Ponderation dans les cellules avec selection dB_sol pour methvit
     ! Determination d'une unique valeur de courant par cellule suivant le
     ! critere choisi
     call un_courant_music(mIcourant,i_couronne,courant_dir,nrj_dir,params)     
enddo

! Calcul du taux de couverture
n_nan=0
do ind_dir=1,params%n_dir
     do ind_prof=1,params%prof
         if(courant_dir(ind_prof,ind_dir)>100.)then
             n_nan=n_nan+1    
         endif
     enddo
enddo
             
n_tot=params%n_dir*params%prof

! taux_hit=1.-(n_nan/real(n_tot,kind=8))
! print*, 'Taux de couverture de la carte radiale:  ',taux_hit

params%not_nan=n_tot-n_nan

return
end subroutine azim_proc


end module azim_module
