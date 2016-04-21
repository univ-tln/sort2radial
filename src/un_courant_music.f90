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


subroutine un_courant_music(mIcourant,i_couronne,courant_dir,nrj_dir,params)  

! Determination d'une unique valeur de courant par cellule suivant le
! critere choisi.
!
! Input:
!  - mIcourant matrice pour une couronne des vitesses en fonction de la
!  position angulaire
!    Taille: (params%n_dir,params%multiv2f,3)
!            (angle,nmax sources,1: vitesse courant - 2: nrj - 3: indice
!            raie )
!  - Obj_par: structure de parametres
! Output:
!  - courant_dir: matrice pour toutes les couronnes et toutes les positions
!  angulaires des vitesses trouvees par DF. Si pas de vitesse trouvee : NaN
!    Taille: (params%prof,params%n_dir)
!            (couronne,angle)
!
! Appel de fonctions:
!

! Initializations
! ---------------
use init_params_module

implicit none

! Parameters structure
type(parameters)::params
real(kind=8),dimension(params%n_dir,params%multiv2f,3)::mIcourant
integer::i_couronne
real(kind=8),dimension(params%prof,params%n_dir)::courant_dir,nrj_dir
real(kind=8),dimension(params%multiv2f,3)::courantsort
real(kind=8),dimension(params%multiv2f)::courant_vit,courant_nrj
integer,allocatable,dimension(:)::yacourant
integer,dimension(1)::ilocal
integer::idir,ind_ya,ierr,nbsolut,ind_multi
real(kind=8)::vlocal,nrjlocal,maxlocal,rap_couronne,seuil_nrjmini
real(kind=8)::ref_couronne,seuil_couronne

! Boucle sur les angles radar
do idir=1,params%n_dir
     ! Extraction des informations
     courantsort=mIcourant(idir,:,:)
     courant_vit=courantsort(:,1)
     courant_nrj=courantsort(:,2)
     
     ind_ya=0
     ! Test si au moins une source dans la cellule
     do ind_multi=1,params%multiv2f
         if(courant_vit(ind_multi)<99999999.)then
             ind_ya=ind_ya+1
         endif
     enddo
     
     allocate(yacourant(ind_ya),stat=ierr)
     nbsolut=ind_ya

     ind_ya=0
     ! Test si au moins une source dans la cellule
     do ind_multi=1,params%multiv2f
         if(courant_vit(ind_multi)<99999999.)then
             ind_ya=ind_ya+1
             yacourant(ind_ya)=ind_multi
         endif
     enddo
     
    ! Il y a au moins une solution dans la cellule!
    if (nbsolut>0)then
         ! Choix de la vitesse attribuée à la cellule en fonction de
         ! On prend la vitesse correspondant a la plus grande source
         ! (si plusieurs de meme valeur, ilocal est la premiere)
         maxlocal=maxval(courant_nrj(yacourant))
         ilocal=maxloc(courant_nrj(yacourant))
         vlocal=courant_vit(yacourant(ilocal(1)))
         nrjlocal=maxlocal      
         courant_dir(i_couronne,idir)=vlocal
         nrj_dir(i_couronne,idir)=nrjlocal
    else
         courant_dir(i_couronne,idir)=999999999.
         nrj_dir(i_couronne,idir)=999999999.
    endif
    
     deallocate(yacourant)

enddo

! Detection de terre selon la fluctuation dynamique dans la couronne
rap_couronne=10.**(params%dB_dynadir/10.)

! Detection de terre selon la fluctuation dynamique dans la couronne
rap_couronne=10.**(params%dB_dynadir/10.)

! Oui mais si il y a de la terre, la valeur moyenne n est pas indicateur
! il faut rester avec le max
ref_couronne=maxval(nrj_dir)
seuil_couronne=ref_couronne/rap_couronne
do idir=1,params%n_dir
     if(nrj_dir(i_couronne,idir)<seuil_couronne)then
         courant_dir(i_couronne,idir)=999999999.
     endif
enddo

seuil_nrjmini=10.**(params%dB_mini/10.)
do idir=1,params%n_dir
     if(nrj_dir(i_couronne,idir)<seuil_nrjmini)then
         courant_dir(i_couronne,idir)=999999999.
     endif
enddo



return
end subroutine un_courant_music
