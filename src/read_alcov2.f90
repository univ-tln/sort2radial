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

subroutine read_alcov2(i_s,ref_mean,ref_max,params)

! Lecture et interpolation de alcov2 en fonction du nombre d'echantillons
! independants params%q_OBS. La determination de alcov2 a ete realisee
! auparavant par simulation Monte Carlo.
!
! Input:
!  - Obj_par: structure de parametres
!  - i_s: nombre de sources recherchees
!
! Output:
!  - ref_mean: reference de covariance moyenne
!  - ref_max: reference de covariance a 90!
!  - Obj_par: structure de parametres mise a jour
!
! Appel de fonctions:
!

! Initializations
! ---------------
use init_params_module

implicit none

! Parameters structure
type(parameters)::params
integer::i_s
real(kind=8)::ref_mean,ref_max
real(kind=8),dimension(7,3,2)::alcov2
integer::i_source,max_obs,max_source,ii_obs
real(kind=8)::i_obs,di_obs


! Test si le fichier a deja ete lu ou non
if(params%read_alcov2)then
     ! Chargement de alcov2
      open(321,file=params%alcov2_file,access='direct',recl=8*7*3*2)
      read(321,rec=1)alcov2
      close(321)
      
     ! Garde en memoire alcov2
     params%alcov2=alcov2

     ! Empeche la relecture du fichier
     params%read_alcov2=.False.
endif

! Determination des parametres
max_obs=7
max_source=3

! Le tableau params%alcov2 ne commence qu'a deux sources
! => calcul indice dans le tableau 
i_source=i_s-1
if(i_source>max_source)then
     i_source=max_source
endif

! Interpolation du tableau sur le nombre d'echantillons params%q_OBS
! Tableau fourni pour des dimensions espace = puissances de
! racine(2)approchees.
! Tableau commence à 4 echantillons = (racine(2)^4 => -3 pour numero de
! ligne
i_obs=(2.*(log10(real(params%q_OBS,kind=8))/log10(2.)))-3.

! Cadrage et interpolation
if(i_obs<1.)then
     i_obs=1
endif
if(i_obs>real(max_obs,kind=8))then
     i_obs=max_obs
endif

! Interpolation
ii_obs=floor(i_obs)
di_obs=i_obs-ii_obs
if(di_obs > 0.01)then
     ref_mean=((1.-di_obs)*params%alcov2(ii_obs,i_source,1))+(di_obs*params%alcov2(ii_obs+1,i_source,1))
     ref_max=((1.-di_obs)*params%alcov2(ii_obs,i_source,2))+(di_obs*params%alcov2(ii_obs+1,i_source,2))
else
     ! Cas valeur entiere
     ref_mean=params%alcov2(ii_obs,i_source,1)
     ref_max=params%alcov2(ii_obs,i_source,2)
endif

return
end subroutine read_alcov2
