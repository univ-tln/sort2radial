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


subroutine multi_courant(dopplertp_m,i_couronne,mIcourant,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
use init_params_module

implicit none

! Parameters structure
type(parameters)::params
real(kind=8),dimension(params%prof,params%n_dir,params%multiv2f,2)::dopplertp_m
integer::i_couronne
real(kind=8),dimension(params%n_dir,params%multiv2f,3)::mIcourant
real(kind=8),dimension(params%n_dir,params%multiv2f,2)::icouronne_courant
real(kind=8)::rapanis,invrapanis,VV
integer::ind_dir,ind_multi

! Rapport correspondant anisotropie 
rapanis=10.**(params%dB_anis/10.)
invrapanis=1./rapanis

! Initialisation
icouronne_courant=dopplertp_m(i_couronne,:,:,:)

! Calcul de la vitesse correspondant aux index
! Avec correction d altitude avec tableaux globaux
do ind_dir=1,params%n_dir
     do ind_multi=1,params%multiv2f
         if(icouronne_courant(ind_dir,ind_multi,1)<999999999.)then
             mIcourant(ind_dir,ind_multi,1)=params%vradz2D(i_couronne,&
                  int(icouronne_courant(ind_dir,ind_multi,1)))
         endif
     enddo 
enddo

! Stockage des energies
mIcourant(:,:,2)=icouronne_courant(:,:,2)

! Le tableau mIcourant garde la trace des indices origine
mIcourant(:,:,3)=icouronne_courant(:,:,1)

! Affaiblissement par l'energie des solutions trop veloces
! Attention! cette solution n'est pas bonne car l'apodisation spectrale dans
! un_courant etale la zone masquee a zero!
do ind_dir=1,params%n_dir
     VV=params%TVMAX(i_couronne,ind_dir)
     do ind_multi=1,params%multiv2f
         if(abs(mIcourant(ind_dir,ind_multi,1))>VV)then
             mIcourant(ind_dir,ind_multi,2)=1.
         endif
     enddo
enddo
 
return
end subroutine multi_courant
