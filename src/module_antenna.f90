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


module antenna_module

! The aim of this module is to compute antenna coordinates and manifold, to apply some corrections, to 
! select the number of "good" antennas and to perform normalization

! Initializations
! -----------------
use init_params_module

implicit none


contains


subroutine proc_antenna(sort_data,var_ant,params)

! The aim of this function is to compute antenna coordinates and manifold, to apply some corrections, to 
! select the number of "good" antennas and to perform normalization

implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),dimension(params%n_chirp,params%n_ant,params%n_range)::sort_data
! Antenna coordinates
real(kind=8),dimension(3,params%n_ant)::coord_ant
! Antenna manifold
complex(kind=8),dimension(params%n_ant,params%nazim)::var_ant
! Antenna power
real(kind=8),dimension(params%n_ant)::ant_hpower,ant_hcorr
! Signal power of the first range cell for all antennas
real(kind=8),dimension(params%n_ant)::ant_DCpower
! Intermediate variables
integer::i,ierr


! Antenna coordinates
call pos_ant(coord_ant,params)

! Antenna manifold
call rep_ant(var_ant,coord_ant,params)

! HF antenna
call hf_ant(sort_data,ant_hpower,ant_hcorr,params)

! Antenna correction
call corr_ant(var_ant,ant_hcorr,params)

! Automatic selection of antennas or not
if(params%select_ant)then
     ! Signal power of the first range cell for all antennas
     call dc_ant(sort_data,ant_DCpower,params)

     ! Antenna selection
     call select_ant(ant_hpower,ant_DCpower,params)
else
     ! All antennas are kept
     params%NAT=params%n_ant
     
     ! Memory allocation
     allocate(params%listant(params%n_ant),stat=ierr)
     
     ! All antennas
     params%listant=(/(i, i=1,params%n_ant, 1)/)
endif
    
return
end subroutine proc_antenna


subroutine pos_ant(coord_ant,params)

! The aim of this function is to compute antenna coordinates

implicit none

! Parameters structure
type(parameters)::params
! Antenna coordinates
real(kind=8),dimension(3,params%n_ant)::coord_ant
! Intermediate variables
real(kind=8)::anglex,ranglex,angley,rangley
real(kind=8),dimension(3,3)::matchange
real(kind=8)::pi


! Initialization
pi = acos(0.0)*2.

! Reads antenna coordinates file
open(ind_file_pos_ant,file=trim(params%pos_ant_file))
read(ind_file_pos_ant,*)coord_ant
close(ind_file_pos_ant)

! Change coordinates
! X rotation angle (in degrees)
anglex   =(360.-params%cap_bous)+90.

! X rotation angle (in radians)
ranglex  =pi*(anglex/180.)

! Y rotation angle (in degrees)
angley   =anglex+90.

! Y rotation angle (in radians)
rangley  =pi*(angley/180.)

! Rotation matrix
matchange(1,1)=cos(ranglex)
matchange(1,2)=cos(rangley)
matchange(1,3)=0.
matchange(2,1)=sin(ranglex)
matchange(2,2)=sin(rangley)
matchange(2,3)=0.
matchange(3,1)=0.
matchange(3,2)=0.
matchange(3,3)=1.

! New antenna coordinates
coord_ant=matmul(matchange, coord_ant)

return
end subroutine pos_ant


subroutine rep_ant(var_ant,coord_ant,params)

! The aim of this function is to compute antenna manifold

implicit none

! Parameters structure
type(parameters)::params
! Antenna coordinates
real(kind=8),dimension(3,params%n_ant)::coord_ant
! Antenna manifold
complex(kind=8),dimension(params%n_ant,params%nazim)::var_ant
! Intermediate variables
complex(kind=8),dimension(params%n_ant)::amp_ant
real(kind=8),dimension(params%nazim)::azimr,tkx,tky
real(kind=8),dimension(2,params%nazim)::vectk
real(kind=8),dimension(3)::pos_centre
real(kind=8),dimension(3,params%n_ant)::pos_ant_relat
real(kind=8),dimension(2,params%n_ant)::pos_antxy
complex(kind=8),dimension(params%n_ant,params%n_ant)::d_amp_ant
real(kind=8)::Fc,B,Q,nu,k
integer::i_ant
real(kind=8)::pi
complex ::j,G


! Initializations
pi = acos(0.0)*2.
j=cmplx(0,1)

! Amplitude complexe propre a chaque antenne calcule ici avec freq_ant
! Resonateur suivant LC amorti
! Boucle sur le nombre d'antennes
do i_ant=1,params%n_ant
     Fc=params%freq_ant(i_ant,1)
     B =params%freq_ant(i_ant,2)
     Q=Fc/B
     nu=params%FE/Fc
     G=1./(1.+j*(Q/nu)*((nu**2.)-1.))
     amp_ant(i_ant)=G
enddo
 
! Module du vecteur d'onde
k=(2.*pi)/params%Lambda

! Tableau des azimuts en radian
! azimr = incidence en radian en sens trigo par rapport a la normale 
! au reseau qui est dirigee vers l axe des y
azimr=pi*params%azimd/180.

! Coordonnees i et j du vecteur k qui est oppose a la direction de pointe azimr
tkx=-cos(azimr)      
tky=-sin(azimr)

! Tableau du vecteur k (appartenant au plan horizontal) par azimut
vectk(1,:)=k*tkx 
vectk(2,:)=k*tky 

! pos_centre centre de gravite de l'antenne 
pos_centre=sum(coord_ant,2)/params%n_ant

! Position relative des antennes par rapport au centre de gravite
do i_ant=1,params%n_ant
     pos_ant_relat(:,i_ant)=coord_ant(:,i_ant)-pos_centre
enddo

! Coordonnees x et y de pos_ant_relat
pos_antxy=pos_ant_relat(1:2,:)

! Reponse des antennes en fonction de la direction = reponse(n_ant,params%nazim)
! Reponse = variete d antennes = signal toutes directions
! Reponse = alpha*exp(-jkx)

! Matrice diagonale amplitudes complexes des antennes selectionnees
! pour faciliter l ecriture du calcul de reponse
d_amp_ant=0.
do i_ant=1,params%n_ant
     d_amp_ant(i_ant,i_ant)=amp_ant(i_ant)
enddo

! Calcul de la variete d'antenne
! La normalisation sera faite dans le calcul de reponse_select
! selon le nombre d antenne selectionnees
! rep_ant=rep_ant/sqrt(params%NAT);
var_ant=matmul(d_amp_ant,exp(-j*(matmul(transpose((pos_antxy)), vectk))))

return
end subroutine rep_ant


subroutine hf_ant(sort_data,ant_hpower,ant_hcorr,params)
! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),dimension(params%n_chirp,params%n_ant,params%n_range)::sort_data
! Antenna power
real(kind=8),dimension(params%n_ant)::ant_hpower
! Antenna correction
real(kind=8),dimension(params%n_ant)::ant_hcorr
! Intermediate variables
real(kind=8),dimension(params%n_ant)::ant_hamp,ant_hnorm,int_mat
real(kind=8)::ant_href,ant_href1,corr_max,corr_min
integer,dimension(1)::ind_min
integer::n_min,n_raie
integer::i_ant


! Puissance haute frequence moyenne sur chaque antenne

! Il faut repousser au de la de la portee possible la bande d analyse
n_min=nint(0.75*params%n_range)
n_raie=params%n_range-n_min+1

! Boucle sur le nombre d'antennes
do i_ant=1,params%n_ant
     ! Moyenne des signaux au carre de n_min a params%n_range
     ant_hpower(i_ant)=sum(sum(abs(sort_data(:,i_ant,n_min:params%n_range))**2,1))/(params%n_chirp*n_raie)
enddo

! Amplitude des signaux
ant_hamp=sqrt(ant_hpower)

! Pour eviter les divergences
do i_ant=1,params%n_ant
     if(ant_hamp(i_ant)==0.)then
         ant_hamp(i_ant)=0.01
     endif
enddo

int_mat=ant_hamp

! S'il y une antenne de reference on la choisit sinon on prend la valeur
! mediane de l'amplitude
if(params%ant_etalon==-1) then
     call quick_select(int_mat,params%n_ant,ant_href)
     if(modulo(params%n_ant,2) .EQ. 0)then
         ind_min=minloc(abs(int_mat-ant_href))
         int_mat(ind_min(1))=999999999
         call quick_select(int_mat,params%n_ant,ant_href1)
         ant_href=(ant_href+ant_href1)/2.
     endif
else
     ant_href=int_mat(params%ant_etalon)
endif

! Amplitude normalisee
ant_hnorm=ant_hamp/ant_href

! On definit les bornes acceptables de la correction
corr_max=10.**(params%corr_max_dB/20.)
corr_min=1./corr_max

! On applique les bornes acceptables
do i_ant=1,params%n_ant
     if(ant_hnorm(i_ant)>corr_max)then
          ant_hnorm(i_ant)=corr_max     
     endif
     if(ant_hnorm(i_ant)<corr_min)then
         ant_hnorm(i_ant)=corr_min     
     endif
enddo

! Facteurs de correction applicables selon ce critere
ant_hcorr=1./ant_hnorm

return
end subroutine hf_ant


subroutine dc_ant(sort_data,ant_DCpower,params)
! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),dimension(params%n_chirp,params%n_ant,params%n_range)::sort_data
! Signal power of the first range cell for all antennas
real(kind=8),dimension(params%n_ant)::ant_DCpower
! Intermediate variables
integer::iant,ichirp

! Initialization
ant_DCpower=0.

! Loop over antennas
do iant=1,params%n_ant
     do ichirp=1,params%n_chirp
         ant_DCpower(iant)=abs(sort_data(ichirp,iant,1)**2.) 
     enddo
     ant_DCpower(iant)=ant_DCpower(iant)/params%n_chirp 
enddo

return
end subroutine dc_ant


subroutine select_ant(ant_hpower,ant_DCpower,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
real(kind=8),dimension(params%n_ant)::ant_hpower,ant_DCpower
! Intermediate variables
logical,dimension(params%n_ant)::ant_hselect,ant_DCselect,ant_select
integer,dimension(params%n_ant)::sort_index
real(kind=8)::ant_moyenne,ant_dyna,ant_mini,ant_maxi,ant_l
integer::iant,NAT,diffant,ierr,i
logical::diffok


! Antenna selection based on dynamic criterion between them
! ----------------------------------------------------------------------------
! If there is no reference antenna, the median value is taken
if(params%ant_etalon==-1)then
     sort_index=(/(i, i=1,params%n_ant, 1)/)
     call ssort(ant_hpower,sort_index,params%n_ant)
     
     if(mod(params%n_ant,2) .EQ. 0)then
         ant_moyenne=(ant_hpower(nint(aint(params%n_ant/2.)))+ant_hpower(nint(aint(params%n_ant/2.))+1))/2.
     else
         ant_moyenne=ant_hpower(nint(aint(params%n_ant/2.))+1)
     endif
else
     ant_moyenne=ant_hpower(params%ant_etalon)
endif

! Linear dynamic value tolerated between antennas
ant_dyna=10.**(params%ant_dyna_dB/10.)

! Maximum and minimum values
ant_mini=ant_moyenne/ant_dyna
ant_maxi=ant_moyenne*ant_dyna

! Loop over antennas
do iant=1,params%n_ant
     ant_l=ant_hpower(iant)

     ! Selection des antennes sur le critere de dynamique entre elles
     ant_hselect(iant)=(ant_l > ant_mini) .AND. (ant_l<ant_maxi)
enddo

! Antenna selection based on the DC level that should remain reasonable
! ---------------------------------------------------------------------------------------------
! If there is no reference antenna, the median value is taken
if(params%ant_etalon==-1)then
     sort_index=(/(i, i=1,params%n_ant, 1)/)
     call ssort(ant_DCpower,sort_index,params%n_ant)
     
     if(mod(params%n_ant,2) .EQ. 0)then
         ant_moyenne=(ant_DCpower(nint(aint(params%n_ant/2.)))+ant_DCpower(nint(aint(params%n_ant/2.))+1))/2.
     else
         ant_moyenne=ant_DCpower(nint(aint(params%n_ant/2.))+1)
     endif
else
     ant_moyenne=ant_DCpower(params%ant_etalon)
endif

! Maximum value
ant_maxi=ant_moyenne*100.

! Loop over antennas
do iant=1,params%n_ant
     ant_l=ant_DCpower(iant)
     ! The DC level that should remain reasonable
     ant_DCselect(iant)=(ant_l<ant_maxi)
enddo

! Tests synthesis
! --------------------
! Vector to cross test
ant_select=ant_hselect .AND. ant_DCselect

! List of selected antennas
NAT=0
do iant=1,params%n_ant
     if(ant_select(iant))then
         NAT=NAT+1
     endif
enddo

! Memory allocation
allocate(params%listant(NAT),stat=ierr)

NAT=0
do iant=1,params%n_ant
     if(ant_select(iant))then
         NAT=NAT+1
         params%listant(NAT)=iant
     endif
enddo

! Number of selected antennas
params%NAT=NAT

! At least 2 consecutives antennas are needed to be able to continue processing
do iant=1,params%NAT-1
     diffant=params%listant(iant+1)-params%listant(iant)

     if(diffant==1)then
         diffok=.FALSE.
         exit
     else
         diffok=.TRUE.
     endif
enddo

if(diffok)then
     print*, 'At least 2 consecutives antennas are needed to be able to continue processing'
endif

return
end subroutine select_ant


subroutine corr_ant(var_ant,ant_hcorr,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
! Antenna manifold
complex(kind=8),dimension(params%n_ant,params%nazim)::var_ant
! Antenna correction
real(kind=8),dimension(params%n_ant)::ant_hcorr
! Intermediate variables
integer::i_ant


! Auto correction en amplitude selon les amplitudes hf calculees sur antennes
if(params%arraycorr==1)then
     do i_ant=1,params%n_ant
         var_ant(i_ant,:)=ant_hcorr(i_ant)*var_ant(i_ant,:)
     enddo

! Les autres méthodes de correction ne sont pas (encore) implementees
else
     print*,'Correction de reponse du reseau pas implementee'
endif

return
end subroutine corr_ant



subroutine norm_ant(var_ant,params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
! Antenna manifold
complex(kind=8),dimension(params%NAT,params%nazim)::var_ant
! Intermediate variables
real(kind=8),dimension(1,params%nazim)::norm2_reponse,norm_reponse
real(kind=8),dimension(params%NAT,params%nazim)::mat_normalise
real(kind=8),dimension(params%NAT,1)::vect_int


! Calcul de la norme
norm2_reponse(1,:)=sum(conjg(var_ant)*var_ant, 1)
norm_reponse=sqrt(norm2_reponse)

! Normalisation pour tous les angles de la reponse du sous reseau
vect_int=1.
mat_normalise=matmul(vect_int,(1./norm_reponse))

! Reponse normalisee
var_ant=mat_normalise*var_ant

return
end subroutine norm_ant




end module antenna_module
