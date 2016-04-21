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

module data_module

! The aim of this module is to load the (U)SORT data files, to carry out the Doppler processing
! and to select independant samples.


! Initializations
! -----------------
use init_params_module

implicit none


contains


subroutine load_data(sort_data,params)

! The aim of this module is to load the (U)SORT data files

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),dimension(params%n_chirp,params%n_ant,params%n_range)::sort_data
! Intermediate variables
character(132)::int_char
integer::ind

! Name of the first data file of the list
int_char=params%list_files(1)

! Checks whether it is a SORT or USORT data file
do ind=len(int_char),1,-1
     if(int_char(ind:ind) .EQ. 'S')then
         exit
     endif
enddo

! Loads SORT or USORT data file
if(int_char(ind-1:ind-1) .EQ. 'U')then
     call load_usort(sort_data,params)
else
     call load_sort(sort_data,params)
endif

end subroutine load_data
   

subroutine load_sort(sort_data,params)

! The aim of this module is to load the SORT data files

implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),dimension(params%n_chirp,params%n_ant,params%n_range)::sort_data
! Intermediate variables to load data
real(kind=4),allocatable,dimension(:)::int_data1
real(kind=4),allocatable,dimension(:)::int_data2
real(kind=4),allocatable,dimension(:,:,:,:)::mat_data1
real(kind=4),allocatable,dimension(:,:,:,:)::mat_data2
real(kind=4),allocatable,dimension(:,:,:,:)::mat_data
integer::size_data,size_data_file,ierr
integer::a,b,c,ind,ind_chirp
complex::j

! Initialization
j=cmplx(0,1)

! Number of elements to be loaded from the SORT file
size_data=params%n_chirp*params%n_ant*params%n_range*2
! Total number of elements of the SORT file including the header
size_data_file=params%n_chirp_file*params%n_ant*params%n_range*2+128

! Memory allocation
allocate(mat_data(2,params%n_chirp,params%n_ant,params%n_range),stat=ierr)
allocate(int_data1(size_data_file),stat=ierr)
allocate(int_data2(size_data_file-128),stat=ierr)
allocate(mat_data1(2,params%n_chirp_file,params%n_ant,params%n_range),stat=ierr)


! Loop over the number of input data files
do ind=1,params%n_files
     ! Loads the entire SORT file
     open(ind_file_file,file=trim(params%list_files(ind)),access='direct',recl=4*size_data_file)
     read(ind_file_file,rec=1)int_data1
     close(ind_file_file)
     
     ! Keeps only valuable data excluding the header
     int_data2=int_data1(129:size_data)

     ! Reshapes data to the correct format
     mat_data1=reshape(int_data2,(/ 2, params%n_chirp_file,params%n_ant,params%n_range/))
     
     ! Fills in the entire matrix
     do ind_chirp=1,params%n_chirp_file
         mat_data(:,params%n_chirp_file*(ind-1)+ind_chirp,:,:)=mat_data1(:,ind_chirp,:,:)
     enddo
     
enddo

! Memory deallocation
deallocate(int_data1)
deallocate(int_data2)
deallocate(mat_data1)

! Tests if reversing the antennas is needed (Brezellec case after 2005)
if(params%rev_ant)then
     ! Memory allocation
     allocate(mat_data2(2,params%n_chirp,params%n_ant,params%n_range),stat=ierr)
     
     ! Temporary saving
     mat_data2=mat_data
 
     ! Loop over antennas
     do a = 1,params%n_ant
         mat_data(:, :, a, :)=mat_data2(:, :, 17-a, :)
     enddo
 
     ! Memory deallocation
     deallocate(mat_data2)
endif

! Computes the complex matrix
do a=1,params%n_range
     do b =1,params%n_ant
         do c=1,params%n_chirp
             sort_data(c,b,a)=dble(mat_data(1,c,b,a))+j*dble(mat_data(2,c,b,a))
         enddo
     enddo
enddo

! Memory deallocation
deallocate(mat_data)

return
end subroutine load_sort



subroutine load_usort(sort_data,params)

! The aim of this module is to load the USORT data files

implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),dimension(params%n_chirp,params%n_ant,params%n_range)::sort_data
! Intermediate variables to load data
real(kind=4),allocatable,dimension(:)::int_data1
real(kind=4),allocatable,dimension(:)::int_data2
real(kind=4),allocatable,dimension(:,:,:,:)::mat_data1
real(kind=4),allocatable,dimension(:,:,:,:)::mat_data2
real(kind=4),allocatable,dimension(:,:,:,:)::mat_data
real(kind=4),dimension(3,params%n_ant)::calib_data_trans
real(kind=4),dimension(params%n_ant,3)::calib_data
integer::size_data,size_data_file,ierr
integer::a,b,c,ind,ind_chirp,ind_range,ind_ant
real(kind=8)::pi
complex::j

! Initializations
j=cmplx(0,1)
pi = acos(0.0)*2.

! Number of elements to be loaded from the USORT file
size_data=params%n_chirp*params%n_ant*params%n_range*2
size_data_file=params%n_chirp_file*params%n_ant*params%n_range*2+128

! Memory allocation
allocate(mat_data(2,params%n_range,params%n_ant,params%n_chirp),stat=ierr)
allocate(int_data1(size_data_file),stat=ierr)
allocate(int_data2(size_data_file-128),stat=ierr)
allocate(mat_data1(2,params%n_range,params%n_ant,params%n_chirp_file),stat=ierr)

! Loop over the number of input data files
do ind=1,params%n_files
    
     ! Loads the entire USORT file
     open(ind_file_file,file=trim(params%list_files(ind)),access='direct',recl=4*size_data_file)
     read(ind_file_file,rec=1)int_data1
     close(ind_file_file)

     ! Keeps only valuable data excluding the header
     int_data2=int_data1(129:size_data_file)

     ! Reshapes data to the correct format
     mat_data1=reshape(int_data2,(/2,params%n_range,params%n_ant,params%n_chirp_file/))
    
     ! Fills in the entire matrix
     do ind_chirp=1,params%n_chirp_file
         mat_data(:,:,:,params%n_chirp_file*(ind-1)+ind_chirp)=mat_data1(:,:,:,ind_chirp)
     enddo
     
enddo

! Memory deallocation
deallocate(int_data1)
deallocate(int_data2)
deallocate(mat_data1)


! Tests if reversing the antennas is needed (Brezellec case after 2005)
if(params%rev_ant)then
     ! Memory allocation
     allocate(mat_data2(2,params%n_range,params%n_ant,params%n_chirp),stat=ierr)
     
     ! Temporary saving
     mat_data2=mat_data
 
     ! Loop over antennas
     do a = 1,params%n_ant
         mat_data(:, :, a, :)=mat_data2(:, :, 17-a, :)
     enddo
 
     ! Memory deallocation
     deallocate(mat_data2)
endif

! Memory allocation
allocate(mat_data2(2,params%n_chirp,params%n_ant,params%n_range),stat=ierr)

! Puts in the correct order
do a=1,params%n_range
     do b =1,params%n_ant
         do c=1,params%n_chirp
             mat_data2(:,c,b,a)=mat_data(:,a,b,c)
         enddo
     enddo
enddo

! Computes the complex matrix
do a=1,params%n_range
     do b =1,params%n_ant
         do c=1,params%n_chirp
             sort_data(c,b,a)=dble(mat_data2(1,c,b,a))+j*dble(mat_data2(2,c,b,a))
         enddo
     enddo
enddo

deallocate(mat_data2)
deallocate(mat_data)

! Loads the calibration file
open(ind_file_file,file=trim(params%calib_file))
read(ind_file_file,*)calib_data_trans
close(ind_file_file)

! Puts in the correct order
calib_data=transpose(calib_data_trans)

! Loops over the entire domain to apply the calibration file
do ind_range=1,params%n_range
     do ind_ant=1,params%n_ant
         do ind_chirp=1,params%n_chirp
             sort_data(ind_chirp, ind_ant, ind_range)=sort_data(ind_chirp, ind_ant, ind_range)*exp(j*calib_data(ind_ant, 3)*pi/180.)
         enddo
     enddo
enddo

return
end subroutine load_usort


subroutine doppler_proc(sort_data,data_chevfft,params)

! The aim of this function is to perform the doppler processing to the input data

implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),dimension(params%n_chirp,params%n_ant,params%n_range)::sort_data
! Data after Doppler processing
complex(kind=8),dimension(params%NPTS,params%n_ant,params%NSUB,params%n_range)::data_chevfft
! Intermediate variables
integer::idebut0,i,i_range,i_sub,i_ant,ierr,i_echant
integer::ind_test,ind_int,ind_all,i_zero
logical::val_test
real(kind=8)::pi
real(kind=8),dimension(params%NPTS)::Tour,wbtime
integer::lensav,lenwrk
real(kind=8),allocatable,dimension(:)::work
real(kind=8),allocatable,dimension(:)::wsave
integer,allocatable,dimension(:)::vect_all,vect_int


! Initializations
! -----------------
pi = acos(0.0)*2.
! Parameters for carrying out FFT processing
lensav=2*params%NPTS+int(log(real(params%NPTS,kind=8)))+4
lenwrk=2*params%NPTS
allocate(wsave(lensav),stat=ierr)
allocate(work(lenwrk),stat=ierr)
call zfft1i(params%NPTS,wsave,lensav,ierr)

! Fills in the matrix regarding Kchevauch and vaca_div parameters
do i_sub=1,params%NSUB
     idebut0=floor(((i_sub-1)*params%NPTS)/real(params%Kchevauch,kind=8))
     data_chevfft(:,:,i_sub,:)=sort_data(idebut0+(/(i, i=1,params%NPTS, 1)/),:,:)
enddo
 
! Computes weigthing parameters to avoid secondary lobes in the Doppler spectrum from 
! 0 to 2*pi-epsilon
Tour=2*pi*(/(i, i=0,params%NPTS-1, 1)/)/params%NPTS

! Blackman window
wbtime=0.35875-(0.48829*cos(Tour))+ (0.14128*cos(2.*Tour))-(0.01168*cos(3.*Tour))

! Loops over the entire domain to apply the window function
do i_range=1,params%n_range
     do i_sub=1,params%NSUB
         do i_ant=1,params%n_ant          
             data_chevfft(:,i_ant,i_sub,i_range)=wbtime*data_chevfft(:,i_ant,i_sub,i_range)
         enddo
     enddo
enddo

! Loops over the entire domain to perform FFT processing
do i_range=1,params%n_range
     do i_sub=1,params%NSUB
         do i_ant=1,params%n_ant
             call zfft1f(params%NPTS,1,data_chevfft(:,i_ant,i_sub,i_range),params%NPTS,wsave,lensav,work,lenwrk,ierr)         
         enddo
     enddo
enddo

! Scaling due to chosen FFT program
data_chevfft=params%NPTS*data_chevfft

! Test is there is a selection of observed samples or not 
if(params%do_selection_OBS)then
     call select_samp(data_chevfft,params)
else
     allocate(params%list_OBS_ok(params%NSUB),stat=ierr)
     params%list_OBS_ok=(/(i, i=1,params%NSUB, 1)/)
endif

! Total number of selected samples
params%n_OBS=size(params%list_OBS_ok)

! Computes the equivalent number of independant observations
! ----------------------------------------------------------------------------------
! Memory allocation
allocate(vect_all(params%NPTS*params%vaca_div),stat=ierr)
allocate(vect_int(params%NPTS),stat=ierr)

! Initialization
ind_test=0
vect_all=0

! Loop over selected sample
do i_echant=1,params%n_OBS
     ! General index of the considered sample
     i_zero=floor(((params%list_OBS_ok(i_echant)-1)*params%NPTS)/real(params%Kchevauch,kind=8))
     
     ! Intermediate vector
     vect_int=i_zero+(/(i, i=1,params%NPTS, 1)/)
     
     ! Performs "unique" equivalent Matlab function (without sorting)
     do ind_int=1,params%NPTS
         val_test=.TRUE.
         do ind_all=1,params%NPTS*params%vaca_div
             if(vect_int(ind_int) .EQ. vect_all(ind_all))then
                 val_TEST=.FALSE.
                 exit
             endif
         enddo
         if(val_test)then
             ind_test=ind_test+1
             vect_all(ind_test)=vect_int(ind_int)
         endif
     enddo
enddo

!  "Truly" independant samples
! /2 due to Blackman windowing (to be confirmed)
params%q_OBS=nint(ind_test/(params%NPTS/2.))

! Minimum between the two parameters
params%q_OBS=min(params%q_OBS,params%n_OBS)

! Empirical choice for dB_snr for selecting Doppler lines depending on the equivalent
! number of independant samples
! Test if automatic scaling or not
if(params%dB_snr0==-1)then
     if(params%q_OBS<6)then
         params%dB_snr=3.
     else if (params%q_OBS<10)then
         params%dB_snr=2.
     else if (params%q_OBS<15)then
         params%dB_snr=1.5
     else
         params%dB_snr=1. 
     endif         
else
     params%dB_snr=params%dB_snr0
endif

return
end subroutine doppler_proc



subroutine select_samp(data_chevfft,params)

! The aim of this function is to determine the list of independant samples

implicit none

! Parameters structure
type(parameters)::params
! Data after Doppler processing
complex(kind=8),dimension(params%NPTS,params%n_ant,params%NSUB,params%n_range)::data_chevfft
! Intermediate variables
integer::n_listant,n_raies,n_porte
integer::i_echant,i_range,i_ant,ierr,i
integer::rangdernier,nmini
complex(kind=8),allocatable,dimension(:)::robs
real(kind=8)::robs_niveau
real(kind=8),dimension(params%NSUB)::niveau_bruit,niveau_bruit_dB,sort_niveau
integer,dimension(params%NSUB)::sort_index


! Initializations
n_listant =size(params%listant)
n_raies =size(params%liste_raies_bruit)
n_porte =size(params%liste_range)

! Memory allocation
allocate(robs(n_raies),stat=ierr)

! Loops over the entire domain
do i_echant=1,params%NSUB
     do i_range=1,n_porte
         do i_ant=1,n_listant
             ! Takes the associated data
             robs=data_chevfft(params%liste_raies_bruit,params%listant(i_ant),i_echant,params%liste_range(i_range))

             ! Sum of energies
             robs_niveau=sum(robs*conjg(robs))
                  
             ! Updates the vector of level of noise
             niveau_bruit(i_echant)=niveau_bruit(i_echant)+robs_niveau
         enddo
     enddo
enddo

! Memory deallocation
deallocate(robs)

! Mean level of noise and conversion to dB
niveau_bruit =niveau_bruit/(n_listant*n_raies*n_porte)
niveau_bruit_dB=10.*log10(niveau_bruit)

! Data sorting
sort_niveau=niveau_bruit_dB
sort_index=(/(i, i=1,params%NSUB, 1)/)
call ssort(sort_niveau,sort_index,params%NSUB)

! Normalization
sort_niveau=sort_niveau-sort_niveau(1)

! Determines the index of the last sample that fills in the condition 
do  i_echant=params%NSUB,1,-1
     if(sort_niveau(i_echant)<params%dB_select_OBS)then
         exit
     endif
enddo
rangdernier=i_echant

! Checks if the number of samples if good enough compared to the number of selected antennas.
! If not, the minimal number of selected samples is forced.
nmini=nint(params%NSUB/2.)
nmini=min(nmini,(n_listant-1))
if((rangdernier)<nmini)then
     rangdernier=nmini
endif

! Memory allocation
allocate(params%list_OBS_ok(rangdernier),stat=ierr)

! List of selected samples
params%list_OBS_ok=sort_index(1:rangdernier)

return
end subroutine select_samp

end module data_module
