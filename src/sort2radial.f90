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


program sort2radial

! DESCRIPTION : the aim of this program is to derive radial current velocities 
! from input SORT or USORT Wera data using the MUSIC algorithm.
!
! INPUT : (U)SORT data files and parameters file
!
! OUTPUT : Radial velocities on polar grid


! Initializations
! -----------------
use init_params_module
use data_module
use antenna_module
use azim_module

implicit none

! Parameters structure
type(parameters)::params
! Sort data
complex(kind=8),allocatable,dimension(:,:,:)   ::sort_data
! Antenna manifold
complex(kind=8),allocatable,dimension(:,:)     ::var_ant
! Data after Doppler processing
complex(kind=8),allocatable,dimension(:,:,:,:)::data_chevfft
! Current map
real(kind=8),allocatable,dimension(:,:)            ::courant_dir
! Intermediate variables
integer :: ierr
complex(kind=8),allocatable,dimension(:,:,:,:)::data_int
complex(kind=8),allocatable,dimension(:,:)     ::var_int


print *,'Running sort2radial...'

! Initialization of the params structure
call init_params(params)

! Data loading
allocate(sort_data(params%n_chirp,params%n_ant,params%n_range),stat=ierr)
call load_data(sort_data,params)

! Antenna processing
allocate(var_int(params%n_ant,params%nazim),stat=ierr)
call proc_antenna(sort_data,var_int,params)

! Updates antenna manifold
allocate(var_ant(params%NAT,params%nazim),stat=ierr)
var_ant=var_int(params%listant,:)
deallocate(var_int)

! Antenna normalization
call norm_ant(var_ant,params)

! Doppler processing
allocate(data_chevfft(params%NPTS,params%n_ant,params%NSUB,params%n_range),stat=ierr)
allocate(data_int(params%NPTS,params%n_ant,params%NSUB,params%n_range),stat=ierr)
call doppler_proc(sort_data,data_chevfft,params)
deallocate(sort_data)

! Updates data
data_int=data_chevfft
deallocate(data_chevfft)
allocate(data_chevfft(params%NPTS,params%n_ant,params%n_OBS,params%n_range),stat=ierr)
data_chevfft=data_int(:,:,params%list_OBS_ok,:)
deallocate(data_int)

! Azimuth processing
allocate(courant_dir(params%prof,params%n_dir),stat=ierr)
call azim_proc(data_chevfft,var_ant,courant_dir,params)

! Map processing
call rad2ascii(courant_dir,params)

print *,'sort2radial done, result in :',trim(params%output_file)

! Exits program
call exit

end program sort2radial
