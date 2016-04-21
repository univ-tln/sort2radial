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


module init_params_module

! The aim of this module is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Definition of the params structure
type parameters
   ! Input parameters
   ! -----------------------
   ! Input filenames
   character(132),allocatable,dimension(:)::list_files ! List of files to be processed
   character(132)                          ::calib_file            ! Filename of the calibration file only needed if USORT data files 
   character(132)                          ::pos_ant_file       ! Coordinates file of the Rx array with the origin centered on the 0 antenna [mod_ant]
   character(132)                          ::mask_file           ! Velocity mask file
   character(132)                          ::grid_file             ! File containing the cartesian grid coordinates
   character(132)                          ::alcov2_file         ! alcov2 file
   ! Output filename
   character(132)                          ::output_file          ! Output filename
   ! Mandatory parameters
   character(132)                          ::site_name          ! Site name
   logical                                       ::rev_ant               ! Test if reversing the antennas is needed (Brezellec case)
   logical                                       ::verbose             ! Verbose (boolean)
   logical                                       ::figures               ! Intermediate visualization figures (boolean)
   real(kind=8)                              ::distmin              ! Minimal range of interest (km)
   real(kind=8)                              ::distmax             ! Maximal range of interest (km)
   real(kind=8)                              ::vmax                 ! Maximal radial velocity expected for the current (m/s) 
   real(kind=8)                              ::dteta                 ! Azimutal increment (degree)

   ! Advanced parameters
   ! ----------------------------
   ! Antenna processing
   real(kind=8),allocatable,dimension(:,:)    ::freq_ant           ! Resonance frequency and bandwidth per antenna (MHz) (used to compute the antenna manifold)
   integer                                                     ::ant_etalon        ! Number of the reference antenna (0 ... n_ant-1). If (-1), the median value of power is taken
   real(kind=8)                                            ::z_ant                ! Altitude of the RX array (meter)
   real(kind=8)                                            ::lcable               ! Length of the cables TX+system+RX (meter) 
   real(kind=8)                                            ::ncable              ! Propagation property of the cable=1/permittivity
   real(kind=8)                                            ::fclow                ! High-pass cut-off frequency of the first order in baseband (Hz)
   logical                                                      ::ya_sym            ! Symetry of the reception array (boolean)
   integer                                                     ::arraycorr         ! Correction of the Rx array
   real(kind=8)                                            ::corr_max_dB    ! Maximal amplitude correction for the antennas (dB)
   logical                                                      ::select_ant        ! Automatic selection of the antennas (boolean)
   real(kind=8)                                            ::ant_dyna_dB     ! Half-dynamic value allowed between antennas (dB)
   ! Doppler processing
   real(kind=8)                                            ::delta_rayon     ! Radial depth of a joint crown (km)
   logical                                                     ::modeg            ! Sliding mode for crowns (boolean)
   integer                                                    ::vaca_div          ! Spliting ratio of the data file for subseries
   integer                                                    ::Kchevauch       !  Kchevauch=1/(shift rate of the successive samples)
   logical                                                     ::do_selection_OBS ! Selection of the observation samples based on a noise criterion out of the Bragg area (boolean)
   real(kind=8)                                           ::dB_select_OBS     ! Dynamic of selection of the samples (subseries of overlapping chirps) in order to get good observations (dB)
   ! Azimuth processing
   integer                                                    ::methselect       ! Choice of the selection method for determining Doppler lines of interest and the level of noise
   real(kind=8)                                           ::valseuil             ! Level of asymptotic noise, useful for  methselect=1, 2 or 3
   real(kind=8)                                           ::dB_dynaraie      !  Maximal dynamic value between the strongest line and the smallest signal accepted for each doppler domain (dB)
   real(kind=8)                                           ::dB_snr0            !  Selection rate value of lines above noise. If -1 : automatic choice (dB)
   real(kind=8)                                           ::dB_snr              !  Chosen value (dB)
   integer                                                    ::nbbosse          ! Number of bosses for Doppler line selection
   real(kind=8)                                           ::dB_new_pic      !  Value to declare a new peak (dB)
   real(kind=8)                                           ::dB_inter_pic      !  Maximal dynamic value between local peaks on the same domain (dB)
   integer                                                    ::toler_con          ! Non connexity tolerancy (>=1) to know the number of point tolerated under the level of noise to determine the number of lines (number of points) 
   real(kind=8)                                           ::dB_anis             !  Maximum anisotropy value between the two domains (dB)
   integer                                                    ::nb_source_max  ! Maximum number of sources considered for each Doppler line of interest
   integer                                                    ::n_liss0             ! Length of the filter (in points) to smooth the MUSIC factor
   logical                                                     ::egalispectre    ! Equalization right-left of the Doppler lines (boolean)
   real(kind=8)                                           ::dB_dynadir       !  Dynamic value for acceptation of the solution per crown on the music map compared to the mean of the crown (dB)
   real(kind=8)                                           ::dB_mini            !  Minimum value of energy for source acceptation (dB)
   ! Miscellaneous parameters
   integer                                                    ::headersize      ! Size of the header of the data file (bytes)
   integer                                                    ::floatsize          ! Float size (bytes)
   integer                                                    ::slave_bit          ! Slave bit (to determine upward/downard chirp)
   real(kind=8)                                           ::c_light              !  Air celerity (m/s)
   real(kind=8)                                           ::dazim               !  Step size (degree)
   real(kind=8)                                           ::azimmin            ! Research aperture (degree)
   real(kind=8)                                           ::azimmax           ! Research aperture (degree)
   real(kind=8)                                           ::kzero               ! Coefficient of the width of the zero Doppler line to be rejected
   real(kind=8)                                           ::kzero_wide      ! Coefficient for larger width of the zero Doppler line to be rejected 
   real(kind=8)                                           ::n_range_min     ! Coefficient to determine the minimal distance to study the level of noise allowing to select the observed samples
   integer                                                    ::multiv2f           ! Maximum number of sources in a grid cell (for memory allocation) 
   real(kind=8)                                           ::beam_width     ! Multiplicative factor to determine empirically the width of the beam (depends on site)
   real(kind=8),dimension(3)                      ::proba_ref         ! Parameters for fct_select_raies
   real(kind=8),dimension(3)                      ::kmargin            ! Parameters for fct_select_raies
   real(kind=8)                                           ::pfalarm             ! Parameters for fct_select_raies
   integer                                                    ::accent              ! According to the form of the distribution, if white: 0 is enough. 
   real(kind=8)                                           ::ksup                 ! Criteria border for the equalization of the spectrum
   real(kind=8)                                           ::kinf                   ! Criteria border for the equalization of the spectrum
   real(kind=8)                                           ::satur_gain         ! Saturation of the gain to compensate the effect of the high-pass filter. Used only if methselect=3
   logical                                                     ::test_mask        ! Test if there is a velocity mask or not

   ! Header parameters
   ! --------------------------
   integer                                                    ::n_chirp_file      ! Number of chirps
   character(132)                                        ::date_str          ! Date 
   character(132)                                        ::heure_str        ! Time 
   character(132)                                        ::nom_camp     ! Name of the campaign
   character(132)                                        ::fm_type          ! Code of the type of data in the header
   real(kind=8)                                           ::FE                   ! Frequency (MHz)
   real(kind=8)                                           ::F                     ! Frequency (Hz)
   character(132)                                        ::annee             ! Year
   real(kind=8)                                           ::resolution        ! Range resolution (km)
   real(kind=8)                                           ::cap_bous        ! Angular position of the position for the Rx linear array from antenna 0 to the last one (degree)
   real(kind=8)                                           ::rate                ! Chirp duration (s)
   integer                                                   ::n_range          ! Nombre de cellules distance
   character(132)                                       ::longitude         ! Longitude position of the radar
   character(132)                                       ::latitude            ! Latitude position of the radar
   integer                                                    ::MT                  ! Nb of samples per chirp
   real(kind=8)                                           ::pwr                 !  Transmitted power
   integer                                                    ::n_ant              ! Number of antennas
   integer                                                    ::md                 ! The "mode bits word" MD is hexadecimally coded
   real(kind=8)                                           ::offsett            !  Offset
   real(kind=8)                                           ::rxoffset          !  RX offset
   character(132)                                       ::comment        ! Comment

   ! Computed parameters
   ! ------------------------------
   integer                                                    ::n_chirp           ! Total number of chirps
   real(kind=8)                                           ::tetamax          !  Definition of the angular sector for display (degree, "trigonometry reference") 
   real(kind=8)                                           ::tetamin           !  Definition of the angular sector for display (degree, "trigonometry reference") 
   real(kind=8)                                           ::Lambda          !  Wavelength (meter)
   real(kind=8)                                           ::slave_state     !  Determination of the direction of the sweep of the chirp
   real(kind=8)                                           ::chirp_direction  !  Determination of the direction of the sweep of the chirp
   integer                                                    ::reg_radial       ! Computation of the number gathered crowns (at least 1)
   integer,allocatable,dimension(:)              ::porte              ! Number of gathered crowns to consider
   integer                                                    ::prof                ! Number of studied crowns
   real(kind=8),allocatable,dimension(:)     ::hdist                 ! Range of the radial distances "as the crow flies" , depends on crowns
   real(kind=8)                                           ::decalage         ! If self.reg_radial is odd, the crown distance is the central gate one
   integer                                                    ::NPTS               ! Number of chirps in a subserie (vacation/vacadiv)
   integer                                                    ::NSUB              ! Number of subseries of measure resulting for a file
   integer                                                    ::raie_pos_max  ! Maximum number of the positive line
   integer                                                    ::raie_neg_min   ! Minimum number of the negative line
   real(kind=8)                                           ::DeltaF              ! Frequential resolution DeltaF = 1/Integration time 
   real(kind=8)                                           ::FB                    ! Bragg frequency in Hz
   integer                                                    ::DIbragg           ! Bragg line for z_ant=0
   integer,dimension(2)                              ::Ind_bragg        ! Number of the Bragg line (two number ordered from min to max)
   real(kind=8),allocatable,dimension(:)     ::lcoefnoise        ! Coeff in dB applicable according the distance to the level of noise (only if methselect=3)
   real(kind=8),allocatable,dimension(:,:)   ::vradz2D         ! Table to take into account the correction of altitude for velocity computation
   integer,allocatable,dimension(:)              ::liste_raies_zero,liste_raies_zero_wide              ! Width of the forbidden area
   integer,allocatable,dimension(:)              ::liste_raies_bruit     ! The zero area of the largest width is not taken into account
   integer,allocatable,dimension(:)              ::liste_range         !  List range for the selection of the sample, base on a level of energy to far distances
   integer                                                    ::n_dir                  ! Number of angular sectors
   real(kind=8),allocatable,dimension(:)     ::tetad                  !  Table of angular sectors
   integer                                                    ::nazim                  ! Number of thin azimuths
   real(kind=8),allocatable,dimension(:)     ::azimd                  !  Table of azimuth in degree
   real(kind=8),allocatable,dimension(:)     ::hors_jeu_sym       ! Calculation of the bench area
   integer,allocatable,dimension(:)              ::listant                  ! List of antennaes
   integer                                                    ::NAT                    ! Number of antennaes
   integer,allocatable,dimension(:)              ::list_OBS_ok         ! List of selected samples
   integer                                                    ::n_OBS                 ! Number of selected samples
   integer                                                    ::q_OBS                 ! Number of truly independent samples
   integer,allocatable,dimension(:)              ::liste_raies_vague
   integer,allocatable,dimension(:)              ::support_vague
   integer                                                    ::demi_zero
   real(kind=8),allocatable,dimension(:,:)    ::TVMAX
   integer,allocatable,dimension(:)              ::idomaine
   integer                                                    ::total_raies
   integer                                                    ::total_sources
   integer                                                    ::not_nan
   character(132)                                        ::map_func
   integer                                                    ::n_files
   logical                                                     ::read_alcov2
   real(kind=8),dimension(7,3,2)                ::alcov2
end type parameters


! Global variables
! ----------------
! Indices for files opening
integer,parameter::ind_list_files=100
integer,parameter::ind_file_file=101
integer,parameter::ind_file_calib=102
integer,parameter::ind_file_pos_ant=103
integer,parameter::ind_file_mask=104
integer,parameter::ind_file_alcov2=105


contains


subroutine init_params(params)

! The aim of this function is to build the params structure that contains
! the main parameters of the simulation

! Initializations
! ---------------
implicit none

! Parameters structure
type(parameters)::params
! Input filenames
character(132)    ::calib_file,pos_ant_file,mask_file,grid_file,alcov2_file          
! Mandatory parameters
character(132)    ::site_name,map_func
logical                 ::rev_ant,verbose
real(kind=8)       ::distmin,distmax,vmax,dteta
integer                ::ind,ierr


! Reading of the input parameters file
read(*,*)
read(*,*)calib_file
read(*,*)
read(*,*)pos_ant_file
read(*,*)
read(*,*)mask_file
read(*,*)
read(*,*)alcov2_file 
read(*,*)
read(*,*)site_name
read(*,*)
read(*,*)rev_ant
read(*,*)
read(*,*)verbose
read(*,*)
read(*,*)distmin
read(*,*)
read(*,*)distmax
read(*,*)
read(*,*)vmax
read(*,*)
read(*,*)dteta
            

! Fills in the params structure
! Input filenames
params%calib_file=calib_file
params%pos_ant_file=pos_ant_file
params%mask_file=mask_file
params%grid_file=grid_file
params%alcov2_file=alcov2_file
! Mandatory parameters
params%site_name=site_name
params%rev_ant=rev_ant
params%verbose=verbose
params%distmin=distmin
params%distmax=distmax
params%vmax=vmax
params%dteta=dteta
params%map_func=map_func

! Determines the number of files to be processed by counting the number of lines
params%n_files=iargc()-1
allocate(params%list_files(params%n_files),stat=ierr)

call getarg(1, params%output_file)

do ind=2,params%n_files+1
     call getarg(ind, params%list_files(ind-1))
enddo

! Advanced parameters needed for the simulation
call advanced_parameters(params) ! Function located in init_params_module

! Reads parameters from header file
call header_parameters(params) ! Function located in init_params_module

! Computes more complex parameters needed for the simulation
call computed_parameters(params) ! Function located in init_params_module

return
end subroutine init_params


subroutine advanced_parameters(params)

! The aim of the this function is to provide advanced parameters for the simulation

! Initializations
implicit none

! Parameters structure
type(parameters)                    ::params


! Antenna processing
! --------------------------
! Number of the reference antenna (0 ... n_ant-1). If (-1), the median value of power is taken [mod_ant]
params%ant_etalon =-1

! Altitude of the RX array (meter). Used to derive radial velocity values [mod_par]
params%z_ant =0.

! Length of the cables TX+system+RX (meter) [mod_par]
params%lcable=(250.+10.+200.)

! Propagation property of the cable=1/permittivity [mod_par]
params%ncable=3./2.

! High-pass cut-off frequency of the first order in baseband (Hz) [mod_par]
params%fclow=67.72

! Symetry of the reception array (boolean) [mod_par, fct_sourc_cour]
params%ya_sym=.TRUE.

! Correction of the Rx array
!   0 = perfect array
!   1 = auto-correction
!   2 = use of a calibration file (not implemented)
!   3 = use of a file of measure of lobe (not implemented)
!   4 = use of an attenuator file (not implemented)
params%arraycorr=1

! Maximal amplitude correction for the antennas (dB) [mod_ant]
params%corr_max_dB=10.

! Automatic selection of the antennas (boolean) [mod_ant]
params%select_ant=.FALSE. 

! Half-dynamic value allowed between antennas (dB). Used only if params%select_ant=.TRUE. [mod_ant]
params%ant_dyna_dB=3.

! Doppler processing
! -------------------------
! Radial depth of a joint crown (km) [mod_par]
params%delta_rayon=1.5

! Sliding mode for crowns (boolean) [mod_par]
params%modeg=.TRUE.

! Spliting ratio of the data file for subseries
! Sample of length = length(vacation)/vaca_div [mod_par, def_vp_noise]
params%vaca_div=2

! Kchevauch=1/(shift rate of the successive samples) [mod_par, mod_data, fct_vp_noise]
params%Kchevauch=2

! Selection of the observation samples based on a noise criterion out of the Bragg area (boolean) [mod_data]
params%do_selection_OBS=.TRUE.

! Dynamic of selection of the samples (subseries of overlapping chirps) in order to get good observations (dB)
! Used only if params%do_selection_OBS=.TRUE. [mod_data]
params%dB_select_OBS=6.


! Azimuth processing
! -------------------------
! Choice of the selection method for determining Doppler lines of interest and the level of noise [fct_select_raies]
!   0 automatic
!   1 automatic selection of the lines and definition of a constant level of noise 
!   2  constant level of noise
!   3 level of noise function of radial ranges 
params%methselect=0

!         if((params%methselect==1) | (params%methselect==2) | (params%methselect==3)):
! Level of asymptotic noise, useful for  params%methselect=1, 2 or 3 (dB) [fct_select_raies]
params%valseuil=21.5

! Maximal dynamic value between the strongest line and the smallest signal accepted for each doppler domain (dB) [fct_select_raies]
params%dB_dynaraie=20.

! Selection rate value of lines above noise. If -1 : automatic choice (dB) [mod_data]
params%dB_snr0=-1

! Number of bosses to select Doppler lines
params%nbbosse=1

! Value to declare a new peak (dB) [fct_select_raies]
params%dB_new_pic=1.

! Maximal dynamic value between local peaks on the same domain (dB) [fct_select_raies]
params%dB_inter_pic=26.

! Non connexity tolerancy (>=1) to know the number of point tolerated under the level of noise to determine the number of lines (number of points) [fct_select_raies]
params%toler_con=1

! Maximum anisotropy value between the two domains [fct_select_raies, fct_multi_courant]
params%dB_anis=20.

! Maximum number of sources considered for each Doppler line of interest [fct_sourc_cour]. 
! For the moment, it is impossible to take more because of alcov2, which is calculated for a maximum number of 4 sources
params%nb_source_max=1

! Length of the filter (in points) to smooth the MUSIC factor [fct_sourc_cour]
params%n_liss0=7

! Equalization right-left of the Doppler lines (boolean) [def_multi_courant]
params%egalispectre=.False.

! Dynamic value for acceptation of the solution per crown on the music map compared to the mean of the crown 
! Generally between 10 and 15, you can take more only if there is a masking phenomenon as in the Frioul case.
! Detection of ground echoes according to the dynamic modification in the crown (dB) [fct_un_courant_music]
params%dB_dynadir=121.

! Minimum value of energy for source acceptation (dB) [fct_un_courant_music]
params%dB_mini=15.

! Read alcov3
params%read_alcov2=.TRUE.


! Miscellaneous parameters
! ----------------------------------
! Size of the header of the data file (bytes) [mod_par]
params%headersize=512

! Float size (bytes) [mod_par]
params%floatsize=4

! Slave bit (to determine upward/downard chirp) [mod_par]
params%slave_bit=128

! Air celerity (m/s) [several methods]
params%c_light=3.e8

! Angular parameters for computing the MUSIC factor
! Total azimutal coverage on 2*Pi for all sites for MUSIC
! Step size (degree) [mod_ant, mod_par]
params%dazim=0.5
! Research aperture (degree)
params%azimmin=0.                        ![var_ant, comp_params]
params%azimmax=360.-params%dazim ! [var_ant]

! Coefficient of the width of the zero Doppler line to be rejected [mod_par, fct_select_raies]
params%kzero=1/128.

! Coefficient for larger width of the zero Doppler line to be rejected [mod_par, fct_select_raies]
params%kzero_wide=4.*params%kzero

! Coefficient to determine the minimal distance to study the level of noise allowing 
! to select the observed samples (distmin=params%n_range_min*params%range) [mod_par]
! Peyras : 1/4, Garchine : 1/2
params%n_range_min=0.25

! Maximum number of sources in a grid cell (for memory allocation) [fct_azim_processing, fct_sourc_cour, fct_multi_courant]
params%multiv2f=600

! Multiplicative factor to determine empirically the width of the beam (depends on site) [sourc_cour]
params%beam_width=2.

! Parameters for fct_select_raies (analogy false alarm probability /detection probability) [fct_select_raies]
params%proba_ref(1)=0.4
params%proba_ref(2)=0.6
params%proba_ref(3)=0.8
params%kmargin(1)=3.7375
params%kmargin(2)=2.3288
params%kmargin(3)=1.5294
params%pfalarm   =0.05

! According to the form of the distribution, if white: 0 is enough. 
! If gaussian, it has to be higher according to false alarm rate (refers previously) (without unity) [select_raies]
params%accent=1

! Criteria border for the equalization of the spectrum [fct_multi_courant]
params%ksup=1.2
params%kinf=1./params%ksup

! Saturation of the gain to compensate the effect of the high-pass filter. Used only if params%methselect=3 [mod_par]
params%satur_gain=0.2

! Test if there is a velocity mask or not [mod_par]
params%test_mask=.FALSE.

return
end subroutine advanced_parameters


subroutine header_parameters(params)

! The aim of the this function is to read parameters from the (U)SORT data files.

! Initializations
implicit none

! Parameters structure
type(parameters)                             ::params
! Header data
character,dimension(512)::header
character(132)                ::int_char


! Then reads the header that contains information
open(ind_file_file,file=trim(params%list_files(1)),access='direct', form = 'unformatted',recl=4*128)
read(ind_file_file,rec=1)header
close(ind_file_file)

! Number of chirps
int_char=header(1)//header(2)//header(3)//header(4)
read(int_char,*)params%n_chirp_file

! Date
params%date_str=header(16)//header(17)//header(18)//header(19)&
     //header(20)//header(21)//header(22)//header(23)//header(24)

! Time
params%heure_str=header(26)//header(27)//header(28)//header(29)&
     //header(30)

! Name of the campaign
params%nom_camp=header(37)//header(38)//header(39)//header(40)&
     //header(41)//header(42)//header(43)//header(44)//header(45)//header(46)//header(47)

! Code of the type of data in the header
params%fm_type=header(49)//header(50)//header(51)//header(52)&
     //header(53)//header(54)//header(55)

! Frequency (MHz)
int_char=header(66)//header(67)//header(68)//header(69)&
     //header(70)//header(71)
read(int_char,*)params%FE

! Frequency (Hz)
params%F=params%FE*1.e6

! Year
params%annee=header(81)//header(82)//header(83)//header(84)

! Range resolution (km)
int_char=header(98)//header(99)//header(100)//header(101)&
     //header(102)//header(103)
read(int_char,*)params%resolution

! Angular position of the position for the Rx linear array from antenna 0 to the last one (degree)
! Measure of the cap with the compass: 0=north, 90=east, etc..... [mod_par, mod_ant]
int_char=header(118)//header(119)//header(120)//header(121)
read(int_char,*)params%cap_bous
!params%cap_bous=340.

! Chirp duration (s)
int_char=header(133)//header(134)//header(135)//header(136)&
     //header(137)//header(138)//header(139)
read(int_char,*)params%rate

! Nombre de cellules distance
int_char=header(151)//header(152)//header(153)

if(int_char .EQ. '    ')then
     int_char=header(159)//header(160)//header(161)
     read(int_char,*)params%n_range
else
     read(int_char,*)params%n_range
endif

! Longitude position radar
params%longitude=header(171)//header(172)//header(173)//header(174)&
     //header(175)//header(176)//header(177)//header(178)

! Latitude position of the radar
params%latitude=header(189)//header(190)//header(191)//header(192)&
     //header(193)//header(194)//header(195)//header(196)

! Nb of samples per chirp
int_char=header(204)//header(205)//header(206)//header(207)//header(208)
read(int_char,*)params%MT

! Transmitted power
int_char=header(215)//header(216)//header(217)//header(218)//header(219)
read(int_char,*)params%pwr

! Number of antennas
int_char=header(224)//header(225)
read(int_char,*)params%n_ant

! The "mode bits word" MD is hexadecimally coded
int_char=header(231)//header(232)
!read(int_char,*)params%md

! Offset
int_char=header(241)//header(242)//header(243)//header(244)
read(int_char,*)params%offsett

! Rxoffset
int_char=header(255)//header(256)//header(257)//header(258)//header(259)&
     //header(260)//header(261)//header(262)//header(263)
read(int_char,*)params%rxoffset

! Comment
!params%comment=header[288:384]

return
end subroutine header_parameters


subroutine computed_parameters(params)

! The aim of the this function is to compute parameters needed for the simulation.

! Initializations
implicit none

! Parameters structure
type(parameters)                    ::params
! Miscealenous parameters
integer                                    ::len_porte
real(kind=8)                            ::R_coupure
integer                                    ::i,icouronne,ierr
integer                                    ::ind,ii,jj
integer                                    ::ind_pos,ind_neg,ind_all
logical                                     ::test_var
complex                                  ::j

! Computed parameters
integer                                                    ::pas_radial       ! Radial step (if gliding mode or not)
real(kind=8)                                           ::dkm                !  Equivalent additional delta range for the cables (km)
integer                                                    ::npzero           ! Number of the first crown to consider. Be careful: +1 for null distance
integer                                                    ::npmax            ! Number of the last crown to consider
integer                                                    ::nb_couronne   ! Effective number of crowns
real(kind=8)                                           ::decalage         ! If self.reg_radial is odd, the crown distance is the central gate one
real(kind=8)                                           ::coef                ! Coefficient of conversion doppler frequency-->speed
real(kind=8)                                           ::ChirpFreq        ! Coefficient of conversion doppler frequency-->speed
integer,allocatable,dimension(:)              ::frindex            ! Ordered numbers with a step in the middle
integer,allocatable,dimension(:)              ::fra_signe         ! Sign chart of frindex, used after for speed
real(kind=8),allocatable,dimension(:,:)   ::fraz2D            ! Table to take into account the correction of altitude
real(kind=8),allocatable,dimension(:)     ::coefzD           ! Table to take into account the correction of altitude
real(kind=8)                                           ::R_icouronne    ! Radial position of the crown
real(kind=8)                                           ::distance_icouronne    ! To have the length of the path in the air, minus the equivalent length of the cables
real(kind=8)                                           ::sin_theta,racsin,FB_icouronne    ! Altitude of the antenna taken into account
real(kind=8),allocatable,dimension(:)     ::fraz_icouronne     ! Intermediary table
real(kind=8)                                           ::a                          ! Computation for the relative gain
integer                                                    ::zero_largeur, demi_zero                    ! Interdiction to search lines near the zero line
integer                                                    ::zero_wide_largeur, demi_zero_wide  ! Width larger than the forbidden area
integer,allocatable,dimension(:)              ::i_demi_zero, i_demi_zero_wide          ! Intermediary table
real(kind=8)                                           ::fmax,delta2                          ! Selection of the lines to evaluate the level of noise outoff the sea spectrum
integer,allocatable,dimension(:)              ::int_seuil,liste_raies_bruit_pos,liste_raies_bruit_neg,liste_raies_bruit  ! Intermediary table
real(kind=8)                                           ::azim_ini_hors_jeu,iazim_ini_hors_jeu ! Calculation of the bench area

! Initialization
j=cmplx(0,1)

! Freq_ant
allocate(params%freq_ant(params%n_ant,2),stat=ierr)
params%freq_ant(:,1)=12.3
params%freq_ant(:,2)=0.5

! Total number of chirps
!params%n_chirp=params%n_chirp_file*nlines
params%n_chirp=params%n_chirp_file*params%n_files

! Definition of the angular sector for display (degree, "trigonometry reference") 
params%tetamax=(180.-params%cap_bous)+60.
params%tetamin=(180.-params%cap_bous)-60.

! Wavelength (meter)
params%Lambda=params%c_light/params%F

! Determination of the direction of the sweep of the chirp
!params%slave_state=(((params%md)% 2**8) & ((params%slave_bit)% 2**8 ))
!params%chirp_direction=1-(params%slave_state==params%slave_bit)

! Computation of the number gathered crowns (at least 1)
params%reg_radial=int(max(floor(params%delta_rayon/params%resolution),1))

! Exact depth of the reconstituted radial crowns (in km)
params%delta_rayon=params%reg_radial*params%resolution

! Radial step (if gliding mode or not)
if(params%modeg)then
     pas_radial=1
else
     pas_radial=params%reg_radial
endif

! Equivalent additional delta range for the cables (km)
dkm=((params%lcable*params%ncable)/2.)/1000.

! Number of the first crown to consider. Be careful: +1 for null distance
npzero=int(nint((params%distmin+dkm)/params%resolution)+1.)
if(npzero<=floor(params%reg_radial/2.))then
     npzero=int(floor(params%reg_radial/2.)+1.)
endif

! Number of the last crown to consider
npmax=int(nint((params%distmax+dkm)/params%resolution)+1.)
if(npmax >= (params%n_range-floor(params%reg_radial/2.)))then
     npmax=params%n_range-ceiling(params%reg_radial/2.)
endif

! Effective number of crowns
nb_couronne=int((npmax-npzero)/pas_radial)+1

! Modification on npmax 
npmax=npzero+(pas_radial*(nb_couronne-1))

! Number of gathered crowns to consider
len_porte=npmax-npzero+1
allocate(params%porte(len_porte),stat=ierr)
params%porte=(/(i, i=npzero,npmax, pas_radial)/)

!  Number of studied crowns
params%prof=size(params%porte)

! Range of the radial distances "as the crow flies" , depends on crowns
allocate(params%hdist(len_porte),stat=ierr)
params%hdist=((params%porte-1)*params%resolution)-dkm

! If params%reg_radial is odd, the crown distance is the central gate one
! If params%reg_radial is par, the distance is shifted from a half-crown
decalage=mod((params%reg_radial-1),2)*0.5*params%resolution

! New range of shifted radial distances
params%hdist=params%hdist-decalage

! Modified distances to the horizontal surface to take in account the altitude params%z_ant of the station
params%hdist=sqrt((params%hdist*params%hdist)-((params%z_ant/1000.)*(params%z_ant/1000.)))

! Number of chirps in a subserie (vacation/vacadiv)
params%NPTS=int(floor((params%n_chirp/dble(params%vaca_div))))

! Number of subseries of measure resulting for a file
! Number of subseries of chirps= "incoherent summation"
params%NSUB=(params%Kchevauch*(params%vaca_div-1))+1

! Maximum number of the positive line
params%raie_pos_max=floor(params%NPTS/2.)

! Minimum number of the negative line
params%raie_neg_min=params%raie_pos_max+1

! Frequential resolution DeltaF = 1/Integration time 
params%DeltaF=1./(params%rate*params%NPTS)

! Bragg frequency in Hz
params%FB=0.1020*sqrt(params%FE)

! Bragg line for z_ant=0
! Ratio params%FB on the step of frequential analysis
params%DIbragg=nint(params%FB/params%DeltaF)

! Number of the Bragg line (two number ordered from min to max)
params%Ind_bragg(1)=params%DIbragg+1
params%Ind_bragg(2)=(params%NPTS)-params%DIbragg+1

! Coefficient of conversion doppler frequency-->speed (?)
coef=(params%c_light/(2*params%FE*1e6))*params%FB

! Used only for method 3 of selection of lines
if(params%methselect==3)then
     !Coeff in dB applicable according the distance to the level of noise
     !variable because of the high-pass filter of the receptor and the smoothing filter
     !low-pass used by cleaning
     !has to be reviewed : different of the observed behaviour (base of the holder)
     allocate(params%lcoefnoise(params%prof),stat=ierr)

     ! Distance R_coupure such as the gain is 3 dB under the max gain, so the doppler noise level is 3dB lower than its max
     R_coupure=params%fclow*params%resolution*params%rate
endif

! Chirp frequency (1/Duration of a chirp)
ChirpFreq=1./(params%rate)

! Ordered numbers
allocate(frindex(params%NPTS),stat=ierr)
frindex=(/(i, i=1,params%NPTS, 1)/)

! Ordered number but with a step in the middle
! Number of doppler frequency sign
! frindex is negative for the high ranks=> calcul modulo
frindex(params%raie_neg_min:params%NPTS)=(/(i, i=params%raie_neg_min,params%NPTS, 1)/)-params%NPTS
frindex=frindex-1

! Sign chart of frindex, used after for speed
allocate(fra_signe (params%NPTS),stat=ierr)
fra_signe =sign(1,frindex)
fra_signe(1)=1

! Table to take into account the correction of altitude
allocate(fraz2D(params%prof,params%NPTS),stat=ierr)
allocate(coefzD(params%prof),stat=ierr)
allocate(params%vradz2D(params%prof,params%NPTS),stat=ierr)
allocate(fraz_icouronne(params%NPTS),stat=ierr)

! Loop on the crown to calculate the associated radial speed to each line of the sprectrum
do icouronne=1, params%prof
     ! Radial position of the crown
     R_icouronne=(params%porte(icouronne)-1)*params%resolution

     ! To have the length of the path in the air, minus the equivalent length of the cables
     distance_icouronne=R_icouronne-dkm

     ! Altitude of the antenna taken into account
     sin_theta=sqrt(1./(1.+(((params%z_ant/1000.)/distance_icouronne)**2.)))
     racsin=sqrt(sin_theta)
     FB_icouronne=params%FB*racsin

     ! Intermediary table
     fraz_icouronne=(ChirpFreq/FB_icouronne)*(frindex/real(params%NPTS,kind=8))

     ! Always to take the altitude into account
     ! DIbragg and Indbragg  also depend on distance and should be generated as table here
     fraz2D(icouronne,:) =fraz_icouronne-fra_signe
     coefzD(icouronne)   =coef/racsin
     params%vradz2D(icouronne,:)=(coef/racsin)*(fraz_icouronne-fra_signe )
          
     ! Used only for method 3 of line selection
     if(params%methselect==3)then
         !Calculation for the relative gain
         a=abs(R_icouronne/(R_icouronne+(j*R_coupure)))
         
         !Saturation to eliminate the effect of the high frequency filtering
         if(a<params%satur_gain)then
             a=params%satur_gain
             params%lcoefnoise(icouronne)=20.*log10(a)
         endif
     endif
enddo


! Interdiction to search lines near the zero line
! Width of the forbidden area
zero_largeur=nint(params%NPTS*params%kzero)
demi_zero   =nint(zero_largeur/2.)
allocate(i_demi_zero(2*demi_zero+1),stat=ierr)
i_demi_zero =(/(i, i=-demi_zero,demi_zero, 1)/)
params%demi_zero=demi_zero
! Numbers of the forbidden area
allocate(params%liste_raies_zero(2*demi_zero+1),stat=ierr)
params%liste_raies_zero=modulo(i_demi_zero,params%NPTS)+1

! Width larger than the forbidden area
zero_wide_largeur=nint(params%NPTS*params%kzero_wide)
demi_zero_wide=nint(zero_wide_largeur/2.)
allocate(i_demi_zero_wide(2*demi_zero_wide+1),stat=ierr)
i_demi_zero_wide =(/(i, i=-demi_zero_wide,demi_zero_wide, 1)/)

! Number of the forbidden area
allocate(params%liste_raies_zero_wide(2*demi_zero_wide+1),stat=ierr)
params%liste_raies_zero_wide=modulo(i_demi_zero_wide,params%NPTS)+1

! Selection of the lines to evaluate the level of noise outoff the sea spectrum i.e. the Bragg peaks (some precautions are taken to not to work on the second order peak)
! definition according vmax of the exclusion area of the lines for the research of sources around the Bragg lines
! fmax doppler corresponds to the defined max speed
fmax=params%vmax*(2*params%F/params%c_light)

! Minimum space betwenn Bragg lines and lines used to estimate the noise (in spectral points)
delta2=max(((1.4*fmax)/params%DeltaF),(0.5*params%DIbragg));

! Determination of the number of the lines used to estimate the noise
allocate(int_seuil(params%NPTS),stat=ierr)
allocate(liste_raies_bruit_pos(params%NPTS),stat=ierr)
allocate(liste_raies_bruit_neg(params%NPTS),stat=ierr)
allocate(liste_raies_bruit(params%NPTS),stat=ierr)

int_seuil=(/(i, i=1,params%NPTS, 1)/)
liste_raies_bruit_pos=0
liste_raies_bruit_neg=0
liste_raies_bruit=0

ind_pos=0
ind_neg=0
ind_all=0

do ind=1,params%NPTS
     if(abs(int_seuil(ind)-params%Ind_bragg(1))>delta2)then
         ind_pos=ind_pos+1
         liste_raies_bruit_pos(ind_pos)=ind
     endif
enddo

do ind=1,params%NPTS
     if(abs(int_seuil(ind)-params%Ind_bragg(2))>delta2)then
         ind_neg=ind_neg+1
         liste_raies_bruit_neg(ind_neg)=ind
     endif
enddo

do ii=1,ind_pos
     do jj=1,ind_neg
         if(liste_raies_bruit_pos(ii) .EQ. liste_raies_bruit_neg(jj))then
             ind_all=ind_all+1
             liste_raies_bruit(ind_all)=liste_raies_bruit_pos(ii)
             exit
         endif
     enddo
enddo

ind_pos=0
! The zero area of the largest width is not taken into account
do ii=1,ind_all
     test_var=.TRUE.
     do jj=1,2*demi_zero_wide+1
         if(liste_raies_bruit(ii) .EQ. params%liste_raies_zero_wide(jj))then
             test_var=.FALSE.
             exit
         endif
     enddo
     if(test_var)then
         ! Intermediate variable used for test purpose only
         ind_pos=ind_pos+1
         liste_raies_bruit_pos(ind_pos)=liste_raies_bruit(ii)
     endif
enddo

allocate(params%liste_raies_bruit(ind_pos),stat=ierr)
do ind =1,ind_pos
     params%liste_raies_bruit(ind)=liste_raies_bruit_pos(ind)
enddo

! List range for the selection of the sample, base on a level of energy to far distances, without effect of holder base noise retrodiffusion de la mer
allocate(params%liste_range(ceiling((1.-params%n_range_min)*params%n_range)),stat=ierr)
params%liste_range=(/(i, i=ceiling(params%n_range_min*params%n_range),params%n_range, 1)/)

! Declaration of the table of angular sectors , which can depends on the chosen method, according to the angular step dteta
params%n_dir=ceiling((params%tetamax-params%tetamin)/params%dteta)+1
allocate(params%tetad(params%n_dir),stat=ierr)
params%tetad=params%tetamin+(params%dteta*(/(i, i=0,params%n_dir-1, 1)/))

!! Preparation of TVMAX, mask table of the max speed
allocate(params%TVMAX(params%prof,params%n_dir),stat=ierr)
params%TVMAX=params%vmax

!
! Number of thin azimuths [var_ant]
params%nazim=ceiling((params%azimmax-params%azimmin)/params%dazim)+1
! Table of azimuth in degree [var_ant]
allocate(params%azimd(params%nazim),stat=ierr)
params%azimd=params%azimmin+(params%dazim*(/(i, i=0,params%nazim-1, 1)/))

! Calculation of the bench area
if(params%ya_sym)then
     azim_ini_hors_jeu  =modulo((360.-params%cap_bous)+90.+180.,360.)
     iazim_ini_hors_jeu =(azim_ini_hors_jeu-params%azimmin)/params%dazim
     allocate(params%hors_jeu_sym(nint(params%nazim/2.)),stat=ierr)
     params%hors_jeu_sym=modulo(nint(iazim_ini_hors_jeu+(/(i, i=0,nint(params%nazim/2.), 1)/)),params%nazim)+1
endif

return
end subroutine computed_parameters


end module init_params_module
