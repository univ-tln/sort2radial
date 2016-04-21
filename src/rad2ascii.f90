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


subroutine rad2ascii(courant_dir,params)

! The aim of this function is to export radial values to ascii file
! the main parameters of the simulation

! Initializations
! ---------------
use init_params_module
implicit none

! Parameters structure
type(parameters)::params
real(kind=8),dimension(params%prof,params%n_dir)::courant_dir
real(kind=8)::LON_Radar,LAT_Radar,h,m,s,lat,lon
real(kind=8),dimension(4,params%not_nan)::mat_int
integer::ind,ind_az,ind_ran,ind_char,ierr,ref_char,durmin,dursec
real(kind=8)::az_mod,duration,mintetad,maxtetad,minazi,maxazi
character,allocatable,dimension(:)::output_file


read(params%longitude(1:2),*)h
read(params%longitude(4:5),*)m
read(params%longitude(7:8),*)s
LON_Radar=-(h+m/60.+s/3600.)

read(params%latitude(1:2),*)h
read(params%latitude(4:5),*)m
read(params%latitude(7:8),*)s
LAT_Radar=(h+m/60.+s/3600.)
mat_int=0.
ind=0

do ind_az=1,params%n_dir
     az_mod=modulo(450.-params%tetad(ind_az),360.)
     do ind_ran=1,params%prof
         if(courant_dir(ind_ran,ind_az)<100.)then
         ! Conversion cartesian reference
             call radial_LatLon_WGS84_r16(az_mod,params%hdist(ind_ran),LAT_Radar,LON_Radar,lat,lon)
             ind=ind+1
             mat_int(1, ind)=lat
             mat_int(2, ind)=lon
             mat_int(3, ind)=courant_dir(ind_ran,ind_az)
             mat_int(4, ind)=az_mod
         endif
     enddo
enddo


! Output
ref_char=0
do ind_char=len(trim(params%output_file)),1,-1
     if(params%output_file(ind_char:ind_char) .NE. '/')then
         ref_char=ref_char+1
     else
         exit
     endif
enddo
allocate(output_file(ref_char),stat=ierr)
do ind_char=ref_char,1,-1
     output_file(ref_char-ind_char+1:ref_char-ind_char+1)=params%output_file(&
          len(trim(params%output_file))-ind_char+1:len(trim(params%output_file))-ind_char+1)
enddo

duration=params%n_chirp*params%rate

durmin=floor(duration/60.)
dursec=nint(modulo(duration,60.))

!ajout temporaire export radial
ind=0
open(ind_file_file,file=trim(params%output_file),action="write",status="replace")
write(ind_file_file,*),"Radial surface current velocities"
write(ind_file_file,'(A,2F10.6)'),trim(params%site_name),LON_Radar,LAT_Radar
write(ind_file_file,*),trim(params%date_str)," ", trim(params%heure_str)
write(ind_file_file,'(2I5)'),durmin,dursec
write(ind_file_file,'(2F12.6,1I5)'),minval(params%hdist(:)),maxval(params%hdist(:)),params%prof
mintetad=minval(params%tetad(:))
maxtetad=maxval(params%tetad(:))
minazi=modulo(450.-maxtetad,360.)
maxazi=minazi+(maxtetad-mintetad)

write(ind_file_file,'(2F12.6,1I5)'),minazi,maxazi,params%n_dir
write(ind_file_file,*),"Iran  Iazi  R(Km)     Azi(deg.) Lat(deg.)     Lon(deg.)     Vel(m/s)"

do ind_az=1,params%n_dir
     az_mod=modulo(450.-params%tetad(ind_az),360.)
     do ind_ran=1,params%prof
         if(courant_dir(ind_ran,ind_az)<100.)then
             ind=ind+1
             write(ind_file_file,'(2I6,1F11.5,1F9.3,2F14.8,1F7.3)'), ind_ran, &
                  params%n_dir-ind_az+1, params%hdist(ind_ran), mat_int(4, ind),&
                  mat_int(1, ind), mat_int(2, ind), mat_int(3, ind)
         endif
     enddo
enddo
close(ind_file_file)


return
end subroutine rad2ascii
