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

subroutine pass2(qlmax_up,nazim,qlmax)

    implicit none

    ! Parameters
    integer,intent(in) :: nazim
    real(kind=8),intent(out)    :: qlmax(nazim)
    real(kind=8),intent(in)   :: qlmax_up(nazim)

    ! Local variables
    integer ind_azim

    do ind_azim=nazim,1,-1
         if(qlmax_up(ind_azim)<qlmax_up(modulo(ind_azim-2, (nazim))+1))then
            qlmax(ind_azim)=0.
         else
            qlmax(ind_azim)=qlmax_up(ind_azim)
         endif
     enddo

end subroutine pass2
