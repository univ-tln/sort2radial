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

subroutine smooth(qf,inter_liss,n_liss,nazim,ql)

     implicit none

     ! Parameters
     integer:: n_liss,nazim
     integer:: inter_liss(n_liss)
     real(kind=8) :: qf(nazim)
     real(kind=8) :: ql(nazim)

     ! Local variables
     integer a_inter
     integer ind_azim,ind_liss
     ql=0.
     !  Boucle sur la taille du vecteur facteur music
     do ind_azim = 1, nazim
         !print*,ind_azim,n_liss
         ! Boucle sur la taille du vecteur filtre
         do ind_liss = 1,n_liss
             !print*,ind_azim,ind_liss
             a_inter=ind_azim+inter_liss(ind_liss)
             if(a_inter > nazim)then
                 a_inter=a_inter-(nazim)
             endif
             if(a_inter < 1)then
                 a_inter=(nazim)+a_inter
             endif

             ql(ind_azim)=ql(ind_azim)+qf(a_inter)
         enddo
         ql(ind_azim)=ql(ind_azim)/n_liss
      enddo

end subroutine smooth
