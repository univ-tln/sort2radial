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


! Compute inverse of complex matrix

      subroutine cinv( A, N, MAXN, Ainv )
      integer*4 N, MAXN
      complex*16 A(MAXN,MAXN), Ainv(MAXN,MAXN)
! Inputs
!   A       Matrix A to be inverted
!   N       Elements used in matrix A (N by N)
!  MAXN     Matrix dimenstions as A(MAXN,MAXN)
! Outputs
!  Ainv     Inverse of matrix A

      integer*4 MAXMAXN
      parameter( MAXMAXN = 200 )
      integer*4 i, j, k, index(MAXMAXN), jPivot, indexJ
      real*8 scale(MAXMAXN), scaleMax, ratio, ratioMax
      complex*16 AA(MAXMAXN,MAXMAXN), B(MAXMAXN,MAXMAXN), coeff, sum

      if( MAXN .gt. MAXMAXN ) then
        write(*,*) 'ERROR in cinv: Matrix too large'
        stop
      endif

      !* Matrix B is initialized to the identity matrix
      do i=1,N
       do j=1,N
         AA(i,j) = A(i,j)  ! Copy matrix so as not to overwrite
         B(i,j) = 0.0
       enddo
       B(i,i) = 1.0
      enddo

      !* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
      do i=1,N
        index(i) = i     ! Initialize row index list
        scaleMax = 0.0
        do j=1,N
          if( abs(AA(i,j)) .gt. scaleMax ) then
            scaleMax = abs(AA(i,j))
          endif
        enddo
        scale(i) = scaleMax
      enddo

      !* Loop over rows k = 1, ..., (N-1)
      do k=1,(N-1)
        !* Select pivot row from max( |a(j,k)/s(j)| )
        ratiomax = 0.0
        jPivot = k
        do i=k,N
          ratio = abs(AA(index(i),k))/scale(index(i))
          if( ratio .gt. ratiomax ) then
            jPivot=i
            ratiomax = ratio
          endif
        enddo
        !* Perform pivoting using row index list
        indexJ = index(k)
        if( jPivot .ne. k ) then     ! Pivot
          indexJ = index(jPivot)
          index(jPivot) = index(k)   ! Swap index jPivot and k
          index(k) = indexJ
        endif
        !* Perform forward elimination
        do i=k+1,N
          coeff = AA(index(i),k)/AA(indexJ,k)
          do j=k+1,N
            AA(index(i),j) = AA(index(i),j) - coeff*AA(indexJ,j)
          enddo
          AA(index(i),k) = coeff
          do j=1,N
            B(index(i),j) = B(index(i),j) - AA(index(i),k)*B(indexJ,j)
          enddo
        enddo
      enddo

      !* Perform backsubstitution
      do k=1,N
        Ainv(N,k) = B(index(N),k)/AA(index(N),N)
        do i=N-1,1,-1
          sum = B(index(i),k)
          do j=i+1,N
            sum = sum - AA(index(i),j)*Ainv(j,k)
          enddo
          Ainv(i,k) = sum/AA(index(i),i)
        enddo
      enddo

      return
      end
