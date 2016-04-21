!***********************************************************************
!*     Calculate eigenvalues and eigenvectors of a Square Hermitian    *
!*     Matrix By Jacobi's Method                                       *
!* ------------------------------------------------------------------- *
!* SAMPLE RUN:                                                         *
!* File Matc4.dat contains:                                            *
!* 4                                                                   *
!* 12   0   -1   -1   2   0   3   3                                    *
!* -1   1   12    0   1  -1  -2   0                                    *
!*  2   0    1    1   8   0  -1  -1                                    *
!*  3  -3   -2    0  -1   1   8   0                                    *
!*                                                                     *
!*  Precision: 1d-10                                                   *
!*  Max. number of iterations: 25                                      *
!*                                                                     *
!*  Eigenvalue  1:   0.400000E+01                                      *
!*  Eigenvector:                                                       *
!* -0.500000E+00 -  0.500000E+00 I                                     *
!* -0.624817E-17 +  0.000000E+00 I                                     *
!*  0.500000E+00 +  0.500000E+00 I                                     *
!*  0.100000E+01 +  0.000000E+00 I                                     *
!*                                                                     *
!*  Eigenvalue  2:   0.800000E+01                                      *
!*  Eigenvector:                                                       *
!* -0.751314E-17 +  0.000000E+00 I                                     *
!* -0.500000E+00 +  0.500000E+00 I                                     *
!*  0.100000E+01 +  0.000000E+00 I                                     *
!* -0.500000E+00 +  0.500000E+00 I                                     *
!*                                                                     *
!*  Eigenvalue  3:   0.120000E+02                                      *
!*  Eigenvector:                                                       *
!*  0.500000E+00 +  0.500000E+00 I                                     *
!*  0.100000E+01 +  0.000000E+00 I                                     *
!*  0.500000E+00 +  0.500000E+00 I                                     *
!*  0.220564E-15 +  0.000000E+00 I                                     *
!*                                                                     *
!*  Eigenvalue  4:   0.160000E+02                                      *
!*  Eigenvector:                                                       *
!*  0.100000E+01 +  0.000000E+00 I                                     *
!* -0.500000E+00 +  0.500000E+00 I                                     *
!* -0.553902E-16 +  0.000000E+00 I                                     *
!*  0.500000E+00 -  0.500000E+00 I                                     *
!*                                                                     *
!* ------------------------------------------------------------------- *
!* Reference: "ALGEBRE Algorithmes et programmes en Pascal             *
!*             By Jean-Louis Jardrin - Dunod BO-PRE 1988" [BIBLI 10]   *
!*                                                                     *
!*                              F90 Release 1.1 By J-P Moreau, Paris.  *
!*                                                                     *
!* Release 1.1 (Nov. 30th, 2007): Added function NORMAL.               * 
!***********************************************************************                                       
!Program TEPHJ
module vp_module

implicit none

integer, parameter :: SIZE = 25

! complex number
TYPE Zcomplex
   REAL*8  R	  ! algebraic form
   REAL*8  I	  !
 END TYPE Zcomplex
  
 integer I,it,J,M,N
 real*8 dta
 Type(Zcomplex) A(SIZE,SIZE), VX(SIZE,SIZE)
 Real*8 R(SIZE)
  
  contains

!  call LME(N,A)      !read complex matrix from text file
!
!  print *,' '
!  print *,'  ** Compute Eigenvalues and Eigenvectors **'
!  print *,'         of a Square Hermitian Matrix'
!  print *,'             By Jacobi''s Method'
!  print *,' '
!  write(*,10,advance='no'); read *, dta
!  write(*,20,advance='no'); read *, M
!  print *,' '
!
!  call EPHJ(dta, M, N, A, it, R, VX)
!  print *,'herhe',M,it
!  call NORMAL(N, R, VX)
!
!  if (it==-1) then
!    print *,'  No convergence%'
!  else                             !display results
!	do J=1, N
!      write(*,30) J, R(J)
!	  print *,'  Eigenvector:'
!	  do I=1, N
!        write(*,40,advance='no') VX(I,J)%R
!        if (VX(I,J)%I >= 0.d0) then 
!		  write(*,41,advance='no')
!        else 
!		  write(*,42,advance='no')
!        end if
!		write(*,50) dabs(VX(I,J)%I)
!      end do
!	  print *,' '
!    end do
!	print *,' '
!  end if
!  stop
!
!10 format('  Precision: ')
!20 format('  Max. number of iterations: ')
!30 format('  Eigenvalue ',I2,':  ',E13.6)
!40 format('  ',E13.6)
!41 format(' + ')
!42 format(' - ')
!50 format(E13.6,' I')  
!
!END  !of main program

  real*8 Function Sqr(a)
    real*8 a
    Sqr = a*a
    return
  end Function Sqr

  !absolute value of a complex number
  real*8 Function CABS1(C)
    type Zcomplex
      real*8  R
      real*8  I
    end type Zcomplex
    type(Zcomplex) C
	real*8 Sqr
    CABS1 = dsqrt(Sqr(C%R)+Sqr(C%I))
	return
  end Function CABS1

  !Add two complex numbers
  Subroutine CADD(c1, c2, c3)
    type Zcomplex
      real*8  R
      real*8  I
    end type Zcomplex
    type(Zcomplex) c1, c2, c3
    c3%R=c1%R+c2%R; c3%I=c1%I+c2%I
	return
  end Subroutine CADD

  !Substract two complex numbers
  Subroutine CDIF(c1, c2, c3)
    type Zcomplex
      real*8  R
      real*8  I
    end type Zcomplex
    type(Zcomplex) c1, c2, c3
    c3%R=c1%R-c2%R; c3%I=c1%I-c2%I
	return
  end Subroutine CDIF

  !Multiply two complex numbers
  Subroutine CMUL(c1, c2, c3)
    type Zcomplex
      real*8  R
      real*8  I
    end type Zcomplex
    type(Zcomplex) c1, c2, c3
    c3%R=c1%R*c2%R - c1%I*c2%I
    c3%I=c1%R*c2%I + c1%I*c2%R
	return
  end Subroutine CMUL  

  !Return conjugate of a complex number
  Subroutine CONJ(c, c1)
    type Zcomplex
      real*8  R
      real*8  I
    end type Zcomplex
    type(Zcomplex) c, c1
    c1%R=c%R; c1%I=-c%I
	return
  end Subroutine CONJ

  !Multiply a complex number by a real number
  Subroutine CPRO(alpha, c, c1)
    type Zcomplex
      real*8  R
      real*8  I
    end type Zcomplex
    real*8 alpha; type(Zcomplex) c, c1
    c1%R=alpha*c%R; c1%I=alpha*c%I
	return
  end Subroutine CPRO


!************************************************************************
!*   Compute all eigenvalues/eigenvectors of a square hermitian matrix  *
!*   using Jacobi's method.                                             *
!* -------------------------------------------------------------------- *
!* Inputs:                                                              *
!*         dta: required precision (double)                             *
!*         M  : maximum number of iterations (integer)                  *
!*         N  : size of matrix                                          *
!*         A  : pointer to hermitian (complex) matrix                   *
!* Outputs:                                                             *
!*         it : flag for convergence (integer) =-1, no convergence      *
!*         R  : vector(1%%N) of eigenvalues (double)                    * 
!*         VX : matrix storing complex eigenvectors in columns          *
!************************************************************************
Subroutine EPHJ(dta, M, N, A, it, R, VX)
  integer, parameter :: SIZE = 25
  type Zcomplex
    real*8  R
    real*8  I
  end type Zcomplex
  real*8 dta, R(SIZE), Sqr, CABS1
  type(Zcomplex) A(SIZE,SIZE), VX(SIZE,SIZE)
  integer I,J,K,K0,L,L0
  real*8 delta,s,s0,t0,t1,w0
  type(Zcomplex) c0,c1,c2,c3,u0,u1,z0,z1

  z0%R=0.d0; z0%I=0.d0; z1%R=1.d0; z1%I=0.d0
  do I=1, N
    do J=1, N
      if (I==J) then
	    VX(I,J)=z1
	  else 
	    VX(I,J)=z0
      end if
    end do
  end do
  it=-1; L=1
  do while (L<=M.and.it.ne.1)
    s=0.d0
    do I=1, N
	  do J=I+1, N
        t0=CABS1(A(I,J))
	    if (t0 > s) then
          s=t0; K0=I; L0=J
		end if
      end do
    end do
    if (s==0.d0) then
	  it=1
    else
      delta=Sqr(A(L0,L0)%R-A(K0,K0)%R)+4.d0*Sqr(CABS1(A(K0,L0)))
      t0=A(L0,L0)%R-A(K0,K0)%R + dsqrt(delta)
      t1=A(L0,L0)%R-A(K0,K0)%R - dsqrt(delta)
      if (dabs(t0) >= dabs(t1)) then
	    w0=t0
	  else 
	    w0=t1
      end if
      s0=dabs(w0)/dsqrt(Sqr(w0)+4.d0*Sqr(CABS1(A(K0,L0))))
      t0=2.d0*s0/w0; call CPRO(t0,A(K0,L0),c0)
      call CONJ(c0,c1)
      do I=1, K0-1
        u0=A(I,K0);
        call CMUL(c0,u0,c2); call CPRO(s0,A(I,L0),c3); call CADD(c2,c3,A(I,K0))
        call CMUL(c1,A(I,L0),c2); call CPRO(s0,u0,c3); call CDIF(c2,c3,A(I,L0))
      end do
      do K=K0+1, L0-1
        u0=A(K0,K)
        call CMUL(c1,u0,c2); call CONJ(A(K,L0),u1); call CPRO(s0,u1,c3)
        call CADD(c2,c3,A(K0,K))
        call CMUL(c1,A(K,L0),c2); call CONJ(u0,u1); call CPRO(s0,u1,c3)
        call CDIF(c2,c3,A(K,L0))
      end do
      do J=L0+1, N
        u0=A(K0,J)
        call CMUL(c1,u0,c2); call CPRO(s0,A(L0,J),c3); call CADD(c2,c3,A(K0,J))
        call CMUL(c0,A(L0,J),c2); call CPRO(s0,u0,c3); call CDIF(c2,c3,A(L0,J))
      end do
      t0=A(K0,K0)%R
      t1=4.d0*Sqr(s0*CABS1(A(K0,L0)))/w0
      A(K0,K0)%R=Sqr(CABS1(c0))*t0 + t1+Sqr(s0)*A(L0,L0)%R
      A(L0,L0)%R=Sqr(s0)*t0 - t1+Sqr(CABS1(c0))*A(L0,L0)%R
      A(K0,L0)=z0
      do I=1, N
        u0=VX(I,K0)
        call CMUL(c0,u0,c2); call CPRO(s0,VX(I,L0),c3); call CADD(c2,c3,VX(I,K0))
        call CMUL(c1,VX(I,L0),c2); call CPRO(s0,u0,c3); call CDIF(c2,c3,VX(I,L0))
      end do
      t0=0.d0
      do I=1, N-1
        do J=I+1, N
          t0 = t0 + Sqr(CABS1(A(I,J)))
        end do
      end do
      s=2.d0*t0
      if (s < dta) then
	    it=1
	  else 
	    L=L+1
      end if
	end if
  end do
  if (it==1) then
    do I=1, N
	  R(I)=A(I,I)%R
    end do
  end if

  return
End Subroutine EPHJ

Subroutine CDIV(AR, AI, BR, BI, ZR, ZI)
!      EFFECT DIVISION Z=(ZR+I*ZI)=(AR+I*AI)/(BR+I*BI)
! -----DO NOT USE IF BR=BI=0
  real*8 AR,AI,BR,BI,ZR,ZI
  real*8 YR,YI,W
  YR=BR
  YI=BI
  if (dabs(YR) > dabs(YI)) then
    W=YI/YR
    YR=W*YI+YR
    ZR=(AR+W*AI)/YR
    ZI=(AI-W*AR)/YR
  else
    W=YR/YI
    YI=W*YR+YI
    ZR=(W*AR+AI)/YI
    ZI=(W*AI-AR)/YI
  end if
  return
End Subroutine CDIV

Subroutine NORMAL(N, R, Z)
!--------------------------------------------------------------
!   Sort eigenvalues by fabsolute values in ascending order and
!   normalize.
!--------------------------------------------------------------
  integer, parameter :: SIZE = 25
  type Zcomplex
    real*8  R
    real*8  I
  end type Zcomplex
  type(Zcomplex) Z(SIZE,SIZE)
  real*8 R(SIZE)
  real*8 Tr(SIZE),Ti(SIZE)
  real*8 ZM,Z1,Z2,VR, CABS1
  type(Zcomplex) V
  integer IM,J,K

!   SORT SOLUTIONS IN ASCENDING ORDER
    do J=2, N
      VR=R(J)
      do K=1, N
        Tr(K)=Z(K,J)%R
        Ti(K)=Z(K,J)%I
      end do

      do I=J-1, 1, -1
        if (dabs(R(I)) <= dabs(VR)) goto 5
        R(I+1)=R(I)
        do K=1, N
          Z(K,I+1)%R=Z(K,I)%R
          Z(K,I+1)%I=Z(K,I)%I
        end do
	  end do
      I=0
5     R(I+1)=VR
      do K=1, N
        Z(K,I+1)%R=Tr(K)
        Z(K,I+1)%I=Ti(K)
      end do
    end do

!   NORMALIZE WITH RESPECT TO BIGGEST ELEMENT
	do J=N, 1, -1
      ZM = 0.d0
      do I=1, N
        V%R=Z(I,J)%R
        V%I=Z(I,J)%I
        Z1=CABS1(V)
        if (Z1 >= ZM) then
          IM = I
          ZM = Z1
        end if
      end do
      Z1 = Z(IM,J)%R
      Z2 = Z(IM,J)%I
      do I=1, N
        call CDIV(Z(I,J)%R,Z(I,J)%I,Z1,Z2,V%R,V%I)
        Z(I,J)%R = V%R
        Z(I,J)%I = V%I
      end do
    end do
    return
End Subroutine NORMAL

! end of file Tephj.f90
end module vp_module
