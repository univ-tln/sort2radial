subroutine bakvec ( nm, n, t, e, m, z, ierr )
!
!*******************************************************************************
!
!! BAKVEC determines eigenvectors by reversing the FIGI transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a nonsymmetric tridiagonal 
!    matrix by back transforming those of the corresponding symmetric 
!    matrix determined by FIGI.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of T and Z, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) T(NM,3), contains the nonsymmetric matrix.  Its subdiagonal 
!    is stored in the positions 2:N of the first column,
!    its diagonal in positions 1:N of the second column,
!    and its superdiagonal in positions 1:N-1 of the third column.  
!    T(1,1) and T(N,3) are arbitrary.
!
!    Input/output, real(kind=8) E(N).  On input, E(2:N) contains the subdiagonal 
!    elements of the symmetric matrix.  E(1) is arbitrary.
!    On output, the contents of E have been destroyed.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) Z(NM,M), contains the eigenvectors.  On output,
!    they have been transformed as requested.
!
!    Output, integer IERR, an error flag.
!    0, for normal return,
!    2*N+I, if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
!    In this case, the symmetric matrix is not similar
!    to the original matrix, and the eigenvectors
!    cannot be found by this program.
!
  implicit none
!
  integer m
  integer n
  integer nm
!
  real(kind=8) e(n)
  integer i
  integer ierr
  integer j
  real(kind=8) t(nm,3)
  real(kind=8) z(nm,m)
!
  ierr = 0

  if ( m == 0 ) then
    return
  end if

  e(1) = real(1.0E+00,kind=8)
  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    if ( e(i) == real(0.0E+00,kind=8) ) then
      if ( t(i,1) /= real(0.0E+00,kind=8) .or. t(i-1,3) /= real(0.0E+00,kind=8) ) then
        ierr = 2 * n + i
        return
      end if
      e(i) = real(1.0E+00,kind=8)
    else
      e(i) = e(i-1) * e(i) / t(i-1,3)
    end if
  end do

  do j = 1, m
    z(2:n,j) = z(2:n,j) * e(2:n)
  end do

  return
end
subroutine balanc ( nm, n, a, low, igh, scale )
!
!*******************************************************************************
!
!! BALANC balances a real(kind=8) matrix before eigenvalue calculations.
!
!
!  Discussion:
!
!    This subroutine balances a real(kind=8) matrix and isolates eigenvalues 
!    whenever possible.
!
!    Suppose that the principal submatrix in rows LOW through IGH
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!
!      SCALE(J) = P(J),    J = 1,...,LOW-1,
!               = D(J,J),  J = LOW,...,IGH,
!               = P(J)     J = IGH+1,...,N.
!
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is returned for LOW if IGH is zero formally.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) A(NM,N), the N by N matrix.  On output,
!    the matrix has been balanced.
!
!    Output, integer LOW, IGH, indicate that A(I,J) is equal to zero if
!    (1) I is greater than J and
!    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!    Output, real(kind=8) SCALE(N), contains information determining the
!    permutations and scaling factors used.
!
  implicit none
!
  integer nm
  integer n
!
  real(kind=8) a(nm,n)
  real(kind=8) b2
  real(kind=8) c
  real(kind=8) f
  real(kind=8) g
  integer i
  integer iexc
  integer igh
  integer j
  integer k
  integer l
  integer low
  integer m
  logical noconv
  real(kind=8) r
  real(kind=8) radix
  real(kind=8) s
  real(kind=8) scale(n)
!
  radix = 16.0E+00

  iexc = 0
  j = 0
  m = 0

  b2 = radix**2
  k = 1
  l = n
  go to 100

20 continue

  scale(m) = j

  if ( j /= m ) then

    do i = 1, l
      call r_swap ( a(i,j), a(i,m) )
    end do

    do i = k, n
      call r_swap ( a(j,i), a(m,i) )
    end do

  end if

50 continue

  if ( iexc == 2 ) go to 130
!
!  Search for rows isolating an eigenvalue and push them down.
!
80 continue

  if ( l == 1 ) then
    low = k
    igh = l
    return
  end if

  l = l - 1

100 continue

  do j = l, 1, -1

     do i = 1, l
       if ( i /= j ) then
         if ( a(j,i) /= real(0.0E+00,kind=8) ) then
           go to 120
         end if
       end if
     end do

     m = l
     iexc = 1
     go to 20

120  continue

  end do

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 continue

  k = k + 1

140 continue

  do j = k, l

    do i = k, l
      if ( i /= j ) then
        if ( a(i,j) /= real(0.0E+00,kind=8) ) then
          go to 170
        end if
      end if
    end do

    m = k
    iexc = 2
    go to 20

170 continue

  end do
!
!  Balance the submatrix in rows K to L.
!
  scale(k:l) = real(1.0E+00,kind=8)
!
!  Iterative loop for norm reduction.
!
  noconv = .true.

  do while ( noconv )

    noconv = .false.

    do i = k, l

      c = real(0.0E+00,kind=8)
      r = real(0.0E+00,kind=8)

      do j = k, l
        if ( j /= i ) then
          c = c + abs ( a(j,i) )
          r = r + abs ( a(i,j) )
        end if
      end do
!
!  Guard against zero C or R due to underflow.
!
      if ( c /= real(0.0E+00,kind=8) .and. r /= real(0.0E+00,kind=8) ) then

        g = r / radix
        f = real(1.0E+00,kind=8)
        s = c + r

        do while ( c < g ) 
          f = f * radix
          c = c * b2
        end do

        g = r * radix

        do while ( c >= g ) 
          f = f / radix
          c = c / b2
        end do
!
!  Balance.
!
        if ( ( c + r ) / f < 0.95E+00 * s ) then

          g = real(1.0E+00,kind=8) / f
          scale(i) = scale(i) * f
          noconv = .true.

          a(i,k:n) = a(i,k:n) * g
          a(1:l,i) = a(1:l,i) * f

        end if

      end if

    end do

  end do

  low = k
  igh = l

  return
end
subroutine balbak ( nm, n, low, igh, scale, m, z )
!
!*******************************************************************************
!
!! BALBAK determines eigenvectors by undoing the BALANC transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a real(kind=8) general matrix by 
!    back transforming those of the corresponding balanced matrix 
!    determined by BALANC.
!
!  Modified:
!
!    18 February 2001
!
!  Reference:
!
!    Parlett and Reinsch,
!    Numerische Mathematik,
!    Volume 13, pages 293-304, 1969.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of Z.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, column indices determined by BALANC.
!
!    Input, real(kind=8) SCALE(N), contains information determining
!    the permutations and scaling factors used by BALANC.
!
!    Input, integer M, the number of columns of Z to be back-transformed.
!
!    Input/output, real(kind=8) Z(LDA,M), contains the real(kind=8) and imaginary parts 
!    of the eigenvectors, which, on return, have been back-transformed.
!
  implicit none
!
  integer m
  integer n
  integer nm
!
  integer i
  integer igh
  integer ii
  integer j
  integer k
  integer low
  real(kind=8) scale(n)
  real(kind=8) z(nm,m)
!
  if ( m <= 0 ) then
    return
  end if

  if ( igh /= low ) then
    do i = low, igh
      z(i,1:m) = z(i,1:m) * scale(i)
    end do
  end if

   do ii = 1, n

     i = ii

     if ( i < low .or. i > igh ) then

       if ( i < low ) then
         i = low - ii
       end if

       k = int ( scale(i) )

       if ( k /= i ) then

         do j = 1, m
           call r_swap ( z(i,j), z(k,j) )
         end do

        end if

      end if

  end do

  return
end
subroutine bandr ( nm, n, mb, a, d, e, e2, matz, z )
!
!*******************************************************************************
!
!! BANDR reduces a symmetric band matrix to symmetric tridiagonal form.
!
!
!  Discussion:
!
!    This subroutine reduces a real(kind=8) symmetric band matrix
!    to a symmetric tridiagonal matrix using and optionally
!    accumulating orthogonal similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MB, is the (half) band width of the matrix, defined as the
!    number of adjacent diagonals, including the principal
!    diagonal, required to specify the non-zero portion of the
!    lower triangle of the matrix.
!
!    Input/output, real(kind=8) A(NM,MB).  On input, contains the lower triangle of 
!    the symmetric band input matrix stored as an N by MB array.  Its lowest 
!    subdiagonal is stored in the last N+1-MB positions of the first column,
!    its next subdiagonal in the last N+2-MB positions of the second column, 
!    further subdiagonals similarly, and finally its principal diagonal in 
!    the N positions of the last column.  Contents of storages not part of 
!    the matrix are arbitrary.  On output, A has been destroyed, except for 
!    its last two columns which contain a copy of the tridiagonal matrix.
!
!    Output, real(kind=8) D(N), the diagonal elements of the tridiagonal matrix.
!
!    Output, real(kind=8) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is set to zero.
!
!    Output, real(kind=8) E2(N), contains the squares of the corresponding elements 
!    of E.  E2 may coincide with E if the squares are not needed.
!
!    Input, logical MATZ, should be set to TRUE if the transformation matrix is
!    to be accumulated, and to FALSE otherwise.
!
!    Output, real(kind=8) Z(NM,N), the orthogonal transformation matrix produced in
!    the reduction if MATZ has been set to TRUE.  Otherwise, Z is not
!    referenced.
!
  implicit none
!
  integer mb
  integer n
  integer nm
!
  real(kind=8) a(nm,mb)
  real(kind=8) b1
  real(kind=8) b2
  real(kind=8) c2
  real(kind=8) d(n)
  real(kind=8) dmin
  real(kind=8) dminrt
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) f1
  real(kind=8) f2
  real(kind=8) g
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer k
  integer kr
  integer l
  integer m1
  logical matz
  integer maxl
  integer maxr
  integer mr
  integer r
  integer r1
  real(kind=8) s2
  real(kind=8) u
  integer ugl
  real(kind=8) z(nm,n)
!
  dmin = real(2.0E+00,kind=8)**(-32)
  dminrt = real(2.0E+00,kind=8)**(-16)
!
!  Initialize the diagonal scaling matrix.
!
  d(1:n) = real(1.0E+00,kind=8)

  if ( matz ) then
    call rmat_ident ( nm, n, z )
  end if

  m1 = mb - 1

  if ( m1 < 1 ) then
    d(1:n) = a(1:n,mb)
    e(1:n) = real(0.0E+00,kind=8)
    e2(1:n) = real(0.0E+00,kind=8)
    return
  end if

  if ( m1 == 1 ) then
    go to 800
  end if

  do k = 1, n - 2

    maxr = min ( m1, n-k )

    do r1 = 2, maxr

      r = maxr + 2 - r1
      kr = k + r
      mr = mb - r
      g = a(kr,mr)
      a(kr-1,1) = a(kr-1,mr+1)
      ugl = k

      do j = kr, n, m1

        j1 = j - 1
        j2 = j1 - 1

        if ( g == real(0.0E+00,kind=8) ) go to 600

        b1 = a(j1,1) / g
        b2 = b1 * d(j1) / d(j)
        s2 = real(1.0E+00,kind=8) / ( real(1.0E+00,kind=8) + b1 * b2 )

        if ( s2 < real(0.5E+00,kind=8) ) then

          b1 = g / a(j1,1)
          b2 = b1 * d(j) / d(j1)
          c2 = real(1.0E+00,kind=8) - s2
          d(j1) = c2 * d(j1)
          d(j) = c2 * d(j)
          f1 = real(2.0E+00,kind=8) * a(j,m1)
          f2 = b1 * a(j1,mb)
          a(j,m1) = -b2 * ( b1 * a(j,m1) - a(j,mb) ) - f2 + a(j,m1)
          a(j1,mb) = b2 * ( b2 * a(j,mb) + f1 ) + a(j1,mb)
          a(j,mb) = b1 * ( f2 - f1 ) + a(j,mb)

          do l = ugl, j2
            i2 = mb - j + l
            u = a(j1,i2+1) + b2 * a(j,i2)
            a(j,i2) = -b1 * a(j1,i2+1) + a(j,i2)
            a(j1,i2+1) = u
          end do

          ugl = j
          a(j1,1) = a(j1,1) + b2 * g

          if ( j /= n ) then

            maxl = min ( m1, n-j1 )

            do l = 2, maxl
              i1 = j1 + l
              i2 = mb - l
              u = a(i1,i2) + b2 * a(i1,i2+1)
              a(i1,i2+1) = -b1 * a(i1,i2) + a(i1,i2+1)
              a(i1,i2) = u
            end do

            i1 = j + m1

            if ( i1 <= n ) then
              g = b2 * a(i1,1)
            end if

          end if

          if ( matz ) then

            do l = 1, n
              u = z(l,j1) + b2 * z(l,j)
              z(l,j) = -b1 * z(l,j1) + z(l,j)
              z(l,j1) = u
            end do

          end if

        else

          u = d(j1)
          d(j1) = s2 * d(j)
          d(j) = s2 * u
          f1 = real(2.0E+00,kind=8) * a(j,m1)
          f2 = b1 * a(j,mb)
          u = b1 * ( f2 - f1 ) + a(j1,mb)
          a(j,m1) = b2 * ( b1 * a(j,m1) - a(j1,mb) ) + f2 - a(j,m1)
          a(j1,mb) = b2 * ( b2 * a(j1,mb) + f1 ) + a(j,mb)
          a(j,mb) = u

          do l = ugl, j2
            i2 = mb - j + l
            u = b2 * a(j1,i2+1) + a(j,i2)
            a(j,i2) = -a(j1,i2+1) + b1 * a(j,i2)
            a(j1,i2+1) = u
          end do

          ugl = j
          a(j1,1) = b2 * a(j1,1) + g

          if ( j /= n ) then

            maxl = min ( m1, n-j1 )

            do l = 2, maxl
              i1 = j1 + l
              i2 = mb - l
              u = b2 * a(i1,i2) + a(i1,i2+1)
              a(i1,i2+1) = -a(i1,i2) + b1 * a(i1,i2+1)
              a(i1,i2) = u
            end do

            i1 = j + m1

            if ( i1 <= n ) then
              g = a(i1,1)
              a(i1,1) = b1 * a(i1,1)
            end if

          end if

          if ( matz ) then

            do l = 1, n
              u = b2 * z(l,j1) + z(l,j)
              z(l,j) = -z(l,j1) + b1 * z(l,j)
              z(l,j1) = u
            end do

          end if

        end if

      end do

600   continue

    end do
!
!  Rescale to avoid underflow or overflow.
!
    if ( mod ( k, 64 ) == 0 ) then

      do j = k, n

        if ( d(j) < dmin ) then

          maxl = max ( 1, mb+1-j )

          a(j,maxl:m1) = dminrt * a(j,maxl:m1)

          if ( j /= n ) then

            maxl = min ( m1, n-j )

            do l = 1, maxl
              i1 = j + l
              i2 = mb - l
              a(i1,i2) = dminrt * a(i1,i2)
            end do

          end if

          if ( matz ) then
            z(1:n,j) = dminrt * z(1:n,j)
          end if

          a(j,mb) = dmin * a(j,mb)
          d(j) = d(j) / dmin
 
        end if

      end do

    end if

  end do
!
!   Form square root of scaling matrix.
!
800 continue

  e(2:n) = sqrt ( d(2:n) )

  if ( matz ) then

    do k = 2, n
      z(1:n,k) = z(1:n,k) * e(k)
    end do

  end if

  u = real(1.0E+00,kind=8)

  do j = 2, n
    a(j,m1) = u * e(j) * a(j,m1)
    u = e(j)
    e2(j) = a(j,m1)**2
    a(j,mb) = d(j) * a(j,mb)
    d(j) = a(j,mb)
    e(j) = a(j,m1)
  end do

  d(1) = a(1,mb)
  e(1) = real(0.0E+00,kind=8)
  e2(1) = real(0.0E+00,kind=8)

  return
end
subroutine bandv ( nm, n, mbw, a, e21, m, w, z, ierr )
!
!*******************************************************************************
!
!! BANDV finds eigenvectors from eigenvalues, for a real(kind=8) symmetric band matrix.
!
!
!  Discussion:
!
!    This subroutine finds those eigenvectors of a real(kind=8) symmetric
!    band matrix corresponding to specified eigenvalues, using inverse
!    iteration.  The subroutine may also be used to solve systems
!    of linear equations with a symmetric or non-symmetric band
!    coefficient matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of A and Z, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MBW, the number of columns of the array A used to store the
!    band matrix.  If the matrix is symmetric, MBW is its (half)
!    band width, denoted MB and defined as the number of adjacent
!    diagonals, including the principal diagonal, required to
!    specify the non-zero portion of the lower triangle of the
!    matrix.  If the subroutine is being used to solve systems
!    of linear equations and the coefficient matrix is not
!    symmetric, it must however have the same number of adjacent
!    diagonals above the main diagonal as below, and in this
!    case, MBW=2*MB-1.
!
!    Input, real(kind=8) A(NM,MBW), the lower triangle of the symmetric band input
!    atrix stored as an N by MB array.  Its lowest subdiagonal is stored 
!    in the last N+1-MB positions of the first column, its next subdiagonal 
!    in the last N+2-MB positions of the second column, further subdiagonals 
!    similarly, and finally its principal diagonal in the N positions of 
!    column MB.  If the subroutine is being used to solve systems of linear
!    equations, and the coefficient matrix is not symmetric, A is
!    N by 2*MB-1 instead, with lower triangle as above and with its first 
!    superdiagonal stored in the first N-1 positions of column MB+1, its 
!    second superdiagonal in the first N-2 positions of column MB+2, further
!    superdiagonals similarly, and finally its highest superdiagonal in 
!    the first N+1-MB positions of the last column.  Contents of storages 
!    not part of the matrix are arbitrary.
!
!    Input, real(kind=8) E21, specifies the ordering of the eigenvalues and contains
!    real(0.0E+00,kind=8) if the eigenvalues are in ascending order, or real(2.0E+00,kind=8) if the 
!    eigenvalues are in descending order.  If the subroutine is being used 
!    to solve systems of linear equations, E21 should be set to real(1.0E+00,kind=8) 
!    if the coefficient matrix is symmetric and to -real(1.0E+00,kind=8) if not.
!
!    Input, integer M, the number of specified eigenvalues or the number of
!    systems of linear equations.
!
!    Input, real(kind=8) W(M), contains the M eigenvalues in ascending or descending 
!    order.  If the subroutine is being used to solve systems of linear
!    equations (A-W(1:M)*I) * X(1:M) = B(1:M), where I is the identity matrix, 
!    W should be set accordingly.
!
!    Input/output, real(kind=8) Z(NM,M).  On input, the constant matrix 
!    columns B(1:M), if the subroutine is used to solve systems of linear 
!    equations.  On output, the associated set of orthogonal eigenvectors.
!    Any vector which fails to converge is set to zero.  If the
!    routine is used to solve systems of linear equations,
!    Z contains the solution matrix columns X(1:M).
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    -R, if the eigenvector corresponding to the R-th eigenvalue fails to
!    converge, or if the R-th system of linear equations is nearly singular.
!
  implicit none
!
  integer mbw
  integer n
  integer nm
!
  real(kind=8) a(nm,mbw)
  real(kind=8) e21
  real(kind=8) eps2
  real(kind=8) eps3
  real(kind=8) eps4
  integer group
  integer i
  integer ierr
  integer ii
  integer ij
  integer ij1
  integer its
  integer j
  integer jj
  integer k
  integer kj
  integer kj1
  integer m
  integer m1
  integer m21
  integer maxj
  integer maxk
  integer mb
  real(kind=8) norm
  real(kind=8) order
  real(kind=8) pythag
  integer r
  real(kind=8) rv(n*(2*mbw-1))
  real(kind=8) rv6(n)
  real(kind=8) u
  real(kind=8) uk
  real(kind=8) v
  real(kind=8) w(m)
  real(kind=8) x0
  real(kind=8) x1
  real(kind=8) xu
  real(kind=8) z(nm,m)
!
  ierr = 0

  if ( m == 0 ) then
    return
  end if

  x0 = real(0.0E+00,kind=8)

  if ( e21 < real(0.0E+00,kind=8) ) then
    mb = ( mbw + 1 ) / 2
  else
    mb = mbw
  end if

  m1 = mb - 1
  m21 = m1 + mb
  order = real(1.0E+00,kind=8) - abs ( e21 )
!
!  Find vectors by inverse iteration.
!
  do r = 1, m

     its = 1
     x1 = w(r)
     if ( r /= 1 ) go to 100
!
!  Compute norm of matrix.
!
     norm = real(0.0E+00,kind=8)

     do j = 1, mb

        jj = mb + 1 - j
        kj = jj + m1
        ij = 1
        v = real(0.0E+00,kind=8)

        do i = jj, n

          v = v + abs ( a(i,j) )

          if ( e21 < real(0.0E+00,kind=8) ) then
            v = v + abs ( a(ij,kj) )
            ij = ij + 1
          end if

        end do

        norm = max ( norm, v )

     end do

     if ( e21 < real(0.0E+00,kind=8) ) then
       norm = real(0.5E+00,kind=8) * norm
     end if
!
!  EPS2 is the criterion for grouping,
!  EPS3 replaces zero pivots and equal roots are modified by eps3,
!  EPS4 is taken very small to avoid overflow.
!
     if ( norm == real(0.0E+00,kind=8) ) then
       norm = real(1.0E+00,kind=8)
     end if

     eps2 = 1.0E-03 * norm * abs ( order)
     eps3 = abs ( norm ) * epsilon ( norm )
     uk = n
     uk = sqrt ( uk )
     eps4 = uk * eps3

80   continue

     group = 0
     go to 120
!
!  Look for close or coincident roots.
!
100  continue

     if ( abs ( x1 - x0 ) >= eps2 ) go to 80
     group = group + 1

     if ( order * ( x1 - x0 ) <= real(0.0E+00,kind=8) ) then
       x1 = x0 + order * eps3
     end if
!
!  Expand matrix, subtract eigenvalue, and initialize vector.
!
120  continue

     do i = 1, n

        ij = i + min ( 0, i-m1 ) * n
        kj = ij + mb * n
        ij1 = kj + m1 * n

        if ( m1 == 0 ) go to 180

        do j = 1, m1

          if ( ij <= m1 ) then
            if ( ij <= 0 ) then
              rv(ij1) = real(0.0E+00,kind=8)
              ij1 = ij1 + n
            end if
          else
            rv(ij) = a(i,j)
          end if

          ij = ij + n
          ii = i + j

          if ( ii <= n ) then

            jj = mb - j

            if ( e21 < real(0.0E+00,kind=8) ) then
              ii = i
              jj = mb + j
            end if

            rv(kj) = a(ii,jj)
            kj = kj + n

          end if

        end do

  180   continue

        rv(ij) = a(i,mb) - x1
        rv6(i) = eps4
        if ( order == real(0.0E+00,kind=8) ) then
          rv6(i) = z(i,r)
        end if

     end do

     if ( m1 /= 0 ) then
!
!  Elimination with interchanges.
!
     do i = 1, n

        ii = i + 1
        maxk = min ( i+m1-1, n )
        maxj = min ( n-i, m21-2 ) * n

        do k = i, maxk

           kj1 = k
           j = kj1 + n
           jj = j + maxj

           do kj = j, jj, n
             rv(kj1) = rv(kj)
             kj1 = kj
           end do

           rv(kj1) = real(0.0E+00,kind=8)

        end do

        if ( i /= n ) then

        u = real(0.0E+00,kind=8)
        maxk = min ( i+m1, n )
        maxj = min ( n-ii, m21-2 ) * n

        do j = i, maxk
          if ( abs ( rv(j) ) >= abs ( u ) ) then
            u = rv(j)
            k = j
          end if
        end do

        j = i + n
        jj = j + maxj

        if ( k /= i ) then

          kj = k

          do ij = i, jj, n
            call r_swap ( rv(ij), rv(kj) )
            kj = kj + n
          end do

          if ( order == real(0.0E+00,kind=8) ) then
            call r_swap ( rv6(i), rv6(k) )
          end if

        end if

        if ( u /= real(0.0E+00,kind=8) ) then

        do k = ii, maxk

           v = rv(k) / u
           kj = k

           do ij = j, jj, n
             kj = kj + n
             rv(kj) = rv(kj) - v * rv(ij)
           end do

           if ( order == real(0.0E+00,kind=8) ) then
             rv6(k) = rv6(k) - v * rv6(i)
           end if

        end do

       end if

      end if

      end do

     end if
!
!  Back substitution.
!
600  continue

     do ii = 1, n

        i = n + 1 - ii
        maxj = min ( ii, m21 )

        if ( maxj /= 1 ) then

          ij1 = i
          j = ij1 + n
          jj = j + (maxj - 2) * n

          do ij = j, jj, n
            ij1 = ij1 + 1
            rv6(i) = rv6(i) - rv(ij) * rv6(ij1)
          end do

        end if

        v = rv(i)
!
!  Error: nearly singular linear system.
!
        if ( abs ( v ) < eps3) then
          if ( order == real(0.0E+00,kind=8) ) then
            ierr = -r
          end if
          v = sign ( eps3, v )
        end if

        rv6(i) = rv6(i) / v

     end do

     xu = real(1.0E+00,kind=8)

     if ( order == real(0.0E+00,kind=8) ) go to 870
!
!  Orthogonalize with respect to previous members of group.
!
     do jj = 1, group

        j = r - group - 1 + jj

        xu = dot_product ( rv6(1:n), z(1:n,j) )

        rv6(1:n) = rv6(1:n) - xu * z(1:n,j)

     end do

     norm = sum ( abs ( rv6(1:n) ) )
!
!  Choose a new starting vector.
!
     if ( norm < 0.1E+00 ) then

       if ( its < n ) then
         its = its + 1
         xu = eps4 / ( uk + real(1.0E+00,kind=8) )
         rv6(1) = eps4
         rv6(2:n) = xu
         rv6(its) = rv6(its) - eps4 * uk
         go to 600
       else
         ierr = -r
         xu = real(0.0E+00,kind=8)
         go to 870
       end if

     end if
!
!   Normalize so that sum of squares is 1 and expand to full order.
!
     u = real(0.0E+00,kind=8)
     do i = 1, n
       u = pythag ( u, rv6(i) )
     end do

     xu = real(1.0E+00,kind=8) / u

 870 continue

     z(1:n,r) = rv6(1:n) * xu

     x0 = x1

  end do

  return
end
subroutine bisect ( n, eps1, d, e, e2, lb, ub, mm, m, w, ind, ierr )
!
!*******************************************************************************
!
!! BISECT computes some eigenvalues of a real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds those eigenvalues of a real(kind=8) symmetric
!    tridiagonal matrix which lie in a specified interval, using bisection.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2, pages 418-439,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) EPS1, is an absolute error tolerance for the computed
!    eigenvalues.  If the input EPS1 is non-positive, it is reset for each 
!    submatrix to a default value, namely, minus the product of the relative 
!    machine precision and the 1-norm of the submatrix.
!
!    Input, real(kind=8) D(N), the diagonal elements of the input matrix.
!
!    Input, real(kind=8) E(N), contains in E(2:N) the subdiagonal elements of the
!    matrix.  E(1) is arbitrary.
!
!    Input/output, real(kind=8) E2(N).  On input, the squares of the corresponding 
!    elements of E.  E2(1) is arbitrary.  On output, elements of E2, 
!    corresponding to elements of E regarded as negligible, have been 
!    replaced by zero, causing the matrix to split into a direct sum of 
!    submatrices.  E2(1) is also set to zero.
!
!    Input, real(kind=8) LB, UB, define the interval to be searched for eigenvalues.
!    If LB is not less than UB, no eigenvalues will be found.
!
!    Input, integer MM, an upper bound for the number of eigenvalues in the 
!    interval.  Warning: if more than MM eigenvalues are determined to lie 
!    in the interval, an error return is made with no eigenvalues found.
!
!    Output, integer M, the number of eigenvalues determined to lie 
!    in (LB,UB).
!
!    Output, real(kind=8) W(M), the eigenvalues in ascending order.
!
!    Output, integer IND(MM), contains in its first M positions the submatrix 
!    indices associated with the corresponding eigenvalues in W:
!    1 for eigenvalues belonging to the first submatrix from the top, 2 for 
!    those belonging to the second submatrix, and so on.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    3*N+1, if M exceeds MM.
!
  implicit none
!
  integer mm
  integer n
!
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) eps1
  integer i
  integer ierr
  integer ii
  integer ind(mm)
  integer isturm
  integer j
  integer k
  integer l
  real(kind=8) lb
  integer m
  integer m1
  integer m2
  integer p
  integer q
  integer r
  real(kind=8) rv4(n)
  real(kind=8) rv5(n)
  integer s
  real(kind=8) t1
  real(kind=8) t2
  integer tag
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) u
  real(kind=8) ub
  real(kind=8) v
  real(kind=8) w(mm)
  real(kind=8) x0
  real(kind=8) x1
  real(kind=8) xu
!
  ierr = 0
  s = 0
  tag = 0
  t1 = lb
  t2 = ub
!
!  Look for small sub-diagonal entries.
!
  e2(1) = real(0.0E+00,kind=8)

  do i = 2, n
    
    tst1 = abs ( d(i) ) + abs ( d(i-1) )
    tst2 = tst1 + abs ( e(i) )

    if ( tst2 <= tst1 ) then
      e2(i) = real(0.0E+00,kind=8)
    end if

  end do
!
!  Determine the number of eigenvalues in the interval.
!
  p = 1
  q = n
  x1 = ub
  isturm = 1
  go to 320

60 continue

  m = s
  x1 = lb
  isturm = 2
  go to 320

80 continue

  m = m - s

  if ( m > mm ) go to 980

  q = 0
  r = 0
!
!  Establish and process next submatrix, refining
!  interval by the Gerschgorin bounds.
!
100 continue

  if ( r == m ) go to 1001

  tag = tag + 1
  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = real(0.0E+00,kind=8)

  do q = p, n

    x1 = u
    u = real(0.0E+00,kind=8)
    v = real(0.0E+00,kind=8)

    if ( q /= n ) then
      u = abs ( e(q+1) )
      v = e2(q+1)
    end if

    xu = min ( d(q) - ( x1 + u ), xu )
    x0 = max ( d(q) + ( x1 + u ), x0 )

    if ( v == real(0.0E+00,kind=8) ) then
      exit
    end if

  end do

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )
  if ( eps1 <= real(0.0E+00,kind=8) ) then
    eps1 = -x1
  end if

  if ( p /= q ) go to 180
!
!  Check for an isolated root within interval.
!
  if ( t1 > d(p) .or. d(p) >= t2 ) go to 940
  m1 = p
  m2 = p
  rv5(p) = d(p)
  go to 900

  180 continue

  x1 = x1 * (q - p + 1)
  lb = max ( t1, xu-x1 )
  ub = min ( t2, x0+x1 )
  x1 = lb
  isturm = 3
  go to 320

  200 continue

  m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320

  220 continue

  m2 = s
  if ( m1 > m2 ) go to 940
!
!  Find roots by bisection.
!
  x0 = ub
  isturm = 5
  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for the K-th eigenvalue.
!
  k = m2

250 continue

     xu = lb

     do ii = m1, k
       i = m1 + k - ii
       if ( xu < rv4(i) ) then
         xu = rv4(i)
         go to 280
       end if
     end do

  280 continue

   x0 = min ( x0, rv5(k) )
!
!  Next bisection step.
!
  300    continue

     x1 = ( xu + x0 ) * real(0.5E+00,kind=8)

     if ( (x0 - xu) <= abs ( eps1 ) ) go to 420

     tst1 = real(2.0E+00,kind=8) * ( abs ( xu ) + abs ( x0 ) )
     tst2 = tst1 + ( x0 - xu )
     if ( tst2 == tst1 ) go to 420
!
!  Sturm sequence.
!
320  continue

     s = p - 1
     u = real(1.0E+00,kind=8)

     do i = p, q

        if ( u == real(0.0E+00,kind=8) ) then
          v = abs ( e(i) ) / epsilon ( v )
          if ( e2(i) == real(0.0E+00,kind=8) ) v = real(0.0E+00,kind=8)
        else
          v = e2(i) / u
        end if

        u = d(i) - x1 - v
        if ( u < real(0.0E+00,kind=8) ) then
          s = s + 1
        end if

     end do

     go to (60,80,200,220,360), isturm
!
!  Refine intervals.
!
  360 continue

     if ( s >= k ) go to 400
     xu = x1

     if ( s < m1 ) then
       rv4(m1) = x1
       go to 300
     end if

  380 continue
     rv4(s+1) = x1
     if ( rv5(s) > x1 ) rv5(s) = x1
     go to 300
400  continue
     x0 = x1
     go to 300
!
!  K-th eigenvalue found.
!
420 continue

  rv5(k) = x1
  k = k - 1
  if ( k >= m1 ) go to 250
!
!  Order eigenvalues tagged with their submatrix associations.
!
900 continue

  s = r
  r = r + m2 - m1 + 1
  j = 1
  k = m1

  do l = 1, r

    if ( j <= s ) then

      if ( k > m2 ) then
        exit
      end if

      if ( rv5(k) >= w(l) ) then
        j = j + 1
        cycle
      end if

      do ii = j, s
        i = l + s - ii
        w(i+1) = w(i)
        ind(i+1) = ind(i)
      end do

    end if

    w(l) = rv5(k)
    ind(l) = tag
    k = k + 1

  end do

940 continue

  if ( q < n ) go to 100

  go to 1001
!
!  Set error: underestimate of number of eigenvalues in interval.
!
980 continue

  ierr = 3 * n + 1

 1001 continue

  lb = t1
  ub = t2

  return
end
subroutine bqr ( nm, n, mb, a, t, r, ierr )
!
!*******************************************************************************
!
!! BQR finds the smallest eigenvalue of a real(kind=8) symmetric band matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalue of smallest magnitude of a real(kind=8) 
!    symmetric band matrix using the QR algorithm with shifts of origin.  
!    Consecutive calls can be made to find further eigenvalues.
!
!    Note that for a subsequent call, N should be replaced by N-1, but
!    MB should not be altered even when it exceeds the current N.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MB, the (half) band width of the matrix, defined as the
!    number of adjacent diagonals, including the principal
!    diagonal, required to specify the non-zero portion of the
!    lower triangle of the matrix.
!
!    Input/output, real(kind=8) A(NM,MB).  On input, A contains the lower triangle 
!    of the symmetric band input matrix stored as an N by MB array.  Its 
!    lowest subdiagonal is stored in the last N+1-MB positions of the first 
!    column, its next subdiagonal in the last N+2-MB positions of the
!    second column, further subdiagonals similarly, and finally its principal 
!    diagonal in the N positions of the last column.  Contents of storages 
!    not part of the matrix are arbitrary.  On a subsequent call, its output
!    contents from the previous call should be passed.  On output, A contains 
!    the transformed band matrix.  The matrix A+T*I derived from the output 
!    parameters is similar to the input A+T*I to within rounding errors.  
!    Its last row and column are null as long as IERR is zero.
!
!    Input/output, real(kind=8) T.  On input, T specifies the shift (of eigenvalues) 
!    applied to the diagonal of A in forming the input matrix.  What is 
!    actually determined is the eigenvalue nearest to T of A+T*I, where I 
!    is the identity matrix.  On a subsequent call, the output value of T 
!    from the previous call should be passed if the next nearest eigenvalue
!    is sought.  On output, T contains the computed eigenvalue of A+T*I,
!    as long as IERR is zero.
!
!    Input/output, real(kind=8) R.  On input for the first call, R should be 
!    specified as zero, and as its output value from the previous call 
!    on a subsequent call.  It is used to determine when the last row and 
!    column of the transformed band matrix can be regarded as negligible.
!    On output, R contains the maximum of its input value and the norm of the
!    last column of the input matrix A.
!
!    Output, integer IERR, error flag.
!    0, normal return.
!    N, if the eigenvalue has not been determined after 30 iterations.
!
  implicit none
!
  integer mb
  integer nm
!
  real(kind=8) a(nm,mb)
  real(kind=8) f
  real(kind=8) g
  integer i
  integer ierr
  integer ii
  integer ik
  integer imult
  integer its
  integer j
  integer jk
  integer jm
  integer k
  integer kj
  integer kj1
  integer kk
  integer km
  integer l
  integer ll
  integer m
  integer m1
  integer m2
  integer m21
  integer m3
  integer m31
  integer m4
  integer mk
  integer mn
  integer mz
  integer n
  integer ni
  real(kind=8) pythag
  real(kind=8) q
  real(kind=8) r
  real(kind=8) rv(2*mb*mb+4*mb-3)
  real(kind=8) s
  real(kind=8) scale
  real(kind=8) t
  real(kind=8) tst1
  real(kind=8) tst2
!
  ierr = 0
  m1 = min ( mb, n )
  m = m1 - 1
  m2 = m + m
  m21 = m2 + 1
  m3 = m21 + m
  m31 = m3 + 1
  m4 = m31 + m2
  mn = m + n
  mz = mb - m1
  its = 0
!
!  Test for convergence.
!
40 continue

  g = a(n,mb)

  if ( m == 0 ) go to 360

  f = real(0.0E+00,kind=8)
  do k = 1, m
    mk = k + mz
    f = f + abs ( a(n,mk) )
  end do

  if ( its == 0 .and. f > r ) then
    r = f
  end if

  tst1 = r
  tst2 = tst1 + f

  if ( tst2 <= tst1 ) go to 360

  if ( its >= 30 ) then
    ierr = n
    return
  end if

  its = its + 1
!
!  Form shift from bottom 2 by 2 minor.
!
  if ( f <= 0.25E+00 * r .or. its >= 5 ) then

    f = a(n,mb-1)

    if ( f /= real(0.0E+00,kind=8) ) then
      q = ( a(n-1,mb) - g ) / ( real(2.0E+00,kind=8) * f )
      s = pythag ( q, real(1.0E+00,kind=8) )
      g = g - f / ( q + sign ( s, q ) )
    end if

    t = t + g

    a(1:n,mb) = a(1:n,mb) - g

  end if

  rv(m31:m4) = real(0.0E+00,kind=8)

  do ii = 1, mn

     i = ii - m
     ni = n - ii

     if ( ni < 0 ) go to 230
!
!  Form column of shifted matrix A-G*I.
!
     l = max ( 1, 2-i )

     rv(1:m3) = real(0.0E+00,kind=8)

     do k = l, m1
       km = k + m
       mk = k + mz
       rv(km) = a(ii,mk)
     end do

     ll = min ( m, ni )

     do k = 1, ll
       km = k + m21
       ik = ii + k
       mk = mb - k
       rv(km) = a(ik,mk)
     end do
!
!  Pre-multiply with Householder reflections.
!
     ll = m2
     imult = 0
!
!  Multiplication procedure.
!
140  continue

     kj = m4 - m1

     do j = 1, ll

        kj = kj + m1
        jm = j + m3

        if ( rv(jm) /= real(0.0E+00,kind=8) ) then

          f = real(0.0E+00,kind=8)

          do k = 1, m1
            kj = kj + 1
            jk = j + k - 1
            f = f + rv(kj) * rv(jk)
          end do

          f = f / rv(jm)
          kj = kj - m1

          do k = 1, m1
            kj = kj + 1
            jk = j + k - 1
            rv(jk) = rv(jk) - rv(kj) * f
          end do

          kj = kj - m1

        end if

     end do

     if ( imult /= 0 ) go to 280
!
!  Householder reflection.
!
     f = rv(m21)
     s = real(0.0E+00,kind=8)
     rv(m4) = real(0.0E+00,kind=8)
     scale = sum ( abs ( rv(m21:m3) ) )

     if ( scale == real(0.0E+00,kind=8) ) go to 210

     do k = m21, m3
       s = s + ( rv(k) / scale )**2
     end do

     s = scale * scale * s
     g = - sign ( sqrt ( s ), f )
     rv(m21) = g
     rv(m4) = s - f * g
     kj = m4 + m2 * m1 + 1
     rv(kj) = f - g

     do k = 2, m1
       kj = kj + 1
       km = k + m2
       rv(kj) = rv(km)
     end do
!
!  Save column of triangular factor R.
!
210  continue

     do k = l, m1
       km = k + m
       mk = k + mz
       a(ii,mk) = rv(km)
     end do

230  continue

     l = max ( 1, m1+1-i )
     if ( i <= 0 ) go to 300
!
!  Perform additional steps.
!
     rv(1:m21) = real(0.0E+00,kind=8)
     ll = min ( m1, ni+m1 )
!
!  Get row of triangular factor R.
!
     do kk = 1, ll
       k = kk - 1
       km = k + m1
       ik = i + k
       mk = mb - k
       rv(km) = a(ik,mk)
     end do
!
!  Post-multiply with Householder reflections.
!
     ll = m1
     imult = 1
     go to 140
!
!  Store column of new a matrix.
!
280  continue

     do k = l, m1
       mk = k + mz
       a(i,mk) = rv(k)
     end do
!
!  Update Householder reflections.
!
300  continue

     if ( l > 1 ) l = l - 1

     kj1 = m4 + l * m1

     do j = l, m2

       jm = j + m3
       rv(jm) = rv(jm+1)

       do k = 1, m1
         kj1 = kj1 + 1
         kj = kj1 - m1
         rv(kj) = rv(kj1)
       end do

     end do

  end do

  go to 40
!
!  Convergence.
!
360 continue

  t = t + g
  a(1:n,mb) = a(1:n,mb) - g

  do k = 1, m1
    mk = k + mz
    a(n,mk) = real(0.0E+00,kind=8)
  end do

  return
end
subroutine cbabk2 ( nm, n, low, igh, scale, m, zr, zi )
!
!*******************************************************************************
!
!! CBABK2 finds eigenvectors by undoing the CBAL transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a complex general
!    matrix by back transforming those of the corresponding
!    balanced matrix determined by CBAL.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of ZI and ZR, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, values determined by CBAL.
!
!    Input, real(kind=8) SCALE(N), information determining the permutations
!    and scaling factors used by CBAL.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) ZR(NM,M), ZI(NM,M).  On input, the real(kind=8) and imaginary 
!    parts, respectively, of the eigenvectors to be back transformed in 
!    their first M columns.  On output, the transformed eigenvectors.
!
  implicit none
!
  integer m
  integer n
  integer nm
!
  integer i
  integer igh
  integer ii
  integer j
  integer k
  integer low
  real(kind=8) s
  real(kind=8) scale(n)
  real(kind=8) zi(nm,m)
  real(kind=8) zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  if ( igh /= low ) then

    do i = low, igh

      s = scale(i)

      zr(i,1:m) = zr(i,1:m) * s
      zi(i,1:m) = zi(i,1:m) * s

    end do

  end if

  do ii = 1, n

    i = ii

    if ( i < low .or. i > igh ) then

      if ( i < low ) then
        i = low - ii
      end if

      k = scale(i)

      if ( k /= i ) then

        do j = 1, m
          call r_swap ( zr(i,j), zr(k,j) )
          call r_swap ( zi(i,j), zi(k,j) )
        end do

      end if

    end if

  end do

  return
end
subroutine cbal ( nm, n, ar, ai, low, igh, scale )
!
!*******************************************************************************
!
!! CBAL balances a complex matrix before eigenvalue calculations.
!
!
!  Discussion:
!
!    This subroutine balances a complex matrix and isolates
!    eigenvalues whenever possible.
!
!    Suppose that the principal submatrix in rows low through igh
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!      SCALE(J) = P(J),    for J = 1,...,LOW-1
!               = D(J,J)       J = LOW,...,IGH
!               = P(J)         J = IGH+1,...,N.
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is returned for IGH if IGH is zero formally.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of AR and AI, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) AR(NM,N), AI(NM,N).  On input, the real(kind=8) and 
!    imaginary parts of the complex matrix to be balanced.  On output,
!    the real(kind=8) and imaginary parts of the balanced matrix.
!
!    Output, integer LOW, IGH, are values such that AR(I,J) and AI(I,J)
!    are zero if I is greater than J and either J=1,...,LOW-1 or 
!    I=IGH+1,...,N.
!
!    Output, real(kind=8) SCALE(N), information determining the
!    permutations and scaling factors used.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  real(kind=8) b2
  real(kind=8) c
  real(kind=8) f
  real(kind=8) g
  integer i
  integer iexc
  integer igh
  integer j
  integer jj
  integer k
  integer l
  integer low
  integer m
  logical noconv
  real(kind=8) r
  real(kind=8) radix
  real(kind=8) s
  real(kind=8) scale(n)
!
  radix = 16.0E+00

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

20 continue

  scale(m) = j

  if ( j /= m ) then

    do i = 1, l
      call r_swap ( ar(i,j), ar(i,m) )
      call r_swap ( ai(i,j), ai(i,m) )
    end do

    do i = k, n
      call r_swap ( ar(j,i), ar(m,i) )
      call r_swap ( ai(j,i), ai(m,i) )
    end do

  end if

  if ( iexc == 2 ) then
    go to 130
  end if
!
!  Search for rows isolating an eigenvalue and push them down.
!
80 continue

  if ( l == 1 ) go to 280
  l = l - 1

100 continue

  do jj = 1, l

     j = l + 1 - jj

     do i = 1, l
       if ( i /= j ) then
         if ( ar(j,i) /= real(0.0E+00,kind=8) .or. ai(j,i) /= real(0.0E+00,kind=8) ) go to 120
       end if
     end do

     m = l
     iexc = 1
     go to 20

120  continue

  end do

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 continue

  k = k + 1

140 continue

   do j = k, l

     do i = k, l
       if ( i /= j ) then
         if ( ar(i,j) /= real(0.0E+00,kind=8) .or. ai(i,j) /= real(0.0E+00,kind=8) ) go to 170
       end if
     end do

     m = k
     iexc = 2
     go to 20
170  continue

  end do
!
!  Now balance the submatrix in rows k to l.
!
  scale(k:l) = real(1.0E+00,kind=8)
!
!  Iterative loop for norm reduction.
!
190 continue

  noconv = .false.

  do i = k, l

    c = real(0.0E+00,kind=8)
    r = real(0.0E+00,kind=8)

    do j = k, l
      if ( j /= i ) then
        c = c + abs ( ar(j,i) ) + abs ( ai(j,i) )
        r = r + abs ( ar(i,j) ) + abs ( ai(i,j) )
      end if
    end do
!
!  Guard against zero C or R due to underflow.
!
     if ( c == real(0.0E+00,kind=8) .or. r == real(0.0E+00,kind=8) ) go to 270

     g = r / radix
     f = real(1.0E+00,kind=8)
     s = c + r

     do while ( c < g )
       f = f * radix
       c = c * b2
     end do

     g = r * radix

     do while  ( c >= g )
       f = f / radix
       c = c / b2
     end do
!
!  Now balance.
!
     if ( ( c + r ) / f < 0.95E+00 * s ) then

       g = real(1.0E+00,kind=8) / f
       scale(i) = scale(i) * f
       noconv = .true.

       ar(i,k:n) = ar(i,k:n) * g
       ai(i,k:n) = ai(i,k:n) * g

       ar(1:l,i) = ar(1:l,i) * f
       ai(1:l,i) = ai(1:l,i) * f

     end if

270  continue

  end do

  if ( noconv ) go to 190

  280 continue

  low = k
  igh = l

  return
end
subroutine cdiv ( ar, ai, br, bi, cr, ci )
!
!*******************************************************************************
!
!! CDIV emulates complex division, using real(kind=8) arithmetic.
!
!
!  Discussion:
!
!    This routine performs complex division:
!
!      (CR,CI) = (AR,AI) / (BR,BI)
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, real(kind=8) AR, AI, the real(kind=8) and imaginary parts of the numerator.
!
!    Input, real(kind=8) BR, BI, the real(kind=8) and imaginary parts of the denominator.
!
!    Output, real(kind=8) CR, CI, the real(kind=8) and imaginary parts of the result.
!
  implicit none
!
  real(kind=8) ai
  real(kind=8) ais
  real(kind=8) ar
  real(kind=8) ars
  real(kind=8) bi
  real(kind=8) bis
  real(kind=8) br
  real(kind=8) brs
  real(kind=8) ci
  real(kind=8) cr
  real(kind=8) s
!
  s = abs ( br ) + abs ( bi )

  ars = ar / s
  ais = ai / s
  brs = br / s
  bis = bi / s

  s = brs**2 + bis**2
  cr = ( ars * brs + ais * bis ) / s
  ci = ( ais * brs - ars * bis ) / s

  return
end
subroutine cg ( nm, n, ar, ai, wr, wi, matz, zr, zi, ierr )
!
!*******************************************************************************
!
!! CG gets eigenvalues and eigenvectors of a complex general matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of EISPACK subroutines 
!    to find the eigenvalues and eigenvectors (if desired)
!    of a complex general matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the arrays, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) AR(NM,N), AI(NM,N).  On input, the real(kind=8) and 
!    imaginary parts of the complex matrix.  On output, AR and AI
!    have been overwritten by other information.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts
!    of the eigenvalues.
!
!    Input, integer MATZ, is 0 if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are to be computed.
!
!    Output, real(kind=8) ZR(NM,N), ZI(NM,N), the real(kind=8) and imaginary parts,
!    respectively, of the eigenvectors, if MATZ is not zero.
!
!    Output, integer IERR, an error completion code described in the 
!    documentation for COMQR and COMQR2.  The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  real(kind=8) fv3(n)
  integer ierr
  integer is1
  integer is2
  integer matz
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) zi(nm,n)
  real(kind=8) zr(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call cbal ( nm, n, ar, ai, is1, is2, fv1 )

  call corth ( nm, n, is1, is2, ar, ai, fv2, fv3 )

  if ( matz == 0 ) then

    call comqr ( nm, n, is1, is2, ar, ai, wr, wi, ierr )

    if ( ierr /= 0 ) then
      return
    end if
  
  else

    call comqr2 ( nm, n, is1, is2, fv2, fv3, ar, ai, wr, wi, zr, zi, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CG - Fatal error!'
      write ( *, '(a)' ) '  Nonzero error return from COMQR2.'
      return
    end if

    call cbabk2 ( nm, n, is1, is2, fv1, n, zr, zi )

  end if

  return
end
subroutine ch ( nm, n, ar, ai, w, matz, zr, zi, ierr )
!
!*******************************************************************************
!
!! CH gets eigenvalues and eigenvectors of a complex Hermitian matrix. .
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of subroutines from the 
!    EISPACK eigensystem package to find the eigenvalues and eigenvectors 
!    of a complex hermitian matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the arrays, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) AR(NM,N), AI(NM,N).  On input, the real(kind=8) and 
!    imaginary parts of the complex matrix.  On output, AR and AI
!    have been overwritten by other information.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Input, integer MATZ, is 0 if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are to be computed.
!
!    Output, real(kind=8) ZR(NM,N), ZI(NM,N), the real(kind=8) and imaginary parts,
!    respectively, of the eigenvectors, if MATZ is not zero.
!
!    Output, integer IERR, an error completion code described in the 
!    documentation for TQLRAT and TQL2.  The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  real(kind=8) fm1(2,n)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  integer ierr
  integer matz
  real(kind=8) w(n)
  real(kind=8) zi(nm,n)
  real(kind=8) zr(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call htridi ( nm, n, ar, ai, w, fv1, fv2, fm1 )

  if ( matz == 0 ) then

    call tqlrat ( n, w, fv2, ierr )

  else

    call rmat_ident ( nm, n, zr )

    call tql2 ( nm, n, w, fv1, zr, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call htribk ( nm, n, ar, ai, fm1, n, zr, zi )

  end if

  return
end
subroutine cinvit ( nm, n, ar, ai, wr, wi, select, mm, m, zr, zi, ierr )
!
!*******************************************************************************
!
!! CINVIT gets eigenvectors from eigenvalues, for a complex Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds those eigenvectors of a complex upper
!    Hessenberg matrix corresponding to specified eigenvalues,
!    using inverse iteration.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the arrays, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) AR(NM,N), AI(NM,N), the real(kind=8) and imaginary parts of 
!    the complex Hessenberg matrix.  
!
!    Input/output, real(kind=8) WR(N), WI(N).  On input, the real(kind=8) and imaginary parts
!    of the eigenvalues of the matrix.  The eigenvalues must be stored in a 
!    manner identical to that of subroutine COMLR, which recognizes possible 
!    splitting of the matrix.  On output, WR may have been altered since 
!    close eigenvalues are perturbed slightly in searching for independent 
!    eigenvectors.
!
!    Input, logical SELECT(N), specifies the eigenvectors to be found.  The
!    eigenvector corresponding to the J-th eigenvalue is specified by 
!    setting SELECT(J) to TRUE.
!
!    Input, integer MM, an upper bound for the number of eigenvectors 
!    to be found.
!
!    Output, integer M, the number of eigenvectors actually found.
!
!    Output, real(kind=8) ZR(NM,MM), ZI(NM,MM), the real(kind=8) and imaginary parts
!    of the eigenvectors.  The eigenvectors are normalized so that the 
!    component of largest magnitude is 1.
!    Any vector which fails the acceptance test is set to zero.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    -(2*N+1), if more than MM eigenvectors have been specified,
!    -K, if the iteration corresponding to the K-th value fails,
!    -(N+K), if both error situations occur.
!
  implicit none
!
  integer mm
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  real(kind=8) eps3
  real(kind=8) growto
  integer i
  integer ierr
  integer ii
  real(kind=8) ilambd
  integer its
  integer j
  integer k
  integer km1
  integer m
  integer mp
  real(kind=8) norm
  real(kind=8) normv
  real(kind=8) pythag
  real(kind=8) rlambd
  real(kind=8) rm1(n,n)
  real(kind=8) rm2(n,n)
  real(kind=8) rv1(n)
  real(kind=8) rv2(n)
  integer s
  logical select(n)
  integer uk
  real(kind=8) ukroot
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) x
  real(kind=8) y
  real(kind=8) zi(nm,mm)
  real(kind=8) zr(nm,mm)
!
  ierr = 0
  uk = 0
  s = 1

  do k = 1, n

    if ( .not. select(k) ) then
      cycle
    end if

    if ( s > mm ) go to 1000

    if ( uk >= k ) go to 200
!
!  Check for possible splitting.
!
     do uk = k, n - 1

       if ( ar(uk+1,uk) == real(0.0E+00,kind=8) .and. ai(uk+1,uk) == real(0.0E+00,kind=8) ) then
         exit
       end if

     end do
!
!  Compute infinity norm of leading UK by UK (Hessenberg) matrix.
!   
     norm = real(0.0E+00,kind=8)
     mp = 1

     do i = 1, uk

       x = real(0.0E+00,kind=8)
       do j = mp, uk
         x = x + pythag ( ar(i,j), ai(i,j) )
       end do

       norm = max ( norm, x )
       mp = i

     end do
!
!  EPS3 replaces zero pivot in decomposition
!  and close roots are modified by EPS3.
!
     if ( norm == real(0.0E+00,kind=8) ) norm = real(1.0E+00,kind=8)
     eps3 = abs ( norm ) * epsilon ( eps3 )
!
!  GROWTO is the criterion for growth.
!
     ukroot = uk
     ukroot = sqrt ( ukroot )
     growto = 0.1E+00 / ukroot

200  continue

     rlambd = wr(k)
     ilambd = wi(k)
     if ( k == 1 ) go to 280
     km1 = k - 1
     go to 240
!
!  Perturb eigenvalue if it is close to any previous eigenvalue.
!
220  continue

     rlambd = rlambd + eps3

240  continue

     do ii = 1, km1
        i = k - ii
        if ( select(i) .and. abs ( wr(i)-rlambd) < eps3 .and. &
            abs ( wi(i)-ilambd) < eps3 ) then
          go to 220
        end if
     end do

     wr(k) = rlambd
!
!  Form upper Hessenberg (ar,ai)-(rlambd,ilambd) * I
!  and initial complex vector.
!
280  continue

     mp = 1

     do i = 1, uk

        do j = mp, uk
          rm1(i,j) = ar(i,j)
          rm2(i,j) = ai(i,j)
        end do

        rm1(i,i) = rm1(i,i) - rlambd
        rm2(i,i) = rm2(i,i) - ilambd
        mp = i
        rv1(i) = eps3

     end do
!
!  Triangular decomposition with interchanges, replacing zero pivots by eps3.
!
     do i = 2, uk

        mp = i - 1

        if ( pythag ( rm1(i,mp), rm2(i,mp) ) > &
             pythag ( rm1(mp,mp),rm2(mp,mp) ) ) then

          do j = mp, uk
            call r_swap ( rm1(i,j), rm1(mp,j) )
            call r_swap ( rm2(i,j), rm2(mp,j) )
          end do

        end if

        if ( rm1(mp,mp) == real(0.0E+00,kind=8) .and. rm2(mp,mp) == real(0.0E+00,kind=8) ) then
          rm1(mp,mp) = eps3
        end if

        call cdiv ( rm1(i,mp), rm2(i,mp), rm1(mp,mp), rm2(mp,mp), x, y )

        if ( x /= real(0.0E+00,kind=8) .or. y /= real(0.0E+00,kind=8) ) then

          do j = i, uk
            rm1(i,j) = rm1(i,j) - x * rm1(mp,j) + y * rm2(mp,j)
            rm2(i,j) = rm2(i,j) - x * rm2(mp,j) - y * rm1(mp,j)
          end do

        end if

     end do

     if ( rm1(uk,uk) == real(0.0E+00,kind=8) .and. rm2(uk,uk) == real(0.0E+00,kind=8) ) then
       rm1(uk,uk) = eps3
     end if

     its = 0
!
!  Back substitution.
!
  660   continue

    do ii = 1, uk

        i = uk + 1 - ii
        x = rv1(i)
        y = real(0.0E+00,kind=8)

        do j = i+1, uk
          x = x - rm1(i,j) * rv1(j) + rm2(i,j) * rv2(j)
          y = y - rm1(i,j) * rv2(j) - rm2(i,j) * rv1(j)
        end do

        call cdiv ( x, y, rm1(i,i), rm2(i,i), rv1(i), rv2(i) )

     end do
!
!  Acceptance test for eigenvector and normalization.
!
     its = its + 1
     norm = real(0.0E+00,kind=8)
     normv = real(0.0E+00,kind=8)

     do i = 1, uk
        x = pythag ( rv1(i), rv2(i) )
        if ( normv < x ) then
          normv = x
          j = i
        end if
        norm = norm + x
     end do

     if ( norm < growto ) go to 840
!
!  Accept vector.
!
     x = rv1(j)
     y = rv2(j)

     do i = 1, uk
       call cdiv ( rv1(i), rv2(i), x, y, zr(i,s), zi(i,s) )
     end do

     if ( uk == n ) go to 940
     j = uk + 1
     go to 900
!
!  Choose a new starting vector.
!
  840    continue

     if ( its < uk ) then

       x = ukroot
       y = eps3 / ( x + real(1.0E+00,kind=8) )

       rv1(1) = eps3
       rv1(2:uk) = y

       j = uk - its + 1
       rv1(j) = rv1(j) - eps3 * x
       go to 660

     end if
!
!  Error: unaccepted eigenvector.
!
  880    continue

     j = 1
     ierr = -k
!
!  Set remaining vector components to zero.
!
900    continue

       zr(j:n,s) = real(0.0E+00,kind=8)
       zi(j:n,s) = real(0.0E+00,kind=8)

940    continue

       s = s + 1

  end do

  go to 1001
!
!  Set error: underestimate of eigenvector space required.
!
 1000 continue
  if ( ierr /= 0 ) ierr = ierr - n
  if ( ierr == 0 ) ierr = -(2 * n + 1)
 1001 continue
  m = s - 1
  return
end
subroutine combak ( nm, low, igh, ar, ai, int, m, zr, zi )
!
!*******************************************************************************
!
!! COMBAK determines eigenvectors by undoing the COMHES transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a complex general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by COMHES.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH = to the order of the matrix.
!
!    Input, real(kind=8) AR(NM,IGH), AI(NM,IGH), the multipliers which were used in the
!    reduction by COMHES in their lower triangles below the subdiagonal.
!
!    Input, integer INT(IGH), information on the rows and columns interchanged 
!    in the reduction by COMHES.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) ZR(NM,M), ZI(NM,M).  On input, the real(kind=8) and imaginary 
!    parts of the eigenvectors to be back transformed.  On output, the real(kind=8) 
!    and imaginary parts of the transformed eigenvectors.
!
  implicit none
!
  integer igh
  integer m
  integer nm
!
  real(kind=8) ai(nm,igh)
  real(kind=8) ar(nm,igh)
  integer i
  integer int(igh)
  integer j
  integer la
  integer low
  integer mm
  integer mp
  real(kind=8) xi
  real(kind=8) xr
  real(kind=8) zi(nm,m)
  real(kind=8) zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  la = igh - 1

  if ( igh - 1 < low + 1 ) then
    return
  end if

  do mm = low + 1, la

     mp = low + igh - mm

     do i = mp+1, igh

        xr = ar(i,mp-1)
        xi = ai(i,mp-1)

        if ( xr /= real(0.0E+00,kind=8) .or. xi /= real(0.0E+00,kind=8) ) then
          zr(i,1:m) = zr(i,1:m) + xr * zr(mp,1:m) - xi * zi(mp,1:m)
          zi(i,1:m) = zi(i,1:m) + xr * zi(mp,1:m) + xi * zr(mp,1:m)
       end if

     end do

     i = int(mp)

     if ( i /= mp ) then

       do j = 1, m
         call r_swap ( zr(i,j), zr(mp,j) )
         call r_swap ( zi(i,j), zi(mp,j) )
       end do

     end if

  end do

  return
end
subroutine comhes ( nm, n, low, igh, ar, ai, int )
!
!*******************************************************************************
!
!! COMHES transforms a complex general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!    Given a complex general matrix, this subroutine
!    reduces a submatrix situated in rows and columns
!    LOW through IGH to upper Hessenberg form by
!    stabilized elementary similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, real(kind=8) AR(NM,N), AI(NM,N).  On input, the real(kind=8) and imaginary 
!    parts of the complex input matrix.  On output, the real(kind=8) and imaginary 
!    parts of the Hessenberg matrix.  The multipliers which were used in the 
!    reduction are stored in the remaining triangles under the
!    Hessenberg matrix.
!
!    Output, integer INT(IGH), information on the rows and columns
!    interchanged in the reduction.
!
  implicit none
!
  integer igh
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  integer i
  integer int(igh)
  integer j
  integer la
  integer low
  integer m
  real(kind=8) xi
  real(kind=8) xr
  real(kind=8) yi
  real(kind=8) yr
!
  la = igh - 1

  do m = low + 1, la

     xr = real(0.0E+00,kind=8)
     xi = real(0.0E+00,kind=8)
     i = m

     do j = m, igh

       if ( abs ( ar(j,m-1) ) + abs ( ai(j,m-1) ) > &
         abs ( xr ) + abs ( xi ) ) then
         xr = ar(j,m-1)
         xi = ai(j,m-1)
         i = j
       end if

     end do

     int(m) = i
!
!  Interchange rows and columns of AR and AI.
!
     if ( i /= m ) then

       do j = m-1, n
         call r_swap ( ar(i,j), ar(m,j) )
         call r_swap ( ai(i,j), ai(m,j) )
       end do

       do j = 1, igh
         call r_swap ( ar(j,i), ar(j,m) )
         call r_swap ( ai(j,i), ai(j,m) )
       end do

     end if

    if ( xr /= real(0.0E+00,kind=8) .or. xi /= real(0.0E+00,kind=8) ) then

      do i = m+1, igh

        yr = ar(i,m-1)
        yi = ai(i,m-1)

        if ( yr /= real(0.0E+00,kind=8) .or. yi /= real(0.0E+00,kind=8) ) then

          call cdiv ( yr, yi, xr, xi, yr, yi )
          ar(i,m-1) = yr
          ai(i,m-1) = yi

          do j = m, n
            ar(i,j) = ar(i,j) - yr * ar(m,j) + yi * ai(m,j)
            ai(i,j) = ai(i,j) - yr * ai(m,j) - yi * ar(m,j)
          end do

          ar(1:igh,m) = ar(1:igh,m) + yr * ar(1:igh,i) - yi * ai(1:igh,i)
          ai(1:igh,m) = ai(1:igh,m) + yr * ai(1:igh,i) + yi * ar(1:igh,i)

        end if

      end do

    end if

  end do

  return
end
subroutine comlr ( nm, n, low, igh, hr, hi, wr, wi, ierr )
!
!*******************************************************************************
!
!! COMLR gets all eigenvalues of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a complex upper Hessenberg 
!    matrix by the modified LR method.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, real(kind=8) HR(NM,N), HI(NM,N).  On input, the real(kind=8) and imaginary 
!    parts of the complex upper Hessenberg matrix.  Their lower triangles 
!    below the subdiagonal contain the multipliers which were used in the 
!    reduction by COMHES if performed.  On output, the upper Hessenberg 
!    portions of HR and HI have been destroyed.  Therefore, they must be 
!    saved before calling COMLR if subsequent calculation of eigenvectors 
!    is to be performed.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts of the 
!    eigenvalues.  If an error exit is made, the eigenvalues should be correct
!    for indices IERR+1,...,N.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the limit of 30*N iterations is exhausted while the J-th 
!      eigenvalue is being sought.
!
  implicit none
!
  integer n
  integer nm
!
  integer en
  integer enm1
  real(kind=8) hi(nm,n)
  real(kind=8) hr(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer l
  integer ll
  integer low
  integer m
  integer mm
  real(kind=8) si
  real(kind=8) sr
  real(kind=8) ti
  real(kind=8) tr
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) xi
  real(kind=8) xr
  real(kind=8) yi
  real(kind=8) yr
  real(kind=8) zzi
  real(kind=8) zzr
!
  ierr = 0
!
!  Store roots isolated by CBAL.
!
  do i = 1, n
    if ( i < low .or. i > igh ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = real(0.0E+00,kind=8)
  ti = real(0.0E+00,kind=8)
  itn = 30 * n
!
!  Search for next eigenvalue.
!
  220 continue

  if ( en < low ) then
    return
  end if

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
  240 continue

  do ll = low, en
     l = en + low - ll
     if ( l == low ) go to 300
     tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) & 
       + abs ( hi(l,l) )
     tst2 = tst1 + abs ( hr(l,l-1) ) + abs ( hi(l,l-1) )
     if ( tst2 == tst1) go to 300
  end do
!
!  Form shift.
!
300 continue

  if ( l == en ) go to 660

  if ( itn == 0 ) then
    ierr = en
    return
  end if

  if ( its == 10 .or. its == 20 ) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
  xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
  if ( xr == real(0.0E+00,kind=8) .and. xi == real(0.0E+00,kind=8) ) go to 340
  yr = (hr(enm1,enm1) - sr) / real(2.0E+00,kind=8)
  yi = (hi(enm1,enm1) - si) / real(2.0E+00,kind=8)
  call csroot ( yr**2-yi**2+xr, real(2.0E+00,kind=8)*yr*yi+xi, zzr, zzi )

  if ( yr * zzr + yi * zzi < real(0.0E+00,kind=8) ) then
    zzr = -zzr
    zzi = -zzi
  end if

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
  320 continue

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = abs ( hi(en,enm1) ) + abs ( hi(enm1,en-2) )

  340 continue

  do i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  xr = abs ( hr(enm1,enm1) ) + abs ( hi(enm1,enm1) )
  yr = abs ( hr(en,enm1) ) + abs ( hi(en,enm1) )
  zzr = abs ( hr(en,en) ) + abs ( hi(en,en) )

  do mm = l, enm1
    m = enm1 + l - mm
    if ( m == l ) then
      exit
    end if
    yi = yr
    yr = abs ( hr(m,m-1) ) + abs ( hi(m,m-1) )
    xi = zzr
    zzr = xr
    xr = abs ( hr(m-1,m-1) ) + abs ( hi(m-1,m-1) )
    tst1 = zzr / yi * (zzr + xr + xi)
    tst2 = tst1 + yr
    if ( tst2 == tst1 ) then
      exit
    end if
  end do
!
!  Triangular decomposition H=L*R.
!
  do i = m+1, en

     xr = hr(i-1,i-1)
     xi = hi(i-1,i-1)
     yr = hr(i,i-1)
     yi = hi(i,i-1)
     if ( abs ( xr ) + abs ( xi ) >= abs ( yr ) + abs ( yi ) ) go to 460
!
!  Interchange rows of HR and HI.
!
     do j = i-1, en
       call r_swap ( hr(i-1,j), hr(i,j) )
       call r_swap ( hi(i-1,j), hi(i,j) )
     end do

     call cdiv ( xr, xi, yr, yi, zzr, zzi )
     wr(i) = real(1.0E+00,kind=8)
     go to 480

460 continue

     call cdiv ( yr, yi, xr, xi, zzr, zzi )
     wr(i) = -real(1.0E+00,kind=8)

480  continue

     hr(i,i-1) = zzr
     hi(i,i-1) = zzi

     do j = i, en
        hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
        hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
     end do

  end do
!
!  Composition R*L=H.
!
  do j = m+1, en

    xr = hr(j,j-1)
    xi = hi(j,j-1)
    hr(j,j-1) = real(0.0E+00,kind=8)
    hi(j,j-1) = real(0.0E+00,kind=8)
!
!  Interchange columns of HR and HI, if necessary.
!
    if ( wr(j) > real(0.0E+00,kind=8) ) then

      do i = l, j
        call r_swap ( hr(i,j-1), hr(i,j) )
        call r_swap ( hi(i,j-1), hi(i,j) )
      end do

    end if

    do i = l, j
      hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
      hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
    end do

  end do

  go to 240
!
!  A root found.
!
  660 continue

  wr(en) = hr(en,en) + tr
  wi(en) = hi(en,en) + ti
  en = enm1
  go to 220
end
subroutine comlr2 ( nm, n, low, igh, int, hr, hi, wr, wi, zr, zi, ierr )
!
!*******************************************************************************
!
!! COMLR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors of a complex 
!    upper Hessenberg matrix by the modified LR method.  The eigenvectors 
!    of a complex general matrix can also be found if COMHES has been used 
!    to reduce this general matrix to Hessenberg form.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input, integer INT(IGH), information on the rows and columns interchanged
!    in the reduction by COMHES, if performed.  If the eigenvectors of the 
!    Hessenberg matrix are desired, set INT(J)=J for these elements.
!
!    Input/output, real(kind=8) HR(NM,N), HI(NM,N).  On input, the real(kind=8) and imaginary
!    parts of the complex upper Hessenberg matrix.  Their lower triangles 
!    below the subdiagonal contain the multipliers which were used in the 
!    reduction by COMHES, if performed.  If the eigenvectors of the Hessenberg
!    matrix are desired, these elements must be set to zero.  On output,
!    the upper Hessenberg portions of HR and HI have been destroyed, but the
!    location HR(1,1) contains the norm of the triangularized matrix.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts of the 
!    eigenvalues.  If an error exit is made, the eigenvalues should be 
!    correct for indices IERR+1,...,N.
!
!    Output, real(kind=8) ZR(NM,N), ZI(NM,N), the real(kind=8) and imaginary parts of the
!    eigenvectors.  The eigenvectors are unnormalized.  If an error exit 
!    is made, none of the eigenvectors has been found.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the limit of 30*N iterations is exhausted while the J-th 
!      eigenvalue is being sought.
!
  implicit none
!
  integer n
  integer nm
!
  integer en
  integer enm1
  real(kind=8) hi(nm,n)
  real(kind=8) hr(nm,n)
  integer i
  integer iend
  integer ierr
  integer igh
  integer ii
  integer int(igh)
  integer itn
  integer its
  integer j
  integer jj
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer nn
  real(kind=8) norm
  real(kind=8) si
  real(kind=8) sr
  real(kind=8) ti
  real(kind=8) tr
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) xi
  real(kind=8) xr
  real(kind=8) yi
  real(kind=8) yr
  real(kind=8) zi(nm,n)
  real(kind=8) zr(nm,n)
  real(kind=8) zzi
  real(kind=8) zzr
!
  ierr = 0
!
!  Initialize the eigenvector matrix.
!
  call rmat_ident ( nm, n, zr )

  zi(1:n,1:n) = real(0.0E+00,kind=8)
!
!  Form the matrix of accumulated transformations from the information left 
!  by COMHES.
!
  iend = igh - low - 1

  do ii = 1, iend

    i = igh - ii

    do k = i+1, igh
      zr(k,i) = hr(k,i-1)
      zi(k,i) = hi(k,i-1)
    end do

    j = int(i)

    if ( i /= j ) then

      do k = i, igh
        zr(i,k) = zr(j,k)
        zi(i,k) = zi(j,k)
        zr(j,k) = real(0.0E+00,kind=8)
        zi(j,k) = real(0.0E+00,kind=8)
      end do

      zr(j,i) = real(1.0E+00,kind=8)

    end if

  end do
!
!  Store roots isolated by CBAL.
!
  do i = 1, n
    if ( i < low .or. i > igh ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = real(0.0E+00,kind=8)
  ti = real(0.0E+00,kind=8)
  itn = 30 * n
!
!  Search for next eigenvalue.
!
220 continue

  if ( en < low ) go to 680

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
  240 continue

  do ll = low, en

     l = en + low - ll

     if ( l == low ) then
       exit
     end if

     tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
       + abs ( hi(l,l) )
     tst2 = tst1 + abs ( hr(l,l-1) ) + abs ( hi(l,l-1) )

     if ( tst2 == tst1 ) then
       exit
     end if

  end do
!
!  Form shift.
!
  if ( l == en ) go to 660
  if ( itn == 0 ) go to 1000
  if ( its == 10 .or. its == 20 ) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
  xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
  if ( xr == real(0.0E+00,kind=8) .and. xi == real(0.0E+00,kind=8) ) go to 340
  yr = (hr(enm1,enm1) - sr) / real(2.0E+00,kind=8)
  yi = (hi(enm1,enm1) - si) / real(2.0E+00,kind=8)
  call csroot ( yr**2-yi**2+xr, real(2.0E+00,kind=8)*yr*yi+xi, zzr, zzi )

  if ( yr * zzr + yi * zzi < real(0.0E+00,kind=8) ) then
    zzr = -zzr
    zzi = -zzi
  end if

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
  320 continue

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = abs ( hi(en,enm1) ) + abs ( hi(enm1,en-2) )

  340 continue

  do i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  xr = abs ( hr(enm1,enm1) ) + abs ( hi(enm1,enm1) )
  yr = abs ( hr(en,enm1) ) + abs ( hi(en,enm1) )
  zzr = abs ( hr(en,en) ) + abs ( hi(en,en) )

  do mm = l, enm1
     m = enm1 + l - mm
     if ( m == l ) then
       exit
     end if
     yi = yr
     yr = abs ( hr(m,m-1) ) + abs ( hi(m,m-1) )
     xi = zzr
     zzr = xr
     xr = abs ( hr(m-1,m-1) ) + abs ( hi(m-1,m-1) )
     tst1 = zzr / yi * (zzr + xr + xi)
     tst2 = tst1 + yr
     if ( tst2 == tst1 ) then
       exit
     end if
  end do
!
!  Triangular decomposition H=L*R.
!
  do i = m+1, en

     xr = hr(i-1,i-1)
     xi = hi(i-1,i-1)
     yr = hr(i,i-1)
     yi = hi(i,i-1)
     if ( abs ( xr ) + abs ( xi) >= abs ( yr ) + abs ( yi ) ) go to 460
!
!  Interchange rows of HR and HI.
!
     do j = i-1, n
       call r_swap ( hr(i-1,j), hr(i,j) )
       call r_swap ( hi(i-1,j), hi(i,j) )
    end do

     call cdiv ( xr, xi, yr, yi, zzr, zzi )
     wr(i) = real(1.0E+00,kind=8)
     go to 480
460  continue

     call cdiv ( yr, yi, xr, xi, zzr, zzi )
     wr(i) = -real(1.0E+00,kind=8)

480  continue

     hr(i,i-1) = zzr
     hi(i,i-1) = zzi

     do j = i, n
       hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
       hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
     end do

  end do
!
!  Composition R*L=H.
!
  do j = m+1, en

     xr = hr(j,j-1)
     xi = hi(j,j-1)
     hr(j,j-1) = real(0.0E+00,kind=8)
     hi(j,j-1) = real(0.0E+00,kind=8)
!
!  Interchange columns of HR, HI, ZR, and ZI.
!
     if ( wr(j) > real(0.0E+00,kind=8) ) then

       do i = 1, j
         call r_swap ( hr(i,j-1), hr(i,j) )
         call r_swap ( hi(i,j-1), hi(i,j) )
       end do

       do i = low, igh
         call r_swap ( zr(i,j-1), zr(i,j) )
         call r_swap ( zi(i,j-1), zi(i,j) )
       end do

    end if

    do i = 1, j
      hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
      hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
    end do
!
!  Accumulate transformations.
!
    do i = low, igh
      zr(i,j-1) = zr(i,j-1) + xr * zr(i,j) - xi * zi(i,j)
      zi(i,j-1) = zi(i,j-1) + xr * zi(i,j) + xi * zr(i,j)
    end do

  end do

  go to 240
!
!  A root found.
!
  660 continue

  hr(en,en) = hr(en,en) + tr
  wr(en) = hr(en,en)
  hi(en,en) = hi(en,en) + ti
  wi(en) = hi(en,en)
  en = enm1
  go to 220
!
!  All roots found.  
!  Backsubstitute to find vectors of upper triangular form.
!
  680 continue

  norm = real(0.0E+00,kind=8)

  do i = 1, n
    do j = i, n
      tr = abs ( hr(i,j) ) + abs ( hi(i,j) )
      if ( tr > norm ) norm = tr
    end do
  end do

  hr(1,1) = norm
  if ( n == 1 ) then
    return
  end if

  if ( norm == real(0.0E+00,kind=8) ) then
    return
  end if

  do nn = 2, n

     en = n + 2 - nn
     xr = wr(en)
     xi = wi(en)
     hr(en,en) = real(1.0E+00,kind=8)
     hi(en,en) = real(0.0E+00,kind=8)
     enm1 = en - 1

     do ii = 1, enm1

        i = en - ii
        zzr = real(0.0E+00,kind=8)
        zzi = real(0.0E+00,kind=8)

        do j = i+1, en
          zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
          zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
        end do

        yr = xr - wr(i)
        yi = xi - wi(i)

        if ( yr == real(0.0E+00,kind=8) .and. yi == real(0.0E+00,kind=8) ) then

          tst1 = norm
          yr = tst1

          do
            yr = 0.01E+00 * yr
            tst2 = norm + yr
            if ( tst2 <=  tst1 ) then
              exit
            end if
          end do

        end if

        call cdiv ( zzr, zzi, yr, yi, hr(i,en), hi(i,en) )
!
!  Overflow control.
!
        tr = abs ( hr(i,en) ) + abs ( hi(i,en) )

        if ( tr /= real(0.0E+00,kind=8) ) then

          tst1 = tr
          tst2 = tst1 + real(1.0E+00,kind=8) / tst1

          if ( tst2 <= tst1 ) then

            hr(i:en,en) = hr(i:en,en) / tr
            hi(i:en,en) = hi(i:en,en) / tr

          end if

        end if

      end do

  end do
!
!  End backsubstitution.
!
  enm1 = n - 1
!
!  Vectors of isolated roots.
!
  do i = 1, n - 1

    if ( i < low .or. i > igh ) then

      zr(i,i+1:n) = hr(i,i+1:n)
      zi(i,i+1:n) = hi(i,i+1:n)

    end if

  end do
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  do jj = low, n - 1

    j = n + low - jj
    m = min ( j, igh )

    do i = low, igh
      zzr = real(0.0E+00,kind=8)
      zzi = real(0.0E+00,kind=8)
      do k = low, m
        zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
        zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
      end do
      zr(i,j) = zzr
      zi(i,j) = zzi
    end do

  end do

  return
!
!  Set error: all eigenvalues have not converged after 30*N iterations.
!
 1000 continue

  ierr = en
  return
end
subroutine comqr ( nm, n, low, igh, hr, hi, wr, wi, ierr )
!
!*******************************************************************************
!
!! COMQR gets eigenvalues of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a complex
!    upper Hessenberg matrix by the QR method.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, real(kind=8) HR(NM,N), HI(NM,N).  On input, the real(kind=8) and imaginary 
!    parts of the complex upper Hessenberg matrix.  Their lower triangles 
!    below the subdiagonal contain information about the unitary 
!    transformations used in the reduction by CORTH, if performed.  On output,
!    the upper Hessenberg portions of HR and HI have been destroyed.  
!    Therefore, they must be saved before calling COMQR if subsequent 
!    calculation of eigenvectors is to be performed.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts of the 
!    eigenvalues.  If an error exit is made, the eigenvalues should be 
!    correct for indices IERR+1,...,N.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the limit of 30*N iterations is exhausted while the J-th 
!       eigenvalue is being sought.
!
  implicit none
!
  integer n
  integer nm
!
  integer en
  integer enm1
  real(kind=8) hi(nm,n)
  real(kind=8) hr(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer l
  integer ll
  integer low
  real(kind=8) norm
  real(kind=8) pythag
  real(kind=8) si
  real(kind=8) sr
  real(kind=8) ti
  real(kind=8) tr
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) xi
  real(kind=8) xr
  real(kind=8) yi
  real(kind=8) yr
  real(kind=8) zzi
  real(kind=8) zzr
!
  ierr = 0
!
!  Create real(kind=8) subdiagonal elements.
!
  l = low + 1

  do i = l, igh

     ll = min ( i+1, igh )

     if ( hi(i,i-1) /= real(0.0E+00,kind=8) ) then

     norm = pythag ( hr(i,i-1), hi(i,i-1) )
     yr = hr(i,i-1) / norm
     yi = hi(i,i-1) / norm
     hr(i,i-1) = norm
     hi(i,i-1) = real(0.0E+00,kind=8)

     do j = i, igh
       si = yr * hi(i,j) - yi * hr(i,j)
       hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
       hi(i,j) = si
     end do

     do j = low, ll
       si = yr * hi(j,i) + yi * hr(j,i)
       hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
       hi(j,i) = si
     end do

    end if

  end do
!
!  Store roots isolated by CBAL.
!
  do i = 1, n
    if ( i < low .or. i > igh ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = real(0.0E+00,kind=8)
  ti = real(0.0E+00,kind=8)
  itn = 30 * n
!
!  Search for next eigenvalue.
!
  220 continue

  if ( en < low ) then
    return
  end if

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
  240 continue

  do ll = low, en
    l = en + low - ll
    if ( l == low ) then
      exit
    end if
    tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
      + abs ( hi(l,l) )
    tst2 = tst1 + abs ( hr(l,l-1) )
    if ( tst2 == tst1 ) then
      exit
    end if
  end do
!
!  Form shift.
!
  if ( l == en ) go to 660

  if ( itn == 0 ) go to 1000

  if ( its == 10 .or. its == 20 ) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1)
  xi = hi(enm1,en) * hr(en,enm1)
  if ( xr == real(0.0E+00,kind=8) .and. xi == real(0.0E+00,kind=8) ) go to 340
  yr = (hr(enm1,enm1) - sr) / real(2.0E+00,kind=8)
  yi = (hi(enm1,enm1) - si) / real(2.0E+00,kind=8)

  call csroot ( yr**2-yi**2+xr, real(2.0E+00,kind=8)*yr*yi+xi, zzr, zzi )

  if ( yr * zzr + yi * zzi < real(0.0E+00,kind=8) ) then
    zzr = -zzr
    zzi = -zzi
  end if

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
320 continue

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = real(0.0E+00,kind=8)

340 continue

  do i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Reduce to triangle (rows).
!
  do i = l+1, en

     sr = hr(i,i-1)
     hr(i,i-1) = real(0.0E+00,kind=8)
     norm = pythag ( pythag ( hr(i-1,i-1), hi(i-1,i-1) ), sr )
     xr = hr(i-1,i-1) / norm
     wr(i-1) = xr
     xi = hi(i-1,i-1) / norm
     wi(i-1) = xi
     hr(i-1,i-1) = norm
     hi(i-1,i-1) = real(0.0E+00,kind=8)
     hi(i,i-1) = sr / norm

     do j = i, en
        yr = hr(i-1,j)
        yi = hi(i-1,j)
        zzr = hr(i,j)
        zzi = hi(i,j)
        hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
        hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
        hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
        hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
    end do

  end do

  si = hi(en,en)

  if ( si /= real(0.0E+00,kind=8) ) then
    norm = pythag ( hr(en,en), si )
    sr = hr(en,en) / norm
    si = si / norm
    hr(en,en) = norm
    hi(en,en) = real(0.0E+00,kind=8)
  end if
!
!  Inverse operation (columns).
!
  do j = l+1, en

     xr = wr(j-1)
     xi = wi(j-1)

     do i = l, j

        yr = hr(i,j-1)
        yi = real(0.0E+00,kind=8)
        zzr = hr(i,j)
        zzi = hi(i,j)
        if ( i /= j ) then
          yi = hi(i,j-1)
          hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
        end if
        hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
        hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
        hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

     end do

  end do

  if ( si /= real(0.0E+00,kind=8) ) then

    do i = l, en
      yr = hr(i,en)
      yi = hi(i,en)
      hr(i,en) = sr * yr - si * yi
      hi(i,en) = sr * yi + si * yr
    end do

  end if

  go to 240
!
!  A root found.
!
660 continue

  wr(en) = hr(en,en) + tr
  wi(en) = hi(en,en) + ti
  en = enm1
  go to 220
!
!  Set error: all eigenvalues have not converged after 30*n iterations.
!
1000 continue

  ierr = en
  return
end
subroutine comqr2 ( nm, n, low, igh, ortr, orti, hr, hi, wr, wi, zr, zi, ierr )
!
!*******************************************************************************
!
!! COMQR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors
!    of a complex upper Hessenberg matrix by the QR
!    method.  The eigenvectors of a complex general matrix
!    can also be found if CORTH has been used to reduce
!    this general matrix to Hessenberg form.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, real(kind=8) ORTR(N), ORTI(N).  On input, information about the 
!    unitary transformations used in the reduction by CORTH, if performed.
!    If the eigenvectors of the Hessenberg matrix are desired, set ORTR(J) and
!    ORTI(J) to real(0.0E+00,kind=8) for these elements.  On output, these arrays
!    have been overwritten.
!
!    Input/output, real(kind=8) HR(NM,N), HI(NM,N).  On input, the real(kind=8) and imaginary 
!    parts of the complex upper Hessenberg matrix.  Their lower triangles 
!    below the subdiagonal contain further information about the 
!    transformations which were used in the reduction by CORTH, if performed.  
!    If the eigenvectors of the Hessenberg matrix are desired, these elements 
!    may be arbitrary.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts of the 
!    eigenvalues.  If an error exit is made, the eigenvalues should be 
!    correct for indices IERR+1,...,N.
!
!    Output, real(kind=8) ZR(NM,N), ZI(NM,N), the real(kind=8) and imaginary parts of the 
!    eigenvectors.  The eigenvectors are unnormalized.  If an error exit 
!    is made, none of the eigenvectors has been found.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the limit of 30*N iterations is exhausted while the J-th 
!      eigenvalue is being sought.
!
  implicit none
!
  integer igh
  integer n
  integer nm
!
  integer en
  integer enm1
  real(kind=8) hi(nm,n)
  real(kind=8) hr(nm,n)
  integer i
  integer iend
  integer ierr
  integer ii
  integer itn
  integer its
  integer j
  integer jj
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer nn
  real(kind=8) norm
  real(kind=8) orti(igh)
  real(kind=8) ortr(igh)
  real(kind=8) pythag
  real(kind=8) si
  real(kind=8) sr
  real(kind=8) ti
  real(kind=8) tr
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) xi
  real(kind=8) xr
  real(kind=8) yi
  real(kind=8) yr
  real(kind=8) zi(nm,n)
  real(kind=8) zr(nm,n)
  real(kind=8) zzi 
  real(kind=8) zzr
!
  ierr = 0
!
!  Initialize eigenvector matrix.
!
  call rmat_ident ( nm, n, zr )

  zi(1:n,1:n) = real(0.0E+00,kind=8)
!
!  Form the matrix of accumulated transformations from the information 
!  left by CORTH.
!
  iend = igh - low - 1
  if ( iend ) 180, 150, 105

105 continue

  do ii = 1, iend

     i = igh - ii
     if ( ortr(i) == real(0.0E+00,kind=8) .and. orti(i) == real(0.0E+00,kind=8) ) go to 140
     if ( hr(i,i-1) == real(0.0E+00,kind=8) .and. hi(i,i-1) == real(0.0E+00,kind=8) ) go to 140
!
!  Norm below is negative of H formed in CORTH.
!
     norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)

     do k = i+1, igh
       ortr(k) = hr(k,i-1)
       orti(k) = hi(k,i-1)
     end do

     do j = i, igh

        sr = real(0.0E+00,kind=8)
        si = real(0.0E+00,kind=8)

        do k = i, igh
          sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
          si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
        end do

        sr = sr / norm
        si = si / norm

        do k = i, igh
          zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
          zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
        end do

      end do

140 continue

  end do
!
!  Create real(kind=8) subdiagonal elements.
!
150 continue

  l = low + 1

  do i = l, igh

     ll = min ( i+1, igh )

     if ( hi(i,i-1) == real(0.0E+00,kind=8) ) go to 170

     norm = pythag ( hr(i,i-1), hi(i,i-1) )
     yr = hr(i,i-1) / norm
     yi = hi(i,i-1) / norm
     hr(i,i-1) = norm
     hi(i,i-1) = real(0.0E+00,kind=8)

     do j = i, n
       si = yr * hi(i,j) - yi * hr(i,j)
       hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
       hi(i,j) = si
     end do

     do j = 1, ll
       si = yr * hi(j,i) + yi * hr(j,i)
       hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
       hi(j,i) = si
     end do

     do j = low, igh
       si = yr * zi(j,i) + yi * zr(j,i)
       zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
       zi(j,i) = si
     end do

170 continue

  end do
!
!  Store roots isolated by CBAL.
!
180 continue

  do i = 1, n
    if ( i < low .or. i > igh) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = real(0.0E+00,kind=8)
  ti = real(0.0E+00,kind=8)
  itn = 30 * n
!
!  Search for next eigenvalue.
!
220 continue

  if ( en < low ) go to 680
  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
240 continue

  do ll = low, en
    l = en + low - ll
    if ( l == low ) then
      exit
    end if
    tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
      + abs ( hi(l,l) )
    tst2 = tst1 + abs ( hr(l,l-1) )
    if ( tst2 == tst1 ) then
      exit
    end if
  end do
!
!  Form shift.
!
  if ( l == en ) go to 660
  if ( itn == 0 ) go to 1000
  if ( its == 10 .or. its == 20 ) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1)
  xi = hi(enm1,en) * hr(en,enm1)
  if ( xr == real(0.0E+00,kind=8) .and. xi == real(0.0E+00,kind=8) ) go to 340
  yr = ( hr(enm1,enm1) - sr ) / real(2.0E+00,kind=8)
  yi = ( hi(enm1,enm1) - si ) / real(2.0E+00,kind=8)

  call csroot ( yr**2-yi**2+xr, real(2.0E+00,kind=8)*yr*yi+xi, zzr, zzi )

  if ( yr * zzr + yi * zzi < real(0.0E+00,kind=8) ) then
    zzr = -zzr
    zzi = -zzi
  end if

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
320 continue

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = real(0.0E+00,kind=8)

340 continue

  do i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Reduce to triangle (rows).
!
  do i = l+1, en

     sr = hr(i,i-1)
     hr(i,i-1) = real(0.0E+00,kind=8)
     norm = pythag ( pythag ( hr(i-1,i-1), hi(i-1,i-1) ), sr )
     xr = hr(i-1,i-1) / norm
     wr(i-1) = xr
     xi = hi(i-1,i-1) / norm
     wi(i-1) = xi
     hr(i-1,i-1) = norm
     hi(i-1,i-1) = real(0.0E+00,kind=8)
     hi(i,i-1) = sr / norm

     do j = i, n
        yr = hr(i-1,j)
        yi = hi(i-1,j)
        zzr = hr(i,j)
        zzi = hi(i,j)
        hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
        hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
        hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
        hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
     end do

  end do

  si = hi(en,en)

  if ( si /= real(0.0E+00,kind=8) ) then

    norm = pythag ( hr(en,en), si )
    sr = hr(en,en) / norm
    si = si / norm
    hr(en,en) = norm
    hi(en,en) = real(0.0E+00,kind=8)

    do j = en+1, n
      yr = hr(en,j)
      yi = hi(en,j)
      hr(en,j) = sr * yr + si * yi
      hi(en,j) = sr * yi - si * yr
    end do

  end if
!
!  Inverse operation (columns).
!
  do j = l+1, en

     xr = wr(j-1)
     xi = wi(j-1)

     do i = 1, j

       yr = hr(i,j-1)
       yi = real(0.0E+00,kind=8)
       zzr = hr(i,j)
       zzi = hi(i,j)

       if ( i /= j ) then
         yi = hi(i,j-1)
         hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
       end if

       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
       hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
       hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

     end do

     do i = low, igh
       yr = zr(i,j-1)
       yi = zi(i,j-1)
       zzr = zr(i,j)
       zzi = zi(i,j)
       zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
       zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
       zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
       zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
     end do

  end do

  if ( si /= real(0.0E+00,kind=8) ) then

    do i = 1, en
      yr = hr(i,en)
      yi = hi(i,en)
      hr(i,en) = sr * yr - si * yi
      hi(i,en) = sr * yi + si * yr
    end do

    do i = low, igh
      yr = zr(i,en)
      yi = zi(i,en)
      zr(i,en) = sr * yr - si * yi
      zi(i,en) = sr * yi + si * yr
    end do

  end if

  go to 240
!
!  A root found.
!
660 continue

  hr(en,en) = hr(en,en) + tr
  wr(en) = hr(en,en)
  hi(en,en) = hi(en,en) + ti
  wi(en) = hi(en,en)
  en = enm1
  go to 220
!
!  All roots found.  
!  Backsubstitute to find vectors of upper triangular form.
!
680 continue

  norm = real(0.0E+00,kind=8)

  do i = 1, n
    do j = i, n
      tr = abs ( hr(i,j) ) + abs ( hi(i,j) )
      norm = max ( norm, tr )
    end do
  end do

  if ( n == 1 ) then
    return
  end if

  if ( norm == real(0.0E+00,kind=8) ) then
    return
  end if

  do nn = 2, n

     en = n + 2 - nn
     xr = wr(en)
     xi = wi(en)
     hr(en,en) = real(1.0E+00,kind=8)
     hi(en,en) = real(0.0E+00,kind=8)
     enm1 = en - 1

     do ii = 1, enm1

        i = en - ii
        zzr = real(0.0E+00,kind=8)
        zzi = real(0.0E+00,kind=8)

        do j = i+1, en
          zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
          zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
        end do

        yr = xr - wr(i)
        yi = xi - wi(i)

        if ( yr == real(0.0E+00,kind=8) .and. yi == real(0.0E+00,kind=8) ) then

           tst1 = norm
           yr = tst1
           do
             yr = 0.01E+00 * yr
             tst2 = norm + yr
             if ( tst2 <= tst1 ) then
               exit
             end if
           end do

        end if

        call cdiv ( zzr, zzi, yr, yi, hr(i,en), hi(i,en) )
!
!  Overflow control.
!
        tr = abs ( hr(i,en) ) + abs ( hi(i,en) )

        if ( tr /= real(0.0E+00,kind=8) ) then

          tst1 = tr
          tst2 = tst1 + real(1.0E+00,kind=8) / tst1

          if ( tst2 <= tst1 ) then

            do j = i, en
              hr(j,en) = hr(j,en)/tr
              hi(j,en) = hi(j,en)/tr
            end do

          end if

       end if

     end do

  end do
!
!  End backsubstitution.
!
  enm1 = n - 1
!
!  Vectors of isolated roots.
!
  do i = 1, n - 1

    if ( i < low .or. i > igh ) then

      do j = i+1, n
        zr(i,j) = hr(i,j)
        zi(i,j) = hi(i,j)
      end do

    end if

  end do
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  do jj = low, n - 1

     j = n + low - jj
     m = min ( j, igh )

     do i = low, igh

        zzr = real(0.0E+00,kind=8)
        zzi = real(0.0E+00,kind=8)
        do k = low, m
          zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
          zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
        end do

        zr(i,j) = zzr
        zi(i,j) = zzi

      end do

  end do

  return
!
!  Set error: all eigenvalues have not converged after 30*n iterations.
!
1000 continue

  ierr = en
  return
end
subroutine cortb ( nm, low, igh, ar, ai, ortr, orti, m, zr, zi )
!
!*******************************************************************************
!
!! CORTB determines eigenvectors by undoing the CORTH transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a complex general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by CORTH.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH to the order of the matrix.
!
!    Input, real(kind=8) AR(NM,IGH), AI(NM,IGH), information about the unitary
!    transformations used in the reduction by CORTH in their strict lower 
!    triangles.
!
!    Input/output, real(kind=8) ORTR(IGH), ORTI(IGH).  On input, further information 
!    about the transformations used in the reduction by CORTH.  On output,
!    ORTR and ORTI have been further altered.
!
!    Input, integer M, the number of columns of ZR and ZI to be back
!    transformed.
!
!    Input/output, real(kind=8) ZR(NM,M), ZI(NM,M).  On input, the real(kind=8) and imaginary 
!    parts of the eigenvectors to be back transformed.  On output, the real(kind=8) 
!    and imaginary parts of the transformed eigenvectors.
!
  implicit none
!
  integer igh
  integer m
  integer nm
!
  real(kind=8) ai(nm,igh)
  real(kind=8) ar(nm,igh)
  real(kind=8) gi
  real(kind=8) gr
  real(kind=8) h
  integer i
  integer j
  integer la
  integer low
  integer mm
  integer mp
  real(kind=8) orti(igh)
  real(kind=8) ortr(igh)
  real(kind=8) zi(nm,m)
  real(kind=8) zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  la = igh - 1

  if ( igh - 1 < low + 1 ) then
    return
  end if

  do mm = low + 1, la

    mp = low + igh - mm

    if ( ar(mp,mp-1) /= real(0.0E+00,kind=8) .or. ai(mp,mp-1) /= real(0.0E+00,kind=8) ) then

      h = ar(mp,mp-1) * ortr(mp) + ai(mp,mp-1) * orti(mp)

      ortr(mp+1:igh) = ar(mp+1:igh,mp-1)
      orti(mp+1:igh) = ai(mp+1:igh,mp-1)

      do j = 1, m

        gr = ( dot_product ( ortr(mp:igh), zr(mp:igh,j) ) &
             + dot_product ( orti(mp:igh), zi(mp:igh,j) ) ) / h

        gi = ( dot_product ( ortr(mp:igh), zi(mp:igh,j) ) &
             - dot_product ( orti(mp:igh), zr(mp:igh,j) ) ) / h

        do i = mp, igh
          zr(i,j) = zr(i,j) + gr * ortr(i) - gi * orti(i)
          zi(i,j) = zi(i,j) + gr * orti(i) + gi * ortr(i)
        end do

      end do

    end if

  end do

  return
end
subroutine corth ( nm, n, low, igh, ar, ai, ortr, orti )
!
!*******************************************************************************
!
!! CORTH transforms a complex general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!    Given a complex general matrix, this subroutine
!    reduces a submatrix situated in rows and columns
!    LOW through IGH to upper Hessenberg form by
!    unitary similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine CBAL.
!    If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, real(kind=8) AR(NM,N), AI(NM,N).  On input, the real(kind=8) and imaginary 
!    parts of the complex input matrix.  On output, the real(kind=8) and imaginary
!    parts of the Hessenberg matrix.  Information about the unitary 
!    transformations used in the reduction is stored in the remaining 
!    triangles under the Hessenberg matrix.
!
!    Output, real(kind=8) ORTR(IGH), ORTI(IGH), further information about the
!    transformations.
!
  implicit none
!
  integer igh
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  real(kind=8) f
  real(kind=8) fi
  real(kind=8) fr
  real(kind=8) g
  real(kind=8) h
  integer i
  integer ii
  integer j
  integer jj
  integer la
  integer m,mp,low
  real(kind=8) orti(igh)
  real(kind=8) ortr(igh)
  real(kind=8) pythag
  real(kind=8) scale
!
  la = igh - 1

  if ( igh - 1 < low + 1 ) then
    return
  end if

  do m = low + 1, la

    h = real(0.0E+00,kind=8)
    ortr(m) = real(0.0E+00,kind=8)
    orti(m) = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)
!
!  Scale column.
!
    do i = m, igh
      scale = scale + abs ( ar(i,m-1) ) + abs ( ai(i,m-1) )
    end do

    if ( scale == real(0.0E+00,kind=8) ) then
      cycle
    end if

    mp = m + igh

    do ii = m, igh
      i = mp - ii
      ortr(i) = ar(i,m-1) / scale
      orti(i) = ai(i,m-1) / scale
      h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
    end do

    g = sqrt ( h )
    f = pythag ( ortr(m), orti(m) )

    if ( f /= real(0.0E+00,kind=8) ) then
      h = h + f * g
      g = g / f
      ortr(m) = ( real(1.0E+00,kind=8) + g ) * ortr(m)
      orti(m) = ( real(1.0E+00,kind=8) + g ) * orti(m)
    else
      ortr(m) = g
      ar(m,m-1) = scale
    end if
!
!  Form (I-(U*Ut)/h) * A.
!
    do j = m, n

      fr = real(0.0E+00,kind=8)
      fi = real(0.0E+00,kind=8)

      do ii = m, igh
        i = mp - ii
        fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
        fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
      end do

      fr = fr / h
      fi = fi / h

      ar(m:igh,j) = ar(m:igh,j) - fr * ortr(m:igh) + fi * orti(m:igh)
      ai(m:igh,j) = ai(m:igh,j) - fr * orti(m:igh) - fi * ortr(m:igh)

    end do
!
!  Form (I-(U*Ut)/h) * A * (I-(U*Ut)/h)
!
    do i = 1, igh

      fr = real(0.0E+00,kind=8)
      fi = real(0.0E+00,kind=8)

      do jj = m, igh
        j = mp - jj
        fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
        fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
      end do

      fr = fr / h
      fi = fi / h

      ar(i,m:igh) = ar(i,m:igh) - fr * ortr(m:igh) - fi * orti(m:igh)
      ai(i,m:igh) = ai(i,m:igh) + fr * orti(m:igh) - fi * ortr(m:igh)

    end do

    ortr(m) = scale * ortr(m)
    orti(m) = scale * orti(m)
    ar(m,m-1) = - g * ar(m,m-1)
    ai(m,m-1) = - g * ai(m,m-1)

  end do

  return
end
subroutine csroot ( xr, xi, yr, yi )
!
!*******************************************************************************
!
!! CSROOT computes the complex square root of a complex quantity.
!
!
!  Discussion:
!
!    The branch of the square function is chosen so that 
!      YR >= real(0.0E+00,kind=8) 
!    and 
!      sign ( YI ) == sign ( XI )
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, real(kind=8) XR, XI, the real(kind=8) and imaginary parts of the quantity
!    whose square root is desired.
!
!    Output, real(kind=8) YR, YI, the real(kind=8) and imaginary parts of the square root.
!
  implicit none
!
  real(kind=8) pythag
  real(kind=8) s
  real(kind=8) ti
  real(kind=8) tr
  real(kind=8) xi
  real(kind=8) xr
  real(kind=8) yi
  real(kind=8) yr
!
  tr = xr
  ti = xi
  s = sqrt ( real(0.5E+00,kind=8) * ( pythag ( tr, ti ) + abs ( tr ) ) )

  if ( tr >= real(0.0E+00,kind=8) ) yr = s
  if ( ti < real(0.0E+00,kind=8) ) s = -s
  if ( tr <= real(0.0E+00,kind=8) ) yi = s

  if ( tr < real(0.0E+00,kind=8) ) then
    yr = real(0.5E+00,kind=8) * ( ti / yi )
  else if ( tr > real(0.0E+00,kind=8) ) then
    yi = real(0.5E+00,kind=8) * ( ti / yr )
  end if

  return
end
subroutine elmbak ( nm, low, igh, a, ind, m, z )
!
!*******************************************************************************
!
!! ELMBAK determines eigenvectors by undoing the ELMHES transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a real(kind=8) general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by ELMHES.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer LOW, IGH, integers determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1 and 
!    IGH equal to the order of the matrix.
!
!    Input, real(kind=8) A(NM,IGH), the multipliers which were used in the
!    reduction by ELMHES in its lower triangle below the subdiagonal.
!
!    Input, integer IND(IGH), information on the rows and columns
!    interchanged in the reduction by ELMHES.
!
!    Input, integer M, the number of columns of Z to be back transformed.
!
!    Input/output, real(kind=8) Z(NM,M).  On input, the real(kind=8) and imaginary parts 
!    of the eigenvectors to be back transformed.  On output, the real(kind=8) and 
!    imaginary parts of the transformed eigenvectors.
!
  implicit none
!
  integer igh
  integer m
  integer nm
!
  real(kind=8) a(nm,igh)
  integer i
  integer ind(igh)
  integer j
  integer la
  integer low
  integer mm
  integer mp
  real(kind=8) x
  real(kind=8) z(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  la = igh - 1

  if ( la < low + 1 ) then
    return
  end if

  do mm = low + 1, la

     mp = low + igh - mm

     do i = mp+1, igh

       x = a(i,mp-1)
       if ( x /= real(0.0E+00,kind=8) ) then
         do j = 1, m
           z(i,j) = z(i,j) + x * z(mp,j)
         end do
       end if

     end do

     i = ind(mp)

     if ( i /= mp ) then

       do j = 1, m
         call r_swap ( z(i,j), z(mp,j) )
       end do

     end if

  end do

  return
end
subroutine elmhes ( nm, n, low, igh, a, ind )
!
!*******************************************************************************
!
!! ELMHES transforms a real(kind=8) general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!    Given a real(kind=8) general matrix, this subroutine reduces a submatrix 
!    situated in rows and columns LOW through IGH to upper Hessenberg 
!    form by stabilized elementary similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    Martin and Wilkinson,
!    ELMHES,
!    Numerische Mathematik,
!    Volume 12, pages 349-368, 1968.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the array A.
!    NM must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine  
!    BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input/output, real(kind=8) A(NM,N).  On input, the matrix to be reduced.
!    On output, the Hessenberg matrix.  The multipliers
!    which were used in the reduction are stored in the
!    remaining triangle under the Hessenberg matrix.
!
!    Output, integer IND(N), contains information on the rows and columns
!    interchanged in the reduction.  Only elements LOW through IGH are used.
!
  implicit none
!
  integer igh
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  integer i
  integer ind(igh)
  integer j
  integer la
  integer low
  integer m
  real(kind=8) x
  real(kind=8) y
!
  la = igh - 1

  do m = low + 1, la

    x = real(0.0E+00,kind=8)
    i = m

    do j = m, igh
      if ( abs ( a(j,m-1) ) > abs ( x ) ) then
        x = a(j,m-1)
        i = j
      end if
    end do

    ind(m) = i
!
!  Interchange rows and columns of the matrix.
!
    if ( i /= m ) then

      do j = m-1, n
        call r_swap ( a(i,j), a(m,j) )
      end do

      do j = 1, igh
        call r_swap ( a(j,i), a(j,m) )
      end do

    end if

    if ( x /= real(0.0E+00,kind=8) ) then

      do i = m+1, igh

        y = a(i,m-1)

        if ( y /= real(0.0E+00,kind=8) ) then

          y = y / x
          a(i,m-1) = y

          do j = m, n
            a(i,j) = a(i,j) - y * a(m,j)
          end do

          a(1:igh,m) = a(1:igh,m) + y * a(1:igh,i)

        end if

      end do

    end if

  end do

  return
end
subroutine eltran ( nm, n, low, igh, a, ind, z )
!
!*******************************************************************************
!
!! ELTRAN accumulates similarity transformations used by ELMHES.
!
!
!  Discussion:
!
!    This subroutine accumulates the stabilized elementary
!    similarity transformations used in the reduction of a
!    real(kind=8) general matrix to upper Hessenberg form by ELMHES.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    Peters and WIlkinson,
!    ELMTRANS,
!    Numerische Mathematik,
!    Volume 16, pages 181-204, 1970.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine  
!    BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input, real(kind=8) A(NM,IGH), the multipliers which were used in the
!    reduction by ELMHES in its lower triangle below the subdiagonal.
!
!    Input, integer IND(IGH), information on the rows and columns
!    interchanged in the reduction by ELMHES.
!
!    Output, real(kind=8) Z(NM,N), the transformation matrix produced in the
!    reduction by ELMHES.
!
  implicit none
!
  integer igh
  integer n
  integer nm
!
  real(kind=8) a(nm,igh)
  integer i
  integer ind(igh)
  integer kl
  integer low
  integer mm
  integer mp
  real(kind=8) z(nm,n)
!
!  Initialize Z to the identity matrix.
!
  call rmat_ident ( nm, n, z )

  kl = igh - low - 1

  if ( kl < 1 ) then
    return
  end if

  do mm = 1, kl

     mp = igh - mm

     do i = mp+1, igh
       z(i,mp) = a(i,mp-1)
     end do

     i = ind(mp)

     if ( i /= mp ) then

       z(mp,mp:igh) = z(i,mp:igh)

       z(i,mp) = real(1.0E+00,kind=8)
       z(i,mp+1:igh) = real(0.0E+00,kind=8)

     end if

  end do

  return
end
subroutine figi ( nm, n, t, d, e, e2, ierr )
!
!*******************************************************************************
!
!! FIGI transforms a real(kind=8) nonsymmetric tridiagonal matrix to symmetric form.
!
!
!  Discussion:
!
!    Given a nonsymmetric tridiagonal matrix such that the products
!    of corresponding pairs of off-diagonal elements are all
!    non-negative, this subroutine reduces it to a symmetric
!    tridiagonal matrix with the same eigenvalues.  If, further,
!    a zero product only occurs when both factors are zero,
!    the reduced matrix is similar to the original matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) T(NM,3) contains the input matrix.  Its subdiagonal is
!    stored in the last N-1 positions of the first column, its diagonal in 
!    the N positions of the second column, and its superdiagonal in the 
!    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.
!
!    Output, real(kind=8) D(N), the diagonal elements of the symmetric matrix.
!
!    Output, real(kind=8) E(N), contains the subdiagonal elements of the symmetric
!    matrix in E(2:N).  E(1) is not set.
!
!    Output, real(kind=8) E2(N), the squares of the corresponding elements of E.
!    E2 may coincide with E if the squares are not needed.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    N+I, if T(I,1) * T(I-1,3) is negative,
!    -(3*N+I), if T(I,1) * T(I-1,3) is zero with one factor non-zero.  In 
!      this case, the eigenvectors of the symmetric matrix are not simply 
!      related to those of T and should not be sought.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  integer i
  integer ierr
  real(kind=8) t(nm,3)
!
  ierr = 0

  do i = 1, n

    if ( i >= 1 ) then

      e2(i) = t(i,1) * t(i-1,3)

      if ( e2(i) < real(0.0E+00,kind=8) ) then

        ierr = n + i
        return

      else if ( e2(i) == real(0.0E+00,kind=8) ) then

        if ( t(i,1) /= real(0.0E+00,kind=8) .or. t(i-1,3) /= real(0.0E+00,kind=8) ) then
          ierr = - 3 * n - i
          return
        end if

        e(i) = real(0.0E+00,kind=8)

      else

        e(i) = sqrt ( e2(i) )

      end if

    end if

    d(i) = t(i,2)

  end do

  return
end
subroutine figi2 ( nm, n, t, d, e, z, ierr )
!
!*******************************************************************************
!
!! FIGI2 transforms a real(kind=8) nonsymmetric tridiagonal matrix to symmetric form.
!
!
!  Discussion:
!
!    Given a nonsymmetric tridiagonal matrix such that the products
!    of corresponding pairs of off-diagonal elements are all
!    non-negative, and zero only when both factors are zero, this
!    subroutine reduces it to a symmetric tridiagonal matrix
!    using and accumulating diagonal similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) T(NM,3) contains the input matrix.  Its subdiagonal is
!    stored in the last N-1 positions of the first column, its diagonal in 
!    the N positions of the second column, and its superdiagonal in the 
!    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.
!
!    Output, real(kind=8) D(N), the diagonal elements of the symmetric matrix.
!
!    Output, real(kind=8) E(N), contains the subdiagonal elements of the symmetric
!    matrix in E(2:N).  E(1) is not set.
!
!    Output, real(kind=8) Z(NM,N), contains the transformation matrix produced in
!    the reduction.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    N+I, if T(I,1) * T(I-1,3) is negative,
!    2*N+I, if T(I,1) * T(I-1,3) is zero with one factor non-zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) h
  integer i
  integer ierr
  real(kind=8) t(nm,3)
  real(kind=8) z(nm,n)
!
  ierr = 0

  do i = 1, n

    z(i,1:n) = real(0.0E+00,kind=8)

    if ( i == 1 ) then

      z(i,i) = real(1.0E+00,kind=8)

    else

      h = t(i,1) * t(i-1,3)

      if ( h < real(0.0E+00,kind=8) ) then

        ierr = n + i
        return

      else if ( h == 0 ) then

        if ( t(i,1) /= real(0.0E+00,kind=8) .or. t(i-1,3) /= real(0.0E+00,kind=8) ) then
          ierr = 2 * n + i
          return
        end if

        e(i) = real(0.0E+00,kind=8)
        z(i,i) = real(1.0E+00,kind=8)

      else

        e(i) = sqrt ( h )
        z(i,i) = z(i-1,i-1) * e(i) / t(i-1,3)

      end if

    end if

    d(i) = t(i,2)

  end do

  return
end
subroutine hqr ( nm, n, low, igh, h, wr, wi, ierr )
!
!*******************************************************************************
!
!! HQR computes all eigenvalues of a real(kind=8) upper Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a real(kind=8)
!    upper Hessenberg matrix by the QR method.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    Martin, Peters, and Wilkinson,
!    HQR,
!    Numerische Mathematik,
!    Volume 14, pages 219-231, 1970.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of H, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, two integers determined by the routine
!    BALANC.  If BALANC is not used, set LOW=1, IGH=N.
!
!    Input/output, real(kind=8) H(NM,N), the N by N upper Hessenberg matrix.
!    Information about the transformations used in the reduction to 
!    Hessenberg form by ELMHES or ORTHES, if performed, is stored
!    in the remaining triangle under the Hessenberg matrix.
!    On output, the information in H has been destroyed.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts of the
!    eigenvalues.  The eigenvalues are unordered, except that complex
!    conjugate pairs of values appear consecutively, with the eigenvalue
!    having positive imaginary part listed first.  If an error exit
!    occurred, then the eigenvalues should be correct for indices
!    IERR+1 through N.
!
!    Output, integer IERR, error flag.
!    0, no error.
!    J, the limit of 30*N iterations was reached while searching for
!      the J-th eigenvalue.
!
  implicit none
!
  integer n
  integer nm
!
  integer en
  integer enm2
  real(kind=8) h(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer na
  real(kind=8) norm
  logical notlas
  real(kind=8) p
  real(kind=8) q
  real(kind=8) r
  real(kind=8) s
  real(kind=8) t
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) w
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) x
  real(kind=8) y
  real(kind=8) zz
!
  ierr = 0
  norm = real(0.0E+00,kind=8)
  k = 1
!
!  Store roots isolated by BALANC and compute matrix norm.
!
  do i = 1, n

    do j = k, n
      norm = norm + abs ( h(i,j) )
    end do

    k = i
    if ( i < low .or. i > igh ) then
      wr(i) = h(i,i)
      wi(i) = real(0.0E+00,kind=8)
    end if

  end do

  en = igh
  t = real(0.0E+00,kind=8)
  itn = 30 * n
!
!  Search for next eigenvalues.
!
60 continue

  if ( en < low ) then
    return
  end if

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for a single small sub-diagonal element.
!
70 continue

  do ll = low, en
    l = en + low - ll
    if ( l == low ) then
      exit
    end if
    s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
    if ( s == real(0.0E+00,kind=8) ) then
      s = norm
    end if
    tst1 = s
    tst2 = tst1 + abs ( h(l,l-1) )
    if ( tst2 == tst1 ) then
      exit
    end if
  end do
!
!  Form shift.
!
  x = h(en,en)
  if ( l == en ) then
    go to 270
  end if

  y = h(na,na)
  w = h(en,na) * h(na,en)

  if ( l == na ) then
    go to 280
  end if

  if ( itn == 0 ) then
    ierr = en
    return
  end if
!
!  Form an exceptional shift.
!
  if ( its == 10 .or. its == 20 ) then

    t = t + x

    do i = low, en
      h(i,i) = h(i,i) - x
    end do

    s = abs ( h(en,na) ) + abs ( h(na,enm2) )
    x = 0.75E+00 * s
    y = x
    w = -0.4375E+00 * s * s

  end if

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  do mm = l, enm2

    m = enm2 + l - mm
    zz = h(m,m)
    r = x - zz
    s = y - zz
    p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
    q = h(m+1,m+1) - zz - r - s
    r = h(m+2,m+1)
    s = abs ( p ) + abs ( q ) + abs ( r )
    p = p / s
    q = q / s
    r = r / s

    if ( m == l ) then
      exit
    end if

    tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
    tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )

    if ( tst2 == tst1 ) then
      exit
    end if

  end do

  do i = m+2, en
    h(i,i-2) = real(0.0E+00,kind=8)
    if ( i /= m+2 ) then
      h(i,i-3) = real(0.0E+00,kind=8)
    end if
  end do
!
!  Double QR step involving rows l to EN and columns M to EN.
!
  do k = m, na

    notlas = k /= na

    if ( k /= m ) then

      p = h(k,k-1)
      q = h(k+1,k-1)

      if ( notlas ) then
        r = h(k+2,k-1)
      else
        r = real(0.0E+00,kind=8)
      end if

      x = abs ( p ) + abs ( q ) + abs ( r )

      if ( x == real(0.0E+00,kind=8) ) then
        cycle
      end if

      p = p / x
      q = q / x
      r = r / x

    end if

    s = sign ( sqrt ( p**2 + q**2 + r**2 ), p )

    if ( k /= m ) then
      h(k,k-1) = - s * x
    else if ( l /= m ) then
      h(k,k-1) = - h(k,k-1)
    end if

    p = p + s
    x = p / s
    y = q / s
    zz = r / s
    q = q / p
    r = r / p

    if ( .not. notlas ) then
!
!  Row modification.
!
      do j = k, n
        p = h(k,j) + q * h(k+1,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
      end do

      j = min ( en, k+3 )
!
!  Column modification.
!
      do i = 1, j
        p = x * h(i,k) + y * h(i,k+1)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
      end do

    else
!
!  Row modification.
!
      do j = k, n
        p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
        h(k+2,j) = h(k+2,j) - p * zz
      end do

      j = min ( en, k+3 )
!
!  Column modification.
!
      do i = 1, j
        p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
        h(i,k+2) = h(i,k+2) - p * r
      end do

    end if

  end do

  go to 70
!
!  One root found.
!
270 continue

  wr(en) = x + t
  wi(en) = real(0.0E+00,kind=8)
  en = na
  go to 60
!
!  Two roots found.
!
280 continue

  p = ( y - x ) / real(2.0E+00,kind=8)
  q = p * p + w
  zz = sqrt ( abs ( q ) )
  x = x + t
!
!  real(kind=8) root, or complex pair.
!
  if ( q >= real(0.0E+00,kind=8) ) then

    zz = p + sign ( zz, p )
    wr(na) = x + zz
    if ( zz == real(0.0E+00,kind=8) ) then
      wr(en) = wr(na)
    else
      wr(en) = x - w / zz
    end if
    wi(na) = real(0.0E+00,kind=8)
    wi(en) = real(0.0E+00,kind=8)

  else

    wr(na) = x + p
    wr(en) = x + p
    wi(na) = zz
    wi(en) = -zz

  end if

  en = enm2
  go to 60
end
subroutine hqr2 ( nm, n, low, igh, h, wr, wi, z, ierr )
!
!*******************************************************************************
!
!! HQR2 computes eigenvalues and eigenvectors of a real(kind=8) upper Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors
!    of a real(kind=8) upper Hessenberg matrix by the qr method.  the
!    eigenvectors of a real(kind=8) general matrix can also be found
!    if ELMHES and ELTRAN or ORTHES and ORTRAN have
!    been used to reduce this general matrix to Hessenberg form
!    and to accumulate the similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of H, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, determined by the balancing routine BALANC.  
!    If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input/output, real(kind=8) H(NM,N), the N by N upper Hessenberg matrix.
!    On output, the information in H has been destroyed.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts of the
!    eigenvalues.  The eigenvalues are unordered, except that complex
!    conjugate pairs of values appear consecutively, with the eigenvalue
!    having positive imaginary part listed first.  If an error exit
!    occurred, then the eigenvalues should be correct for indices
!    IERR+1 through N.
!
!    Input/output, real(kind=8) Z(NM,N).  On input, the transformation matrix 
!    produced by ELTRAN after the reduction by ELMHES, or by ORTRAN after the
!    reduction by ORTHES, if performed.  If the eigenvectors of the Hessenberg 
!    matrix are desired, Z must contain the identity matrix.  On output,
!    Z contains the real(kind=8) and imaginary parts of the eigenvectors.
!    If the I-th eigenvalue is real(kind=8), the I-th column of Z contains its 
!    eigenvector.  If the I-th eigenvalue is complex with positive imaginary 
!    part, the I-th and (I+1)-th columns of Z contain the real(kind=8) and imaginary 
!    parts of its eigenvector.  The eigenvectors are unnormalized.  If an
!    error exit is made, none of the eigenvectors has been found.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the limit of 30*N iterations is exhausted while the J-th 
!      eigenvalue is being sought.
!
  implicit none
!
  integer nm
  integer n
!
  integer en
  integer enm2
  real(kind=8) h(nm,n)
  integer i
  integer ierr
  integer igh
  integer ii
  integer itn
  integer its
  integer j
  integer jj
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer na
  integer nn
  real(kind=8) norm
  logical notlas
  real(kind=8) p
  real(kind=8) q
  real(kind=8) r
  real(kind=8) ra
  real(kind=8) s
  real(kind=8) sa
  real(kind=8) t
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) vi
  real(kind=8) vr
  real(kind=8) w
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) x
  real(kind=8) y
  real(kind=8) z(nm,n)
  real(kind=8) zz
!
  ierr = 0
  norm = real(0.0E+00,kind=8)
  k = 1
!
!  Store roots isolated by BALANC and compute the matrix norm.
!
  do i = 1, n

    do j = k, n
      norm = norm + abs ( h(i,j) )
    end do

    k = i
    if ( i < low .or. i > igh ) then
      wr(i) = h(i,i)
      wi(i) = real(0.0E+00,kind=8)
    end if

  end do

  en = igh
  t = real(0.0E+00,kind=8)
  itn = 30 * n
!
!  Search for next eigenvalues.
!
60 continue

  if ( en < low ) then
    go to 340
  end if

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for single small sub-diagonal element.
!
70 continue

  do ll = low, en

    l = en + low - ll

    if ( l == low ) then
      exit
    end if

    s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
    if ( s == real(0.0E+00,kind=8) ) then
      s = norm
    end if

    tst1 = s
    tst2 = tst1 + abs ( h(l,l-1) )

    if ( tst2 == tst1 ) then
      exit
    end if

  end do
!
!  Form shift.
!
  x = h(en,en)
  if ( l == en ) then
    go to 270
  end if

  y = h(na,na)
  w = h(en,na) * h(na,en)

  if ( l == na ) then
    go to 280
  end if

  if ( itn == 0 ) then
    ierr = en
    return
  end if
!
!  Form exceptional shift.
!
  if ( its == 10 .or. its == 20 ) then

    t = t + x

    do i = low, en
      h(i,i) = h(i,i) - x
    end do

    s = abs ( h(en,na) ) + abs ( h(na,enm2) )
    x = 0.75E+00 * s
    y = x
    w = -0.4375E+00 * s * s

  end if

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  do mm = l, enm2
     m = enm2 + l - mm
     zz = h(m,m)
     r = x - zz
     s = y - zz
     p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
     q = h(m+1,m+1) - zz - r - s
     r = h(m+2,m+1)
     s = abs ( p ) + abs ( q ) + abs ( r )
     p = p / s
     q = q / s
     r = r / s
     if ( m == l ) then
       exit
     end if

     tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
     tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) ) 
     if ( tst2 == tst1 ) then
       exit
     end if

  end do

  do i = m+2, en
    h(i,i-2) = real(0.0E+00,kind=8)
    if ( i /= m+2 ) then
      h(i,i-3) = real(0.0E+00,kind=8)
    end if
  end do
!
!  Double QR step involving rows L to EN and columns M to EN.
!
  do k = m, na

     notlas = k /= na

     if ( k /= m ) then

       p = h(k,k-1)
       q = h(k+1,k-1)
       r = real(0.0E+00,kind=8)
       if ( notlas ) then
         r = h(k+2,k-1)
       end if

       x = abs ( p ) + abs ( q ) + abs ( r )
       if ( x == real(0.0E+00,kind=8) ) then
         cycle
       end if

       p = p / x
       q = q / x
       r = r / x

     end if

     s = sign ( sqrt ( p**2 + q**2 + r**2 ), p )

     if ( k /= m ) then
       h(k,k-1) = - s * x
     else if ( l /= m ) then
       h(k,k-1) = -h(k,k-1)
     end if

     p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p
     if ( notlas ) go to 225
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
     end do
!
!  Accumulate transformations.
!
     do i = low, igh
       p = x * z(i,k) + y * z(i,k+1)
       z(i,k) = z(i,k) - p
       z(i,k+1) = z(i,k+1) - p * q
     end do

     go to 255

225  continue
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
       h(k+2,j) = h(k+2,j) - p * zz
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
       h(i,k+2) = h(i,k+2) - p * r
     end do
!
!  Accumulate transformations.
!
     do i = low, igh
        p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
        z(i,k) = z(i,k) - p
        z(i,k+1) = z(i,k+1) - p * q
        z(i,k+2) = z(i,k+2) - p * r
     end do

255 continue

260 continue

  end do

  go to 70
!
!  One root found.
!
270 continue

  h(en,en) = x + t
  wr(en) = h(en,en)
  wi(en) = real(0.0E+00,kind=8)
  en = na
  go to 60
!
!  Two roots found.
!
280 continue

  p = ( y - x ) / real(2.0E+00,kind=8)
  q = p * p + w
  zz = sqrt ( abs ( q ) )
  h(en,en) = x + t
  x = h(en,en)
  h(na,na) = y + t

  if ( q < real(0.0E+00,kind=8) ) go to 320
!
!  real(kind=8) pair.
!
  zz = p + sign ( zz, p )
  wr(na) = x + zz
  wr(en) = wr(na)

  if ( zz /= real(0.0E+00,kind=8) ) then
    wr(en) = x - w / zz
  end if

  wi(na) = real(0.0E+00,kind=8)
  wi(en) = real(0.0E+00,kind=8)
  x = h(en,na)
  s = abs ( x ) + abs ( zz )
  p = x / s
  q = zz / s
  r = sqrt ( p**2 + q**2 )
  p = p / r
  q = q / r
!
!  Row modification.
!
  do j = na, n
    zz = h(na,j)
    h(na,j) = q * zz + p * h(en,j)
    h(en,j) = q * h(en,j) - p * zz
  end do
!
!  Column modification.
!
  do i = 1, en
    zz = h(i,na)
    h(i,na) = q * zz + p * h(i,en)
    h(i,en) = q * h(i,en) - p * zz
  end do
!
!  Accumulate transformations.
!
  do i = low, igh
    zz = z(i,na)
    z(i,na) = q * zz + p * z(i,en)
    z(i,en) = q * z(i,en) - p * zz
  end do

  go to 330
!
!  Complex pair
!
320 continue

  wr(na) = x + p
  wr(en) = x + p
  wi(na) = zz
  wi(en) = -zz

330 continue

  en = enm2
  go to 60
!
!  All roots found.  
!  Backsubstitute to find vectors of upper triangular form.
!
340 continue

  if ( norm == real(0.0E+00,kind=8) ) then
    return
  end if

  do nn = 1, n

     en = n + 1 - nn
     p = wr(en)
     q = wi(en)
     na = en - 1
     if ( q ) 710, 600, 800
!
!  real(kind=8) vector
!
600  continue

     m = en
     h(en,en) = real(1.0E+00,kind=8)

     if ( na == 0 ) go to 800

     do ii = 1, na

        i = en - ii
        w = h(i,i) - p
        r = dot_product ( h(i,m:en), h(m:en,en) )

        if ( wi(i) < real(0.0E+00,kind=8) ) then
          zz = w
          s = r
          go to 700
        end if

        m = i
        if ( wi(i) /= real(0.0E+00,kind=8) ) go to 640
        t = w

        if ( t == real(0.0E+00,kind=8) ) then

          tst1 = norm
          t = tst1

          do
            t = 0.01E+00 * t
            tst2 = norm + t
            if ( tst2 <= tst1 ) then
              exit
            end if
          end do

        end if

        h(i,en) = -r / t
        go to 680
!
!  Solve real(kind=8) equations.
!
640     continue

        x = h(i,i+1)
        y = h(i+1,i)
        q = ( wr(i) - p ) * ( wr(i) - p) + wi(i) * wi(i)
        t = ( x * s - zz * r ) / q
        h(i,en) = t

        if ( abs ( x ) > abs ( zz ) ) then
          h(i+1,en) = (-r - w * t) / x
        else
          h(i+1,en) = (-s - y * t) / zz
        end if
!
!  Overflow control.
!
680     continue

        t = abs ( h(i,en) )

        if ( t /= real(0.0E+00,kind=8) ) then

          tst1 = t
          tst2 = tst1 + real(1.0E+00,kind=8) / tst1

          if ( tst2 <= tst1 ) then
            h(i:en,en) = h(i:en,en) / t
          end if

        end if

700   continue

    end do
!
!  End real(kind=8) vector
!
     go to 800
!
!  Complex vector
!
710  continue

     m = na
!
!  Last vector component chosen imaginary, so that the eigenvector 
!  matrix is triangular.
!
     if ( abs ( h(en,na) ) > abs ( h(na,en) ) ) then

       h(na,na) = q / h(en,na)
       h(na,en) = -(h(en,en) - p) / h(en,na)

     else

       call cdiv ( real(0.0E+00,kind=8), -h(na,en), h(na,na)-p, q, h(na,na), h(na,en) )

     end if

     h(en,na) = real(0.0E+00,kind=8)
     h(en,en) = real(1.0E+00,kind=8)
     enm2 = na - 1

     do ii = 1, enm2

        i = na - ii
        w = h(i,i) - p
        ra = dot_product ( h(i,m:en), h(m:en,na) )
        sa = dot_product ( h(i,m:en), h(m:en,en) )

        if ( wi(i) < real(0.0E+00,kind=8) ) then
          zz = w
          r = ra
          s = sa
        end if

         m = i

        if ( wi(i) == real(0.0E+00,kind=8) ) then
          call cdiv ( -ra, -sa, w, q, h(i,na), h(i,en) )
          go to 790
        end if
!
!  Solve complex equations.
!
        x = h(i,i+1)
        y = h(i+1,i)
        vr = ( wr(i) - p ) * ( wr(i) - p ) + wi(i) * wi(i) - q * q
        vi = ( wr(i) - p ) * real(2.0E+00,kind=8) * q

        if ( vr == real(0.0E+00,kind=8) .and. vi == real(0.0E+00,kind=8) ) then

          tst1 = norm * ( abs ( w ) + abs ( q ) + abs ( x ) &
            + abs ( y ) + abs ( zz ) )
          vr = tst1

          do
            vr = 0.01E+00 * vr
            tst2 = tst1 + vr
            if ( tst2 <= tst1 ) then
              exit
            end if
          end do

        end if

        call cdiv ( x*r-zz*ra+q*sa, x*s-zz*sa-q*ra, vr, vi, h(i,na), h(i,en) )

        if ( abs ( x ) > abs ( zz ) + abs ( q ) ) then
          h(i+1,na) = ( -ra - w * h(i,na) + q * h(i,en) ) / x
          h(i+1,en) = ( -sa - w * h(i,en) - q * h(i,na) ) / x
        else
          call cdiv ( -r-y*h(i,na), -s-y*h(i,en), zz, q, h(i+1,na), h(i+1,en) )
        end if
!
!  Overflow control.
!
790     continue

        t = max ( abs ( h(i,na) ), abs ( h(i,en) ) )

        if ( t /= real(0.0E+00,kind=8) ) then
          tst1 = t
          tst2 = tst1 + real(1.0E+00,kind=8) / tst1
          if ( tst2 <= tst1 ) then
            h(i:en,na) = h(i:en,na) / t
            h(i:en,en) = h(i:en,en) / t
          end if
        end if

795     continue

      end do
!
!  End complex vector.
!
800 continue

  end do
!
!  End back substitution.
!
!  Vectors of isolated roots.
!
  do i = 1, n

    if ( i < low .or. i > igh ) then
      z(i,i:n) = h(i,i:n)
    end if

  end do
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  do jj = low, n

     j = n + low - jj
     m = min ( j, igh )

     do i = low, igh
       z(i,j) = dot_product ( z(i,low:m), h(low:m,j) )
     end do

  end do

  return
end
subroutine htrib3 ( nm, n, a, tau, m, zr, zi )
!
!*******************************************************************************
!
!! HTRIB3 determines eigenvectors by undoing the HTRID3 transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a complex hermitian
!    matrix by back transforming those of the corresponding
!    real(kind=8) symmetric tridiagonal matrix determined by HTRID3.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, is the order of the matrix.
!
!    Input, real(kind=8) A(NM,N), contains information about the unitary 
!    transformations used in the reduction by HTRID3.
!
!    Input, real(kind=8) TAU(2,N), contains further information about the 
!    transformations.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) ZR(NM,M), ZI(NM,M).  On input, ZR contains the
!    eigenvectors to be back transformed.  On output, ZR and ZI contain
!    the real(kind=8) and imaginary parts of the transformed eigenvectors.
!
  implicit none
!
  integer m
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) h
  integer i
  integer j
  integer k
  integer l
  real(kind=8) s
  real(kind=8) si
  real(kind=8) tau(2,n)
  real(kind=8) zi(nm,m)
  real(kind=8) zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if
!
!  Transform the eigenvectors of the real(kind=8) symmetric tridiagonal matrix 
!  to those of the hermitian tridiagonal matrix.
!
  do k = 1, n
    do j = 1, m
      zi(k,j) = -zr(k,j) * tau(2,k)
      zr(k,j) = zr(k,j) * tau(1,k)
    end do
  end do
!
!  Recover and apply the Householder matrices.
!
  do i = 2, n

    l = i - 1
    h = a(i,i)

    if ( h /= real(0.0E+00,kind=8) ) then

      do j = 1, m

        s = real(0.0E+00,kind=8)
        si = real(0.0E+00,kind=8)

        do k = 1, l
          s = s + a(i,k) * zr(k,j) - a(k,i) * zi(k,j)
          si = si + a(i,k) * zi(k,j) + a(k,i) * zr(k,j)
        end do

        s = ( s / h ) / h
        si = ( si / h ) / h

        zr(1:l,j) = zr(1:l,j) - s * a(i,1:l) - si * a(1:l,i)
        zi(1:l,j) = zi(1:l,j) - si * a(i,1:l) + s * a(1:l,i)

      end do

    end if

  end do

  return
end
subroutine htribk ( nm, n, ar, ai, tau, m, zr, zi )
!
!*******************************************************************************
!
!! HTRIBK determines eigenvectors by undoing the HTRIDI transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a complex hermitian
!    matrix by back transforming those of the corresponding
!    real(kind=8) symmetric tridiagonal matrix determined by HTRIDI.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) AR(NM,N), AI(NM,N), contain information about the unitary 
!    transformations used in the reduction by HTRIDI in their
!    full lower triangles, except for the diagonal of AR.
!
!    Input, real(kind=8) TAU(2,N), contains further information about the 
!    transformations.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) ZR(NM,M), ZI(NM,M).  On input, ZR contains the
!    eigenvectors to be back transformed.  On output, ZR and ZI contain
!    the real(kind=8) and imaginary parts of the transformed eigenvectors.
!
  implicit none
!
  integer m
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  real(kind=8) h
  integer i
  integer j
  integer k
  integer l
  real(kind=8) s
  real(kind=8) si
  real(kind=8) tau(2,n)
  real(kind=8) zi(nm,m)
  real(kind=8) zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if
!
!  Transform the eigenvectors of the real(kind=8) symmetric tridiagonal matrix to 
!  those of the hermitian tridiagonal matrix.
!
  do k = 1, n
    do j = 1, m
      zi(k,j) = -zr(k,j) * tau(2,k)
      zr(k,j) = zr(k,j) * tau(1,k)
    end do
  end do
!
!  Recover and apply the Householder matrices.
!
  do i = 2, n

    l = i - 1
    h = ai(i,i)

    if ( h /= real(0.0E+00,kind=8) ) then

      do j = 1, m

        s = real(0.0E+00,kind=8)
        si = real(0.0E+00,kind=8)
        do k = 1, l
          s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
          si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
        end do

        s = ( s / h ) / h
        si = ( si / h ) / h

        zr(1:l,j) = zr(1:l,j) - s * ar(i,1:l) - si * ai(i,1:l)
        zi(1:l,j) = zi(1:l,j) - si * ar(i,1:l) + s * ai(i,1:l)

      end do

    end if

  end do

  return
end
subroutine htrid3 ( nm, n, a, d, e, e2, tau )
!
!*******************************************************************************
!
!! HTRID3 tridiagonalizes a complex hermitian packed matrix.
!
!
!  Discussion:
!
!    This subroutine reduces a complex hermitian matrix, stored as
!    a single square array, to a real(kind=8) symmetric tridiagonal matrix
!    using unitary similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) A(NM,N).  On input, the lower triangle of the complex 
!    hermitian input matrix.  The real(kind=8) parts of the matrix elements are stored
!    in the full lower triangle of A, and the imaginary parts are stored in 
!    the transposed positions of the strict upper triangle of A.  No storage 
!    is required for the zero imaginary parts of the diagonal elements.
!    On output, A contains information about the unitary transformations
!    used in the reduction.
!
!    Output, real(kind=8) D(N), the diagonal elements of the the tridiagonal matrix.
!
!    Output, real(kind=8) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is set to zero.
!
!    Output, real(kind=8) E2(N), the squares of the corresponding elements of E.
!    E2 may coincide with E if the squares are not needed.
!
!    Output, real(kind=8) TAU(2,N), contains further information about the 
!    transformations.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) f
  real(kind=8) fi
  real(kind=8) g
  real(kind=8) gi
  real(kind=8) h
  real(kind=8) hh
  integer i
  integer ii
  integer j
  integer k
  integer l
  real(kind=8) pythag
  real(kind=8) scale
  real(kind=8) si
  real(kind=8) tau(2,n)
!
  tau(1,n) = real(1.0E+00,kind=8)
  tau(2,n) = real(0.0E+00,kind=8)

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)

    if ( l < 1 ) then
      e(i) = real(0.0E+00,kind=8)
      e2(i) = real(0.0E+00,kind=8)
      go to 290
    end if
!
!  Scale row.
!
     do k = 1, l
       scale = scale + abs ( a(i,k) ) + abs ( a(k,i) )
     end do

     if ( scale == real(0.0E+00,kind=8) ) then
       tau(1,l) = real(1.0E+00,kind=8)
       tau(2,l) = real(0.0E+00,kind=8)
       e(i) = real(0.0E+00,kind=8)
       e2(i) = real(0.0E+00,kind=8)
       go to 290
     end if
   
      do k = 1, l
        a(i,k) = a(i,k) / scale
        a(k,i) = a(k,i) / scale
        h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i)
     end do

     e2(i) = scale * scale * h
     g = sqrt ( h )
     e(i) = scale * g
     f = pythag ( a(i,l), a(l,i) )
!
!  Form next diagonal element of matrix T.
!
     if ( f /= real(0.0E+00,kind=8) ) then

       tau(1,l) = ( a(l,i) * tau(2,i) - a(i,l) * tau(1,i) ) / f
       si = ( a(i,l) * tau(2,i) + a(l,i) * tau(1,i) ) / f
       h = h + f * g
       g = real(1.0E+00,kind=8) + g / f
       a(i,l) = g * a(i,l)
       a(l,i) = g * a(l,i)

       if ( l == 1 ) go to 270

     else

       tau(1,l) = -tau(1,i)
       si = tau(2,i)
       a(i,l) = g

     end if

     f = real(0.0E+00,kind=8)

     do j = 1, l

        g = real(0.0E+00,kind=8)
        gi = real(0.0E+00,kind=8)
!
!  Form element of A*U.
!
        do k = 1, j-1
          g = g + a(j,k) * a(i,k) + a(k,j) * a(k,i)
          gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k)
        end do

        g = g + a(j,j) * a(i,j)
        gi = gi - a(j,j) * a(j,i)

        do k = j+1, l
          g = g + a(k,j) * a(i,k) - a(j,k) * a(k,i)
          gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k)
        end do
!
!  Form element of P.
!
        e(j) = g / h
        tau(2,j) = gi / h
        f = f + e(j) * a(i,j) - tau(2,j) * a(j,i)

     end do

     hh = f / ( h + h )
!
!  Form reduced A.
!
     do j = 1, l

        f = a(i,j)
        g = e(j) - hh * f
        e(j) = g
        fi = -a(j,i)
        gi = tau(2,j) - hh * fi
        tau(2,j) = -gi
        a(j,j) = a(j,j) - real(2.0E+00,kind=8) * ( f * g + fi * gi )

        do k = 1, j-1
          a(j,k) = a(j,k) - f * e(k) - g * a(i,k) + fi * tau(2,k) + gi * a(k,i)
          a(k,j) = a(k,j) - f * tau(2,k) - g * a(k,i) - fi * e(k) - gi * a(i,k)
        end do

     end do

270  continue

     a(i,1:l) = scale * a(i,1:l)
     a(1:l,i) = scale * a(1:l,i)
     tau(2,l) = -si

290  continue

     d(i) = a(i,i)
     a(i,i) = scale * sqrt ( h )

  end do

  return
end
subroutine htridi ( nm, n, ar, ai, d, e, e2, tau )
!
!*******************************************************************************
!
!! HTRIDI tridiagonalizes a complex hermitian matrix.
!
!
!  Discussion:
!
!    This subroutine reduces a complex hermitian matrix to a real(kind=8) symmetric 
!    tridiagonal matrix using unitary similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) AR(NM,N), AI(NM,N).  On input, the real(kind=8) and imaginary 
!    parts, respectively, of the complex hermitian input matrix.
!    Only the lower triangle of the matrix need be supplied.
!    On output, information about the unitary transformations used in the 
!    reduction in their full lower triangles.  Their strict upper triangles 
!    and the diagonal of AR are unaltered.
!
!    Output, real(kind=8) D(N), the diagonal elements of the the tridiagonal matrix.
!
!    Output, real(kind=8) E(N), the subdiagonal elements of the tridiagonal
!    matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, real(kind=8) E2(N), the squares of the corresponding elements of E.
!    E2 may coincide with E if the squares are not needed.
!
!    Output, real(kind=8) TAU(2,N), contains further information about the 
!    transformations.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) ai(nm,n)
  real(kind=8) ar(nm,n)
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) f
  real(kind=8) fi
  real(kind=8) g
  real(kind=8) gi
  real(kind=8) h
  real(kind=8) hh
  integer i
  integer ii
  integer j
  integer k
  integer l
  real(kind=8) pythag
  real(kind=8) scale
  real(kind=8) si
  real(kind=8) tau(2,n)
!
  tau(1,n) = real(1.0E+00,kind=8)
  tau(2,n) = real(0.0E+00,kind=8)

  do i = 1, n
    d(i) = ar(i,i)
  end do

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)

    if ( l < 1 ) then
      e(i) = real(0.0E+00,kind=8)
      e2(i) = real(0.0E+00,kind=8)
      go to 290
    end if
!
!  Scale row.
!
    do k = 1, l
      scale = scale + abs ( ar(i,k) ) + abs ( ai(i,k) )
    end do

    if ( scale == real(0.0E+00,kind=8) ) then
      tau(1,l) = real(1.0E+00,kind=8)
      tau(2,l) = real(0.0E+00,kind=8)
      e(i) = real(0.0E+00,kind=8)
      e2(i) = real(0.0E+00,kind=8)
      go to 290
    end if

    ar(i,1:l) = ar(i,1:l) / scale
    ai(i,1:l) = ai(i,1:l) / scale

    do k = 1, l
      h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
    end do

    e2(i) = scale * scale * h
    g = sqrt ( h )
    e(i) = scale * g
    f = pythag ( ar(i,l), ai(i,l) )
!
!  Form next diagonal element of matrix T.
!
    if ( f /= real(0.0E+00,kind=8) ) then
      tau(1,l) = ( ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i) ) / f
      si = ( ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i) ) / f
      h = h + f * g
      g = real(1.0E+00,kind=8) + g / f
      ar(i,l) = g * ar(i,l)
      ai(i,l) = g * ai(i,l)
      if ( l == 1 ) go to 270
    else
      tau(1,l) = -tau(1,i)
      si = tau(2,i)
      ar(i,l) = g
    end if

    f = real(0.0E+00,kind=8)

    do j = 1, l

      g = real(0.0E+00,kind=8)
      gi = real(0.0E+00,kind=8)
!
!  Form element of A*U.
!
      do k = 1, j
        g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
        gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
      end do

      do k = j+1, l
        g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
        gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
      end do
!
!  Form element of P.
!
      e(j) = g / h
      tau(2,j) = gi / h
      f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)

    end do

    hh = f / ( h + h )
!
!  Form the reduced A.
!
    do j = 1, l

      f = ar(i,j)
      g = e(j) - hh * f
      e(j) = g
      fi = - ai(i,j)
      gi = tau(2,j) - hh * fi
      tau(2,j) = -gi

      do k = 1, j
        ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) + fi * tau(2,k) &
          + gi * ai(i,k)
        ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) - fi * e(k) &
          - gi * ar(i,k)
      end do

    end do

270 continue

    ar(i,1:l) = scale * ar(i,1:l)
    ai(i,1:l) = scale * ai(i,1:l)
    tau(2,l) = -si

290 continue

    hh = d(i)
    d(i) = ar(i,i)
    ar(i,i) = hh
    ai(i,i) = scale * sqrt ( h )
  
  end do

  return
end
subroutine imtql1 ( n, d, e, ierr )
!
!*******************************************************************************
!
!! IMTQL1 computes all eigenvalues of a symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric
!    tridiagonal matrix by the implicit QL method.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) D(N).  On input, the diagonal elements of the matrix.
!    On output, the eigenvalues in ascending order.  If an error exit is 
!    made, the eigenvalues are correct and ordered for indices 1,2,...IERR-1, 
!    but may not be the smallest eigenvalues.
!
!    Input/output, real(kind=8) E(N).  On input, the subdiagonal elements of the matrix
!    in its last N-1 positions.  E(1) is arbitrary.  On output, E has been
!    overwritten.
!
!    Output, integer IERR, error flag.
!    0, normal return,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  implicit none
!
  integer n
!
  real(kind=8) b
  real(kind=8) c
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) f
  real(kind=8) g
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer m
  integer mml
  real(kind=8) p
  real(kind=8) pythag
  real(kind=8) r
  real(kind=8) s
  real(kind=8) tst1
  real(kind=8) tst2
!
  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do
  e(n) = real(0.0E+00,kind=8)

  do l = 1, n

    j = 0
!
!  Look for a small sub-diagonal element.
!
105 continue

    do m = l, n

      if ( m == n ) then
        exit
      end if

      tst1 = abs ( d(m) ) + abs ( d(m+1) )
      tst2 = tst1 + abs ( e(m) )

      if ( tst2 == tst1 ) then
        exit
      end if

    end do

    p = d(l)
    if ( m == l ) go to 215

    if ( j >= 30 ) then
      ierr = l
      return
    end if

    j = j + 1
!
!  Form shift.
!
    g = ( d(l+1) - p ) / ( real(2.0E+00,kind=8) * e(l) )
    r = pythag ( g, real(1.0E+00,kind=8) )
    g = d(m) - p + e(l) / ( g + sign ( r, g ) )
    s = real(1.0E+00,kind=8)
    c = real(1.0E+00,kind=8)
    p = real(0.0E+00,kind=8)
    mml = m - l

    do ii = 1, mml

      i = m - ii
      f = s * e(i)
      b = c * e(i)
      r = pythag ( f, g )
      e(i+1) = r

      if ( r == real(0.0E+00,kind=8) ) go to 210

      s = f / r
      c = g / r
      g = d(i+1) - p
      r = ( d(i) - g ) * s + real(2.0E+00,kind=8) * c * b
      p = s * r
      d(i+1) = g + p
      g = c * r - b

    end do

    d(l) = d(l) - p
    e(l) = g
    e(m) = real(0.0E+00,kind=8)
    go to 105
!
!  Recover from underflow.
!
210 continue

    d(i+1) = d(i+1) - p
    e(m) = real(0.0E+00,kind=8)
    go to 105
!
!  Order eigenvalues.
!
215 continue

    do ii = 2, l
      i = l + 2 - ii
      if ( p >= d(i-1) ) then
        go to 270
      end if
      d(i) = d(i-1)
    end do

    i = 1

270 continue

    d(i) = p

  end do

  return
end
subroutine imtql2 ( nm, n, d, e, z, ierr )
!
!*******************************************************************************
!
!! IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors
!    of a symmetric tridiagonal matrix by the implicit QL method.
!    The eigenvectors of a full symmetric matrix can also
!    be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) D(N).  On input, the diagonal elements of the input 
!    matrix.  On output, the eigenvalues in ascending order.  If an
!    error exit is made, the eigenvalues are correct but
!    unordered for indices 1,2,...,IERR-1.
!
!    Input/output, real(kind=8) E(N).  On input, the subdiagonal elements of the 
!    input matrix in E(2:N).  E(1) is arbitrary.  On output, E is overwritten.
!
!    Input/output, real(kind=8) Z(NM,N).  On input, the transformation matrix 
!    produced in the reduction by TRED2, if performed.  If the eigenvectors
!    of the tridiagonal matrix are desired, Z must contain the identity matrix.
!    On output, Z contains orthonormal eigenvectors of the symmetric
!    tridiagonal (or full) matrix.  If an error exit is made, Z contains 
!    the eigenvectors associated with the stored eigenvalues.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) b
  real(kind=8) c
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) f
  real(kind=8) g
  integer i
  integer ierr
  integer ii
  integer j
  integer k
  integer l
  integer m
  integer mml
  real(kind=8) p
  real(kind=8) pythag
  real(kind=8) r
  real(kind=8) s
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) z(nm,n)
!
  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do
  e(n) = real(0.0E+00,kind=8)

  do l = 1, n

     j = 0
!
!  Look for a small sub-diagonal element.
!
105  continue

     do m = l, n

       if ( m == n ) then
         exit
       end if

       tst1 = abs ( d(m) ) + abs ( d(m+1) )
       tst2 = tst1 + abs ( e(m) )

       if ( tst2 == tst1 ) then
         exit
       end if

     end do

     p = d(l)

     if ( m == l ) go to 240

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift.
!
     g = ( d(l+1) - p ) / ( real(2.0E+00,kind=8) * e(l) )
     r = pythag ( g, real(1.0E+00,kind=8) )
     g = d(m) - p + e(l) / ( g + sign ( r, g ) )
     s = real(1.0E+00,kind=8)
     c = real(1.0E+00,kind=8)
     p = real(0.0E+00,kind=8)
     mml = m - l

     do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)
        r = pythag ( f, g )
        e(i+1) = r

        if ( r == real(0.0E+00,kind=8) ) go to 210

        s = f / r
        c = g / r
        g = d(i+1) - p
        r = ( d(i) - g ) * s + real(2.0E+00,kind=8) * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
!
!  Form vector.
!
        do k = 1, n
          f = z(k,i+1)
          z(k,i+1) = s * z(k,i) + c * f
          z(k,i) = c * z(k,i) - s * f
        end do

     end do

     d(l) = d(l) - p
     e(l) = g
     e(m) = real(0.0E+00,kind=8)
     go to 105
!
!  Recover from underflow.
!
210  continue

     d(i+1) = d(i+1) - p
     e(m) = real(0.0E+00,kind=8)
     go to 105

240 continue

  end do
!
!  Order eigenvalues and eigenvectors.
!
  do ii = 2, n

     i = ii - 1
     k = i
     p = d(i)

     do j = ii, n
       if ( d(j) < p ) then
         k = j
         p = d(j)
       end if
     end do

     if ( k /= i ) then

       d(k) = d(i)
       d(i) = p

       do j = 1, n
         call r_swap ( z(j,i), z(j,k) )
       end do

    end if

  end do

  return
end
subroutine imtqlv ( n, d, e, e2, w, ind, ierr )
!
!*******************************************************************************
!
!! IMTQLV computes all eigenvalues of a real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric tridiagonal
!    matrix by the implicit QL method and associates with them
!    their corresponding submatrix indices.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) D(N), the diagonal elements of the input matrix.
!
!    Input, real(kind=8) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, real(kind=8) E2(N).  On input, the squares of the corresponding 
!    elements of E.  E2(1) is arbitrary.  On output, elements of E2 
!    corresponding to elements of E regarded as negligible have been 
!    replaced by zero, causing the matrix to split into a direct sum of 
!    submatrices.  E2(1) is also set to zero.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.  If an
!    error exit is made, the eigenvalues are correct and ordered for 
!    indices 1,2,...IERR-1, but may not be the smallest eigenvalues.
!
!    Output, integer IND(N), the submatrix indices associated with the
!    corresponding eigenvalues in W: 1 for eigenvalues belonging to the 
!    first submatrix from the top, 2 for those belonging to the second 
!    submatrix, and so on.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  implicit none
!
  integer n
!
  real(kind=8) b
  real(kind=8) c
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) f
  real(kind=8) g
  integer i
  integer ierr
  integer ii
  integer ind(n)
  integer j
  integer k
  integer l
  integer m
  integer mml
  real(kind=8) p
  real(kind=8) pythag
  real(kind=8) r
  real(kind=8) rv1(n)
  real(kind=8) s
  integer tag
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) w(n)
!
  ierr = 0
  k = 0
  tag = 0
  w(1:n) = d(1:n)
  e2(1) = real(0.0E+00,kind=8)
  rv1(1:n-1) = e(2:n)
  rv1(n) = real(0.0E+00,kind=8)

  do l = 1, n

    j = 0
!
!  Look for a small sub-diagonal element.
!
105 continue

     do m = l, n

       if ( m == n ) then
         exit
       end if

       tst1 = abs ( w(m) ) + abs ( w(m+1) )
       tst2 = tst1 + abs ( rv1(m) )

       if ( tst2 == tst1 ) then
         exit
       end if
!
!  Guard against underflowed element of E2.
!
       if ( e2(m+1) == real(0.0E+00,kind=8) ) go to 125

     end do

120  continue

     if ( m <= k ) go to 130

     if ( m /= n ) e2(m+1) = real(0.0E+00,kind=8)

125  continue

     k = m
     tag = tag + 1

130  continue

     p = w(l)

     if ( m == l ) go to 215

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift.
!
     g = ( w(l+1) - p ) / ( real(2.0E+00,kind=8) * rv1(l) )
     r = pythag ( g, real(1.0E+00,kind=8) )
     g = w(m) - p + rv1(l) / (g + sign ( r, g ) )
     s = real(1.0E+00,kind=8)
     c = real(1.0E+00,kind=8)
     p = real(0.0E+00,kind=8)
     mml = m - l

     do ii = 1, mml
       i = m - ii
       f = s * rv1(i)
       b = c * rv1(i)
       r = pythag ( f, g )
       rv1(i+1) = r

       if ( r == real(0.0E+00,kind=8) ) go to 210

       s = f / r
       c = g / r
       g = w(i+1) - p
       r = ( w(i) - g ) * s + real(2.0E+00,kind=8) * c * b
       p = s * r
       w(i+1) = g + p
       g = c * r - b
     end do

     w(l) = w(l) - p
     rv1(l) = g
     rv1(m) = real(0.0E+00,kind=8)
     go to 105
!
!  Recover from underflow.
!
210  continue

     w(i+1) = w(i+1) - p
     rv1(m) = real(0.0E+00,kind=8)
     go to 105
!
!  Order the eigenvalues.
!
215  continue

     do ii = 2, l
        i = l + 2 - ii
        if ( p >= w(i-1) ) go to 270
        w(i) = w(i-1)
        ind(i) = ind(i-1)
     end do

     i = 1

  270   continue

     w(i) = p
     ind(i) = tag

  end do

  return
end
subroutine invit ( nm, n, a, wr, wi, select, mm, m, z, ierr )
!
!*******************************************************************************
!
!! INVIT computes eigenvectors given eigenvalues, for a real(kind=8) upper Hessenberg matrix.
!
!
!  Discussion:
!
!    This subroutine finds those eigenvectors of a real(kind=8) upper Hessenberg 
!    matrix corresponding to specified eigenvalues, using inverse iteration.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) A(NM,N), the Hessenberg matrix.
!
!    Input/output, real(kind=8) WR(N), WI(N).  On input, the real(kind=8) and imaginary 
!    parts, respectively, of the eigenvalues of the matrix.  The eigenvalues 
!    must be stored in a manner identical to that of subroutine HQR,
!    which recognizes possible splitting of the matrix.  On output,
!    WR may have been altered since close eigenvalues are perturbed
!    slightly in searching for independent eigenvectors.
!
!    Input/output, logical SELECT(N).  On input, specifies the eigenvectors 
!    to be found.  The eigenvector corresponding to the J-th eigenvalue is
!    specified by setting SELECT(J) to TRUE.  On output, SELECT may have been 
!    altered.  If the elements corresponding to a pair of conjugate complex 
!    eigenvalues were each initially set to TRUE, the program resets the 
!    second of the two elements to FALSE.
!
!    Input, integer MM, an upper bound for the number of columns required 
!    to store the eigenvectors to be found.  Note that two columns are 
!    required to store the eigenvector corresponding to a complex eigenvalue.
!
!    Input, integer M, the number of columns actually used to store
!    the eigenvectors.
!
!    Output, real(kind=8) Z(NM,MM), the real(kind=8) and imaginary parts of the eigenvectors.
!    If the next selected eigenvalue is real(kind=8), the next column
!    of Z contains its eigenvector.  If the eigenvalue is complex, the next 
!    two columns of Z contain the real(kind=8) and imaginary parts of its eigenvector.  
!    The eigenvectors are normalized so that the component of largest 
!    magnitude is 1.  Any vector which fails the acceptance test is set to zero.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    -(2*N+1), if more than MM columns of Z are necessary to store the 
!      eigenvectors corresponding to the specified eigenvalues.
!    -K, if the iteration corresponding to the K-th value fails,
!    -(N+K), if both error situations occur.
!
  implicit none
!
  integer nm
  integer n
!
  real(kind=8) a(nm,n)
  real(kind=8) eps3
  real(kind=8) growto
  integer i
  integer ierr
  integer ii
  real(kind=8) ilambd
  integer ip
  integer its
  integer j
  integer k
  integer km1
  integer l
  integer m
  integer mm
  integer mp
  integer n1
  real(kind=8) norm
  real(kind=8) normv
  integer ns
  real(kind=8) pythag
  real(kind=8) rlambd
  real(kind=8) rm1(n,n)
  real(kind=8) rv1(n)
  real(kind=8) rv2(n)
  integer s
  logical select(n)
  real(kind=8) t
  integer uk
  real(kind=8) ukroot
  real(kind=8) w
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) x
  real(kind=8) y
  real(kind=8) z(nm,mm)
!
  ierr = 0
  uk = 0
  s = 1
!
!  The value of IP is:
!
!   0, real(kind=8) eigenvalue;
!   1, first of conjugate complex pair;
!  -1, second of conjugate complex pair.
!
  ip = 0
  n1 = n - 1

  do k = 1, n

     if ( wi(k) /= real(0.0E+00,kind=8) .and. ip >= 0 ) then
       ip = 1
       if ( select(k) .and. select(k+1) ) select(k+1) = .false.
     end if

     if ( .not. select(k) ) go to 960

     if ( wi(k) /= real(0.0E+00,kind=8) ) then
       s = s + 1
     end if

     if ( s > mm ) go to 1000

     if ( uk >= k ) go to 200
!
!  Check for possible splitting.
!
     do uk = k, n
       if ( uk == n ) then
         exit
       end if
       if ( a(uk+1,uk) == real(0.0E+00,kind=8) ) then
         exit
       end if
     end do
!
!  Compute infinity norm of leading UK by UK (Hessenberg) matrix.
!
     norm = real(0.0E+00,kind=8)
     mp = 1

     do i = 1, uk

       x = sum ( abs ( a(i,mp:uk) ) )
       norm = max ( norm, x )
       mp = i

     end do
!
!  EPS3 replaces zero pivot in decomposition and close roots are modified 
!  by EPS3.
!
     if ( norm == real(0.0E+00,kind=8) ) then
       norm = real(1.0E+00,kind=8)
     end if

     eps3 = abs ( norm ) * epsilon ( eps3 )
!
!  GROWTO is the criterion for the growth.
!
     ukroot = uk
     ukroot = sqrt ( ukroot )
     growto = 0.1E+00 / ukroot

200  continue

     rlambd = wr(k)
     ilambd = wi(k)
     if ( k == 1 ) go to 280
     km1 = k - 1
     go to 240
!
!  Perturb eigenvalue if it is close to any previous eigenvalue.
!
220 continue

     rlambd = rlambd + eps3

240  continue

     do ii = 1, km1
       i = k - ii
       if ( select(i) .and. abs ( wr(i) - rlambd ) < eps3 .and. &
          abs ( wi(i) - ilambd ) < eps3 ) then
        go to 220
       end if
     end do

     wr(k) = rlambd
!
!  Perturb conjugate eigenvalue to match.
!
     wr(k+ip) = rlambd
!
!  Form upper Hessenberg A - rlambd*I (transposed) and initial real(kind=8) vector.
!
280  continue

     mp = 1

     do i = 1, uk

        rm1(mp:uk,i) = a(i,mp:uk)

        rm1(i,i) = rm1(i,i) - rlambd
        mp = i
        rv1(i) = eps3

     end do

     its = 0

     if ( ilambd /= real(0.0E+00,kind=8) ) go to 520
!
!  real(kind=8) eigenvalue.
!
!  Triangular decomposition with interchanges, replacing zero pivots by eps3.
!
     do i = 2, uk

        mp = i - 1

        if ( abs ( rm1(mp,i) ) > abs ( rm1(mp,mp) ) ) then

          do j = mp, uk
            call r_swap ( rm1(j,i), rm1(j,mp) )
          end do

        end if

        if ( rm1(mp,mp) == real(0.0E+00,kind=8) ) then
          rm1(mp,mp) = eps3
        end if

        x = rm1(mp,i) / rm1(mp,mp)

        if ( x /= real(0.0E+00,kind=8) ) then
          rm1(i:uk,i) = rm1(i:uk,i) - x * rm1(i:uk,mp)
        end if

      end do

      if ( rm1(uk,uk) == real(0.0E+00,kind=8) ) then
        rm1(uk,uk) = eps3
      end if
!
!  Back substitution for real(kind=8) vector.
!
440   continue

      do ii = 1, uk

        i = uk + 1 - ii
        y = rv1(i)

        do j = i+1, uk
          y = y - rm1(j,i) * rv1(j)
        end do

        rv1(i) = y / rm1(i,i)

     end do

     go to 740
!
!  Complex eigenvalue.
!
!  Triangular decomposition with interchanges,
!  replacing zero pivots by EPS3.  
!  Store imaginary parts in upper triangle starting at (1,3)
!
520  continue

     ns = n - s
     z(1,s-1) = -ilambd
     z(1,s) = real(0.0E+00,kind=8)

     if ( n /= 2 ) then
       rm1(1,3) = -ilambd
       z(1,s-1) = real(0.0E+00,kind=8)
       rm1(1,4:n) = real(0.0E+00,kind=8)
     end if

     do i = 2, uk

        mp = i - 1
        w = rm1(mp,i)

        if ( i < n ) then
          t = rm1(mp,i+1)
        else if ( i == n ) then
          t = z(mp,s-1)
        end if

        x = rm1(mp,mp) * rm1(mp,mp) + t * t

        if ( w * w <= x ) go to 580

        x = rm1(mp,mp) / w
        y = t / w
        rm1(mp,mp) = w

        if ( i < n ) then
          rm1(mp,i+1) = real(0.0E+00,kind=8)
        else if ( i == n ) then
          z(mp,s-1) = real(0.0E+00,kind=8)
        end if

        do j = i, uk

          w = rm1(j,i)
          rm1(j,i) = rm1(j,mp) - x * w
          rm1(j,mp) = w

          if ( j >= n1 ) then
            l = j - ns
            z(i,l) = z(mp,l) - y * w
            z(mp,l) = real(0.0E+00,kind=8)
          else
            rm1(i,j+2) = rm1(mp,j+2) - y * w
            rm1(mp,j+2) = real(0.0E+00,kind=8)
          end if

        end do

        rm1(i,i) = rm1(i,i) - y * ilambd

        if ( i >= n1 ) then
          l = i - ns
          z(mp,l) = -ilambd
          z(i,l) = z(i,l) + x * ilambd
        else
          rm1(mp,i+2) = -ilambd
          rm1(i,i+2) = rm1(i,i+2) + x * ilambd
        end if

        go to 640

580     continue

        if ( x == real(0.0E+00,kind=8) ) then
          rm1(mp,mp) = eps3
          if ( i < n ) then
            rm1(mp,i+1) = real(0.0E+00,kind=8)
          else if ( i == n ) then
            z(mp,s-1) = real(0.0E+00,kind=8)
          end if
          t = real(0.0E+00,kind=8)
          x = eps3**2
        end if

        w = w / x
        x = rm1(mp,mp) * w
        y = -t * w

        do j = i, uk
          if ( j >= n1 ) then
            l = j - ns
            t = z(mp,l)
            z(i,l) = -x * t - y * rm1(j,mp)
          else
            t = rm1(mp,j+2)
            rm1(i,j+2) = -x * t - y * rm1(j,mp)
          end if
          rm1(j,i) = rm1(j,i) - x * rm1(j,mp) + y * t
        end do

        if ( i >= n1 ) then
          l = i - ns
          z(i,l) = z(i,l) - ilambd
        else
          rm1(i,i+2) = rm1(i,i+2) - ilambd
        end if

640    continue

     end do

     if ( uk >= n1 ) then
       l = uk - ns
       t = z(uk,l)
     else
       t = rm1(uk,uk+2)
     end if

     if ( rm1(uk,uk) == real(0.0E+00,kind=8) .and. t == real(0.0E+00,kind=8) ) then
       rm1(uk,uk) = eps3
     end if
!
!  Back substitution for complex vector.
!
660  continue

     do ii = 1, uk

        i = uk + 1 - ii
        x = rv1(i)
        y = real(0.0E+00,kind=8)

        do j = i+1, uk

          if ( j >= n1 ) then
            t = z(i,j-ns)
          else
            t = rm1(i,j+2)
          end if

          x = x - rm1(j,i) * rv1(j) + t * rv2(j)
          y = y - rm1(j,i) * rv2(j) - t * rv1(j)

        end do

        if ( i >= n1 ) then
          t = z(i,i-ns)
        else
          t = rm1(i,i+2)
        end if

       call cdiv ( x, y, rm1(i,i), t, rv1(i), rv2(i) )

     end do
!
!  Acceptance test for real(kind=8) or complex eigenvector and normalization.
!
740  continue

     its = its + 1
     norm = real(0.0E+00,kind=8)
     normv = real(0.0E+00,kind=8)

     do i = 1, uk
       if ( ilambd == real(0.0E+00,kind=8) ) then
         x = abs ( rv1(i) )
       else
         x = pythag ( rv1(i), rv2(i) )
       end if
       if ( normv < x )  then
         normv = x
         j = i
       end if
       norm = norm + x
     end do

     if ( norm < growto ) go to 840
!
!  Accept vector.
!
     x = rv1(j)
     if ( ilambd == real(0.0E+00,kind=8) ) then
       x = real(1.0E+00,kind=8) / x
     else
       y = rv2(j)
     end if

     do i = 1, uk
       if ( ilambd == real(0.0E+00,kind=8) ) then
         z(i,s) = rv1(i) * x
       else
         call cdiv ( rv1(i), rv2(i), x, y, z(i,s-1), z(i,s) )
       end if
     end do

     if ( uk == n ) go to 940
     j = uk + 1
     go to 900
!
!  Choose a new starting vector.
!
840  continue

     if ( its >= uk ) go to 880

     x = ukroot
     y = eps3 / ( x + real(1.0E+00,kind=8) )

     rv1(1) = eps3
     rv1(2:uk) = y
     
     j = uk - its + 1
     rv1(j) = rv1(j) - eps3 * x
     if ( ilambd == real(0.0E+00,kind=8) ) go to 440
     go to 660
!
!  Set error: unaccepted eigenvector.
!
880  continue

     j = 1
     ierr = -k
!
!  Set remaining vector components to zero.
!
900  continue

     do i = j, n
       z(i,s) = real(0.0E+00,kind=8)
       if ( ilambd /= real(0.0E+00,kind=8) ) z(i,s-1) = real(0.0E+00,kind=8)
     end do

940  continue

     s = s + 1

960  continue

     if ( ip == (-1) ) then
       ip = 0
     end if

     if ( ip == 1 ) then
       ip = -1
     end if

  end do

  go to 1001
!
!  Set error: underestimate of eigenvector space required.
!
1000 continue

  if ( ierr /= 0 ) then
    ierr = ierr - n
  end if

  if ( ierr == 0 ) then
    ierr = -(2 * n + 1)
  end if

1001 continue

  m = s - 1 - abs ( ip )

  return
end
subroutine minfit ( nm, m, n, a, w, ip, b, ierr )
!
!*******************************************************************************
!
!! MINFIT solves the least squares problem, for a real(kind=8) overdetermined linear system.
!
!
!  Discussion:
!
!    This subroutine is part of an algorithm for solving general linear
!    systems of the form A*X=B.  
!
!    It determines the singular value decomposition 
!      A = U * S * V' 
!    of a real(kind=8) M by N rectangular matrix, forming U' * B 
!    rather than U.  Householder bidiagonalization and a variant of the 
!    QR algorithm are used.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!    NM must be at least as large as the maximum of M and N.
!
!    Input, integer M, the number of rows of A and B.
!
!    Input, integer N, the number of columns of A, and the order of V.
!
!    Input/output, real(kind=8) A(NM,N). On input, the rectangular coefficient matrix.
!    On output, A has been overwritten by the orthogonal matrix V of the
!    decomposition in its first N rows and columns.  If an error exit is made, 
!    the columns of V corresponding to indices of correct singular values 
!    should be correct.
!
!    Output, real(kind=8) W(N), the singular values of A.  These are the diagonal 
!    elements of S.  They are unordered.  If an error exit is made, the 
!    singular values should be correct for indices IERR+1, IERR+2,...,N.
!
!    Input, integer IP, is the number of columns of B.  IP can be zero.
!
!    Input/output, real(kind=8) B(NM,IP).  On input, the constant column matrix,
!    On output, B has been overwritten by U'*B.  If an error exit is made,
!    the rows of U'*B corresponding to indices of correct singular values 
!    should be correct.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    K, if the K-th singular value has not been determined after 30 iterations.
!
  implicit none
!
  integer ip
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) b(nm,ip)
  real(kind=8) c
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  integer i
  integer i1
  integer ierr
  integer ii
  integer its
  integer j
  integer k
  integer k1
  integer kk
  integer l
  integer l1
  integer ll
  integer m
  integer m1
  real(kind=8) pythag
  real(kind=8) rv1(n)
  real(kind=8) s
  real(kind=8) scale
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) w(n)
  real(kind=8) x
  real(kind=8) y
  real(kind=8) z
!
  ierr = 0
!
!  Householder reduction to bidiagonal form.
!
  g = real(0.0E+00,kind=8)
  scale = real(0.0E+00,kind=8)
  x = real(0.0E+00,kind=8)

  do i = 1, n

    l = i + 1
    rv1(i) = scale * g
    g = real(0.0E+00,kind=8)
    s = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)

    if ( i <= m ) then

      scale = sum ( abs ( a(i:m,i) ) )

      if ( scale /= real(0.0E+00,kind=8) ) then

        a(i:m,i) = a(i:m,i) / scale

        s = s + sum ( a(i:m,i)**2 )

        f = a(i,i)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        a(i,i) = f - g

        do j = l, n

          s = dot_product ( a(i:m,i), a(i:m,j) )

          f = s / h
          a(i:m,j) = a(i:m,j) + f * a(i:m,i)

        end do

        do j = 1, ip

          s = dot_product ( a(i:m,i), b(i:m,j) )

          b(i:m,j) = b(i:m,j) + s * a(i:m,i) / h

        end do

        a(i:m,i) = scale * a(i:m,i)

      end if

    end if

    w(i) = scale * g
    g = real(0.0E+00,kind=8)
    s = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)

    if ( i <= m .and. i /= n ) then

      do k = l, n
        scale = scale + abs ( a(i,k) )
      end do

      if ( scale /= real(0.0E+00,kind=8) ) then

        a(i,l:n) = a(i,l:n) / scale

        s = s + sum ( a(i,l:n)**2 )

        f = a(i,l)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        a(i,l) = f - g
        rv1(l:n) = a(i,l:n) / h

        do j = l, m

          s = dot_product ( a(j,l:n), a(i,l:n) )

          a(j,l:n) = a(j,l:n) + s * rv1(l:n)

        end do

        a(i,l:n) = scale * a(i,l:n)

      end if

    end if

    x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

  end do
!
!  Accumulation of right-hand transformations.
!
  do ii = 1, n

    i = n + 1 - ii

    if ( i /= n ) then

      if ( g /= real(0.0E+00,kind=8) ) then

        a(l:n,i) = ( a(i,l:n) / a(i,l) ) / g

        do j = l, n

          s = dot_product ( a(i,l:n), a(l:n,j) )

          a(l:n,j) = a(l:n,j) + s * a(l:n,i)

        end do

      end if

      a(i,l:n) = real(0.0E+00,kind=8)
      a(l:n,i) = real(0.0E+00,kind=8)

    end if

    a(i,i) = real(1.0E+00,kind=8)
    g = rv1(i)
    l = i

  end do

  if ( m < n .and. ip /= 0 ) then
    m1 = m + 1
    b(m+1:n,1:ip) = real(0.0E+00,kind=8)
  end if
!
!  Diagonalization of the bidiagonal form.
!
  tst1 = x

  do kk = 1, n

    k1 = n - kk
    k = k1 + 1
    its = 0
!
!  Test for splitting.
!
520 continue

    do ll = 1, k

      l1 = k - ll
      l = l1 + 1
      tst2 = tst1 + abs ( rv1(l) )

      if ( tst2 == tst1 ) then
        go to 565
      end if

      tst2 = tst1 + abs ( w(k-ll) )

      if ( tst2 == tst1 ) then
        exit
      end if

    end do
!
!  Cancellation of RV1(l) if l greater than 1.
!
540 continue

    c = real(0.0E+00,kind=8)
    s = real(1.0E+00,kind=8)

    do i = l, k

      f = s * rv1(i)
      rv1(i) = c * rv1(i)
      tst2 = tst1 + abs ( f)

      if ( tst2 == tst1 ) then
        exit
      end if

      g = w(i)
      h = pythag ( f, g )
      w(i) = h
      c = g / h
      s = -f / h

      do j = 1, ip
        y = b(l1,j)
        z = b(i,j)
        b(l1,j) = y * c + z * s
        b(i,j) = -y * s + z * c
      end do

    end do
!
!  Test for convergence.
!
565 continue

    z = w(k)

    if ( l == k ) go to 650
!
!  Shift from bottom 2 by 2 minor.
!
     if ( its >= 30 ) then
       ierr = k
       return
     end if

     its = its + 1
     x = w(l)
     y = w(k1)
     g = rv1(k1)
     h = rv1(k)
     f = real(0.5E+00,kind=8) * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y )
     g = pythag ( f, real(1.0E+00,kind=8) )
     f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + sign ( g, f ) ) - h )
!
!  Next QR transformation.
!
     c = real(1.0E+00,kind=8)
     s = real(1.0E+00,kind=8)

     do i1 = l, k1

        i = i1 + 1
        g = rv1(i)
        y = w(i)
        h = s * g
        g = c * g
        z = pythag ( f, h )
        rv1(i1) = z
        c = f / z
        s = h / z
        f = x * c + g * s
        g = -x * s + g * c
        h = y * s
        y = y * c

        do j = 1, n
          x = a(j,i1)
          z = a(j,i)
          a(j,i1) = x * c + z * s
          a(j,i) = -x * s + z * c
        end do

        z = pythag ( f, h )
        w(i1) = z

        if ( z /= real(0.0E+00,kind=8) ) then
          c = f / z
          s = h / z
        end if

        f = c * g + s * y
        x = -s * g + c * y

        do j = 1, ip
          y = b(i1,j)
          z = b(i,j)
          b(i1,j) = y * c + z * s
          b(i,j) = -y * s + z * c
        end do

     end do

     rv1(l) = real(0.0E+00,kind=8)
     rv1(k) = f
     w(k) = x
     go to 520
!
!  Convergence.
!
650 continue

    if ( z < real(0.0E+00,kind=8) ) then
      w(k) = - z
      a(1:n,k) = - a(1:n,k)
    end if

  end do

  return
end
subroutine ortbak ( nm, low, igh, a, ort, m, z )
!
!*******************************************************************************
!
!! ORTBAK determines eigenvectors by undoing the ORTHES transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a real(kind=8) general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by ORTHES.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer LOW, IGH, are determined by the balancing routine BALANC.  
!    If BALANC has not been used, set LOW = 1 and IGH equal to the order of 
!    the matrix.
!
!    Input, real(kind=8) A(NM,IGH), contains information about the orthogonal 
!    transformations used in the reduction by ORTHES in its strict 
!    lower triangle.
!
!    Input/output, real(kind=8) ORT(IGH), contains further information about the 
!    transformations used in the reduction by ORTHES.  On output, ORT
!    has been altered.
!
!    Input, integer M, the number of columns of Z to be back transformed.
!
!    Input/output, real(kind=8) Z(NM,N).  On input, the real(kind=8) and imaginary parts of 
!    the eigenvectors to be back transformed in the first M columns.  On
!    output, the real(kind=8) and imaginary parts of the transformed eigenvectors.
!
  implicit none
!
  integer igh
  integer m
  integer nm
!
  real(kind=8) a(nm,igh)
  real(kind=8) g
  integer i
  integer j
  integer low
  integer mp
  real(kind=8) ort(igh)
  real(kind=8) z(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  do mp = igh - 1, low + 1, -1

    if ( a(mp,mp-1) /= real(0.0E+00,kind=8) ) then

      ort(mp+1:igh) = a(mp+1:igh,mp-1)

      do j = 1, m

        g = dot_product ( ort(mp:igh), z(mp:igh,j) )

        g = ( g / ort(mp) ) / a(mp,mp-1)

        do i = mp, igh
          z(i,j) = z(i,j) + g * ort(i)
        end do

      end do

    end if

  end do

  return
end
subroutine orthes ( nm, n, low, igh, a, ort )
!
!*******************************************************************************
!
!! ORTHES transforms a real(kind=8) general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!    Given a real(kind=8) general matrix, this subroutine reduces a submatrix 
!    situated in rows and columns LOW through IGH to upper Hessenberg form by
!    orthogonal similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine BALANC.  
!    If BALANC has not been used, set LOW = 1 and IGH = N.
!
!    Input/output, real(kind=8) A(NM,N).  On input, the matrix.  On output,
!    the Hessenberg matrix.  Information about the orthogonal transformations 
!    used in the reduction is stored in the remaining triangle under the
!    Hessenberg matrix.
!
!    Output, real(kind=8) ORT(IGH), contains further information about the 
!    transformations.
!
  implicit none
!
  integer igh
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  integer i
  integer ii
  integer j
  integer jj
  integer la
  integer low
  integer m
  integer mp
  real(kind=8) ort(igh)
  real(kind=8) scale
!
  la = igh - 1

  do m = low + 1, la

     h = real(0.0E+00,kind=8)
     ort(m) = real(0.0E+00,kind=8)
     scale = real(0.0E+00,kind=8)
!
!  Scale the column.
!
     do i = m, igh
       scale = scale + abs ( a(i,m-1) )
     end do

     if ( scale /= real(0.0E+00,kind=8) ) then

     mp = m + igh

     do ii = m, igh
       i = mp - ii
       ort(i) = a(i,m-1) / scale
       h = h + ort(i) * ort(i)
     end do

     g = - sign ( sqrt ( h ), ort(m) )
     h = h - ort(m) * g
     ort(m) = ort(m) - g
!
!  Form (I-(U*Ut)/h) * A.
!
     do j = m, n

        f = real(0.0E+00,kind=8)

        do ii = m, igh
          i = mp - ii
          f = f + ort(i) * a(i,j)
        end do

        f = f / h

        do i = m, igh
          a(i,j) = a(i,j) - f * ort(i)
        end do

     end do
!
!  Form (I-(u*ut)/h) * A * (I-(u*ut)/h).
!
     do i = 1, igh

        f = real(0.0E+00,kind=8)
        do jj = m, igh
          j = mp - jj
          f = f + ort(j) * a(i,j)
        end do

        a(i,m:igh) = a(i,m:igh) - f * ort(m:igh) / h

     end do

     ort(m) = scale * ort(m)
     a(m,m-1) = scale * g

    end if

  end do

  return
end
subroutine ortran ( nm, n, low, igh, a, ort, z )
!
!*******************************************************************************
!
!! ORTRAN accumulates similarity transformations generated by ORTHES.
!
!
!  Discussion:
!
!    This subroutine accumulates the orthogonal similarity
!    transformations used in the reduction of a real(kind=8) general
!    matrix to upper Hessenberg form by ORTHES.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing routine BALANC.  
!    If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input, real(kind=8) A(NM,IGH), contains information about the orthogonal 
!    transformations used in the reduction by ORTHES in its strict lower 
!    triangle.
!
!    Input/output, real(kind=8) ORT(IGH), contains further information about the 
!    transformations used in the reduction by ORTHES.  On output, ORT
!    has been further altered.
!
!    Output, real(kind=8) Z(NM,N), contains the transformation matrix produced in the
!    reduction by ORTHES.
!
  implicit none
!
  integer igh
  integer n
  integer nm
!
  real(kind=8) a(nm,igh)
  real(kind=8) g
  integer i
  integer j
  integer kl
  integer low
  integer mm
  integer mp
  real(kind=8) ort(igh)
  real(kind=8) z(nm,n)
!
!  Initialize Z to the identity matrix.
!
  call rmat_ident ( nm, n, z )

  kl = igh - low - 1

  if ( kl < 1 ) then
    return
  end if

  do mm = 1, kl

    mp = igh - mm

    if ( a(mp,mp-1) /= real(0.0E+00,kind=8) ) then

      ort(mp+1:igh) = a(mp+1:igh,mp-1)

      do j = mp, igh

        g = real(0.0E+00,kind=8)
        do i = mp, igh
          g = g + ort(i) * z(i,j)
        end do

        g = ( g / ort(mp) ) / a(mp,mp-1)

        z(mp:igh,j) = z(mp:igh,j) + g * ort(mp:igh)

      end do

    end if

  end do

  return
end
function pythag ( a, b )
!
!*******************************************************************************
!
!! PYTHAG computes SQRT ( A**2 + B**2 ) carefully.
!
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A**2 + B**2 )
!
!    is reasonably accurate, but can fail if, for example, A**2 is larger 
!    than the machine overflow.  The formula can lose most of its accuracy 
!    if the sum of the squares is very large or very small.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Modified:
!
!    02 March 2000
!
!  Parameters:
!
!    Input, real(kind=8) A, B, the two legs of a right triangle.
!
!    Output, real(kind=8) PYTHAG, the length of the hypotenuse.
!
  implicit none
!
  real(kind=8) a
  real(kind=8) b
  real(kind=8) p
  real(kind=8) pythag
  real(kind=8) r
  real(kind=8) s
  real(kind=8) t
  real(kind=8) u
!
  p = max ( abs ( a ), abs ( b ) )

  if ( p /= real(0.0E+00,kind=8) ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = real(4.0E+00,kind=8) + r

      if ( t == real(4.0E+00,kind=8) ) then
        exit
      end if

      s = r / t
      u = real(1.0E+00,kind=8) + real(2.0E+00,kind=8) * s
      p = u * p
      r = ( s / u )**2 * r

    end do

  end if

  pythag = p

  return
end
subroutine qzhes ( nm, n, a, b, matz, z )
!
!*******************************************************************************
!
!! QZHES carries out transformations for a generalized eigenvalue problem.
!
!
!  Discussion:
!
!    This subroutine is the first step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This subroutine accepts a pair of real(kind=8) general matrices and
!    reduces one of them to upper Hessenberg form and the other
!    to upper triangular form using orthogonal transformations.
!    it is usually followed by QZIT, QZVAL and, possibly, QZVEC.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices.
!
!    Input/output, real(kind=8) A(NM,N).  On input, the first real(kind=8) general matrix.
!    On output, A has been reduced to upper Hessenberg form.  The elements
!    below the first subdiagonal have been set to zero.
!
!    Input/output, real(kind=8) B(NM,N).  On input, a real(kind=8) general matrix.
!    On output, B has been reduced to upper triangular form.  The elements
!    below the main diagonal have been set to zero.
!
!    Input, logical MATZ, should be TRUE if the right hand transformations
!    are to be accumulated for later use in computing eigenvectors.
!
!    Output, real(kind=8) Z(NM,N), contains the product of the right hand 
!    transformations if MATZ is TRUE.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) b(nm,n)
  integer i
  integer j
  integer k
  integer l
  integer l1
  integer lb
  logical matz
  integer nk1
  integer nm1
  real(kind=8) r
  real(kind=8) rho
  real(kind=8) s
  real(kind=8) t
  real(kind=8) u1
  real(kind=8) u2
  real(kind=8) v1
  real(kind=8) v2
  real(kind=8) z(nm,n)
!
!  Set Z to the identity matrix.
!
  if ( matz ) then
    call rmat_ident ( nm, n, z )
  end if
!
!  Reduce B to upper triangular form.
!
  if ( n <= 1 ) then
    return
  end if

  nm1 = n - 1

  do l = 1, n-1

    l1 = l + 1

    s = sum ( abs ( b(l+1:n,l) ) )

    if ( s /= real(0.0E+00,kind=8) ) then

      s = s + abs ( b(l,l) )
      b(l:n,l) = b(l:n,l) / s

      r = sqrt ( sum ( b(l:n,l)**2 ) )
      r = sign ( r, b(l,l) )
      b(l,l) = b(l,l) + r
      rho = r * b(l,l)

      do j = l+1, n

        t = dot_product ( b(l:n,l), b(l:n,j) )

        b(l:n,j) = b(l:n,j) - t * b(l:n,l) / rho

      end do

      do j = 1, n

        t = dot_product ( b(l:n,l), a(l:n,j) )

        a(l:n,j) = a(l:n,j) - t * b(l:n,l) / rho

      end do

      b(l,l) = -s * r
      b(l+1:n,l) = real(0.0E+00,kind=8)

    end if

  end do
!
!  Reduce A to upper Hessenberg form, while keeping B triangular.
!
  if ( n == 2 ) then
    return
  end if

  do k = 1, n-2

     nk1 = nm1 - k

     do lb = 1, nk1

        l = n - lb
        l1 = l + 1
!
!  Zero A(l+1,k).
!
        s = abs ( a(l,k) ) + abs ( a(l1,k) )

        if ( s /= real(0.0E+00,kind=8) ) then

        u1 = a(l,k) / s
        u2 = a(l1,k) / s
        r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
        v1 =  -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1

        do j = k, n
          t = a(l,j) + u2 * a(l1,j)
          a(l,j) = a(l,j) + t * v1
          a(l1,j) = a(l1,j) + t * v2
        end do

        a(l1,k) = real(0.0E+00,kind=8)

        do j = l, n
          t = b(l,j) + u2 * b(l1,j)
          b(l,j) = b(l,j) + t * v1
          b(l1,j) = b(l1,j) + t * v2
        end do
!
!  Zero B(l+1,l).
!
        s = abs ( b(l1,l1) ) + abs ( b(l1,l) )

        if ( s /= 0.0 ) then

        u1 = b(l1,l1) / s
        u2 = b(l1,l) / s
        r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
        v1 =  -( u1 + r ) / r
        v2 = -u2 / r
        u2 = v2 / v1

        do i = 1, l1
          t = b(i,l1) + u2 * b(i,l)
          b(i,l1) = b(i,l1) + t * v1
          b(i,l) = b(i,l) + t * v2
        end do

        b(l1,l) = real(0.0E+00,kind=8)

        do i = 1, n
          t = a(i,l1) + u2 * a(i,l)
          a(i,l1) = a(i,l1) + t * v1
          a(i,l) = a(i,l) + t * v2
        end do

        if ( matz ) then

          do i = 1, n
            t = z(i,l1) + u2 * z(i,l)
            z(i,l1) = z(i,l1) + t * v1
            z(i,l) = z(i,l) + t * v2
          end do

        end if

        end if

      end if

    end do

  end do

  return
end
subroutine qzit ( nm, n, a, b, eps1, matz, z, ierr )
!
!*******************************************************************************
!
!! QZIT carries out iterations to solve a generalized eigenvalue problem.
!
!
!  Discussion:
!
!    This subroutine is the second step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This subroutine accepts a pair of real(kind=8) matrices, one of them
!    in upper Hessenberg form and the other in upper triangular form.
!    It reduces the Hessenberg matrix to quasi-triangular form using
!    orthogonal transformations while maintaining the triangular form
!    of the other matrix.  It is usually preceded by QZHES and
!    followed by QZVAL and, possibly, QZVEC.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices.
!
!    Input/output, real(kind=8) A(NM,N).  On input, a real(kind=8) upper Hessenberg matrix.
!    On output, A has been reduced to quasi-triangular form.  The elements
!    below the first subdiagonal are still zero and no two consecutive 
!    subdiagonal elements are nonzero.
!
!    Input/output, real(kind=8) B(NM,N).  On input, a real(kind=8) upper triangular matrix.
!    On output, B is still in upper triangular form, although its elements
!    have been altered.  The location B(N,1) is used to store EPS1 times 
!    the norm of B for later use by QZVAL and QZVEC.
!
!    Input, real(kind=8) EPS1, a tolerance used to determine negligible elements.
!    EPS1 = real(0.0E+00,kind=8) (or negative) may be input, in which case an element 
!    will be neglected only if it is less than roundoff error times the 
!    norm of its matrix.  If the input EPS1 is positive, then an element 
!    will be considered negligible if it is less than EPS1 times the norm 
!    of its matrix.  A positive value of EPS1 may result in faster execution,
!    but less accurate results.
!
!    Input, logical MATZ, should be TRUE if the right hand transformations
!    are to be accumulated for later use in computing eigenvectors.
!
!    Input/output, real(kind=8) Z(NM,N).  If MATZ is FALSE, Z is not referenced.
!    Otherwise, on input, the transformation matrix produced in the reduction 
!    by QZHES, if performed, or else the identity matrix.  On output, Z 
!    contains the product of the right hand transformations for both steps.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    J, if the limit of 30*N iterations is exhausted while the J-th 
!      eigenvalue is being sought.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) a1
  real(kind=8) a11
  real(kind=8) a12
  real(kind=8) a2
  real(kind=8) a21
  real(kind=8) a22
  real(kind=8) a3
  real(kind=8) a33
  real(kind=8) a34
  real(kind=8) a43
  real(kind=8) a44
  real(kind=8) ani
  real(kind=8) anorm
  real(kind=8) b(nm,n)
  real(kind=8) b11
  real(kind=8) b12
  real(kind=8) b22
  real(kind=8) b33
  real(kind=8) b34
  real(kind=8) b44
  real(kind=8) bni
  real(kind=8) bnorm
  integer en
  integer enm2
  integer enorn
  real(kind=8) ep
  real(kind=8) eps1
  real(kind=8) epsa
  real(kind=8) epsb
  integer i
  integer ierr
  integer ish
  integer itn
  integer its
  integer j
  integer k
  integer k1
  integer k2
  integer km1
  integer l
  integer l1
  integer ld
  integer ll
  integer lm1
  integer lor1
  logical matz
  integer na
  logical notlas
  real(kind=8) r
  real(kind=8) s
  real(kind=8) sh
  real(kind=8) t
  real(kind=8) u1
  real(kind=8) u2
  real(kind=8) u3
  real(kind=8) v1
  real(kind=8) v2
  real(kind=8) v3
  real(kind=8) z(nm,n)
!
  ierr = 0
!
!  Compute EPSA and EPSB.
!
  anorm = real(0.0E+00,kind=8)
  bnorm = real(0.0E+00,kind=8)

  do i = 1, n

    ani = real(0.0E+00,kind=8)
    if ( i /= 1 ) ani = abs ( a(i,i-1) )
    bni = real(0.0E+00,kind=8)

    do j = i, n
      ani = ani + abs ( a(i,j) )
      bni = bni + abs ( b(i,j) )
    end do

    anorm = max ( anorm, ani )
    bnorm = max ( bnorm, bni )

  end do

  if ( anorm == real(0.0E+00,kind=8) ) anorm = real(1.0E+00,kind=8)
  if ( bnorm == real(0.0E+00,kind=8) ) bnorm = real(1.0E+00,kind=8)
  ep = eps1

  if ( ep > real(0.0E+00,kind=8) ) go to 50
!
!  Use roundoff level if EPS1 is 0.
!
  ep = epsilon ( ep )

50 continue

  epsa = ep * anorm
  epsb = ep * bnorm
!
!  Reduce A to quasi-triangular form, while keeping B triangular.
!
  lor1 = 1
  enorn = n
  en = n
  itn = 30 * n
!
!  Begin QZ step.
!
60 continue

  if ( en <= 2 ) go to 1001
  if (.not. matz) enorn = en
  its = 0
  na = en - 1
  enm2 = na - 1

70 continue

  ish = 2
!
!  Check for convergence or reducibility.
!
  do ll = 1, en
    lm1 = en - ll
    l = lm1 + 1
    if ( l == 1 ) go to 95
    if ( abs ( a(l,lm1) ) <= epsa ) then
      exit
    end if
  end do

90 continue

  a(l,lm1) = real(0.0E+00,kind=8)
  if ( l < na ) go to 95
!
!  1-by-1 or 2-by-2 block isolated.
!
  en = lm1
  go to 60
!
!  Check for small top of B.
!
95 continue

  ld = l

100 continue

  l1 = l + 1
  b11 = b(l,l)

  if ( abs ( b11 ) > epsb ) go to 120

  b(l,l) = real(0.0E+00,kind=8)
  s = abs ( a(l,l) ) + abs ( a(l1,l) )
  u1 = a(l,l) / s
  u2 = a(l1,l) / s
  r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
  v1 = - ( u1 + r ) / r
  v2 = -u2 / r
  u2 = v2 / v1

  do j = l, enorn
    t = a(l,j) + u2 * a(l1,j)
    a(l,j) = a(l,j) + t * v1
    a(l1,j) = a(l1,j) + t * v2
    t = b(l,j) + u2 * b(l1,j)
    b(l,j) = b(l,j) + t * v1
    b(l1,j) = b(l1,j) + t * v2
  end do

  if ( l /= 1 ) then
    a(l,lm1) = -a(l,lm1)
  end if
  lm1 = l
  l = l1
  go to 90

120 continue

  a11 = a(l,l) / b11
  a21 = a(l1,l) / b11
  if ( ish == 1 ) go to 140
!
!  Iteration strategy.
!
  if ( itn == 0 ) go to 1000
  if ( its == 10 ) go to 155
!
!  Determine type of shift.
!
  b22 = b(l1,l1)
  if ( abs ( b22 ) < epsb ) b22 = epsb
  b33 = b(na,na)
  if ( abs ( b33 ) < epsb ) b33 = epsb
  b44 = b(en,en)
  if ( abs ( b44 ) < epsb ) b44 = epsb
  a33 = a(na,na) / b33
  a34 = a(na,en) / b44
  a43 = a(en,na) / b33
  a44 = a(en,en) / b44
  b34 = b(na,en) / b44
  t = real(0.5E+00,kind=8) * (a43 * b34 - a33 - a44)
  r = t * t + a34 * a43 - a33 * a44

  if ( r < real(0.0E+00,kind=8) ) go to 150
!
!  Determine single shift zeroth column of A.
!
  ish = 1
  r = sqrt ( r )
  sh = -t + r
  s = -t - r
  if ( abs ( s - a44 ) < abs ( sh - a44 ) ) sh = s
!
!  Look for two consecutive small sub-diagonal elements of A.
! 
  do ll = ld, enm2
    l = enm2 + ld - ll
    if ( l == ld ) then
      exit
    end if
    lm1 = l - 1
    l1 = l + 1
    t = a(l,l)
    if ( abs ( b(l,l) ) > epsb ) t = t - sh * b(l,l)
    if ( abs ( a(l,lm1) ) <= abs ( t / a(l1,l) ) * epsa ) go to 100
  end do

140 continue

  a1 = a11 - sh
  a2 = a21

  if ( l /= ld ) then
    a(l,lm1) = -a(l,lm1)
  end if

  go to 160
!
!  Determine double shift zeroth column of A.
!
150 continue

  a12 = a(l,l1) / b22
  a22 = a(l1,l1) / b22
  b12 = b(l,l1) / b22
  a1 = ( ( a33 - a11 ) * ( a44 - a11 ) - a34 * a43 + a43 * b34 * a11 ) &
    / a21 + a12 - a11 * b12
  a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34
  a3 = a(l1+1,l1) / b22
  go to 160
!
!  Ad hoc shift.
!
155 continue

  a1 = real(0.0E+00,kind=8)
  a2 = real(1.0E+00,kind=8)
  a3 = 1.1605E+00

  160 continue
  its = its + 1
  itn = itn - 1
  if ( .not. matz ) lor1 = ld
!
!  Main loop.
!
  do k = l, na

     notlas = k /= na .and. ish == 2
     k1 = k + 1
     k2 = k + 2
     km1 = max ( k-1, l )
     ll = min ( en, k1+ish )

     if ( notlas ) go to 190
!
!  Zero A(k+1,k-1).
!
     if ( k /= l ) then
       a1 = a(k,km1)
       a2 = a(k1,km1)
     end if

     s = abs ( a1 ) + abs ( a2 )

     if ( s == real(0.0E+00,kind=8) ) go to 70

     u1 = a1 / s
     u2 = a2 / s
     r = sign ( sqrt ( u1**2 + u1**2 ), u1 )
     v1 = -( u1 + r ) / r
     v2 = -u2 / r
     u2 = v2 / v1

     do j = km1, enorn
       t = a(k,j) + u2 * a(k1,j)
       a(k,j) = a(k,j) + t * v1
       a(k1,j) = a(k1,j) + t * v2
       t = b(k,j) + u2 * b(k1,j)
       b(k,j) = b(k,j) + t * v1
       b(k1,j) = b(k1,j) + t * v2
     end do

     if ( k /= l ) then
       a(k1,km1) = real(0.0E+00,kind=8)
     end if

     go to 240
!
!  Zero A(k+1,k-1) and A(k+2,k-1).
!
190  continue

     if ( k /= l ) then
       a1 = a(k,km1)
       a2 = a(k1,km1)
       a3 = a(k2,km1)
     end if

     s = abs ( a1 ) + abs ( a2 ) + abs ( a3 )

     if ( s == real(0.0E+00,kind=8) ) go to 260

     u1 = a1 / s
     u2 = a2 / s
     u3 = a3 / s
     r = sign ( sqrt ( u1**2 + u2**2 + u3**2 ), u1 )
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     v3 = -u3 / r
     u2 = v2 / v1
     u3 = v3 / v1

     do j = km1, enorn
       t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
       a(k,j) = a(k,j) + t * v1
       a(k1,j) = a(k1,j) + t * v2
       a(k2,j) = a(k2,j) + t * v3
       t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
       b(k,j) = b(k,j) + t * v1
       b(k1,j) = b(k1,j) + t * v2
       b(k2,j) = b(k2,j) + t * v3
     end do

     if ( k /= l ) then
       a(k1,km1) = real(0.0E+00,kind=8)
       a(k2,km1) = real(0.0E+00,kind=8)
     end if
!
!  Zero B(k+2,k+1) and B(k+2,k).
!
     s = abs ( b(k2,k2) ) + abs ( b(k2,k1) ) + abs ( b(k2,k) )
     if ( s == real(0.0E+00,kind=8) ) go to 240
     u1 = b(k2,k2) / s
     u2 = b(k2,k1) / s
     u3 = b(k2,k) / s
     r = sign ( sqrt ( u1**2 + u2**2 + u3**2 ), u1 )
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     v3 = -u3 / r
     u2 = v2 / v1
     u3 = v3 / v1

     do i = lor1, ll
       t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
       a(i,k2) = a(i,k2) + t * v1
       a(i,k1) = a(i,k1) + t * v2
       a(i,k) = a(i,k) + t * v3
       t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
       b(i,k2) = b(i,k2) + t * v1
       b(i,k1) = b(i,k1) + t * v2
       b(i,k) = b(i,k) + t * v3
     end do

     b(k2,k) = real(0.0E+00,kind=8)
     b(k2,k1) = real(0.0E+00,kind=8)

     if ( matz ) then

       do i = 1, n
         t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
         z(i,k2) = z(i,k2) + t * v1
         z(i,k1) = z(i,k1) + t * v2
         z(i,k) = z(i,k) + t * v3
       end do

     end if
!
!  Zero B(k+1,k).
!
240  continue

     s = abs ( b(k1,k1) ) + abs ( b(k1,k) )

     if ( s /= real(0.0E+00,kind=8) ) then

     u1 = b(k1,k1) / s
     u2 = b(k1,k) / s
     r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
     v1 = -( u1 + r ) / r
     v2 = -u2 / r
     u2 = v2 / v1

     do i = lor1, ll
       t = a(i,k1) + u2 * a(i,k)
       a(i,k1) = a(i,k1) + t * v1
       a(i,k) = a(i,k) + t * v2
       t = b(i,k1) + u2 * b(i,k)
       b(i,k1) = b(i,k1) + t * v1
       b(i,k) = b(i,k) + t * v2
     end do

     b(k1,k) = real(0.0E+00,kind=8)

     if ( matz ) then

       do i = 1, n
         t = z(i,k1) + u2 * z(i,k)
         z(i,k1) = z(i,k1) + t * v1
         z(i,k) = z(i,k) + t * v2
       end do

     end if

   end if

260  continue

  end do

  go to 70
!
!  Set error: not all eigenvalues have converged after 30*N iterations.
!
1000 continue

  ierr = en
!
!  Save EPSB for use by QZVAL and QZVEC.
!
 1001 continue

  if ( n > 1 ) then
    b(n,1) = epsb
  end if

  return
end
subroutine qzval ( nm, n, a, b, alfr, alfi, beta, matz, z )
!
!*******************************************************************************
!
!! QZVAL computes eigenvalues for a generalized eigenvalue problem.
!
!
!  Discussion:
!
!    This subroutine is the third step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This subroutine accepts a pair of real(kind=8) matrices, one of them
!    in quasi-triangular form and the other in upper triangular form.
!    It reduces the quasi-triangular matrix further, so that any
!    remaining 2-by-2 blocks correspond to pairs of complex
!    eigenvalues, and returns quantities whose ratios give the
!    generalized eigenvalues.  It is usually preceded by QZHES
!    and QZIT and may be followed by QZVEC.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices.
!
!    Input/output, real(kind=8) A(NM,N).  On input, a real(kind=8) upper quasi-triangular 
!    matrix.  On output, A has been reduced further to a quasi-triangular 
!    matrix in which all nonzero subdiagonal elements correspond to
!    pairs of complex eigenvalues.
!
!    Input/output, real(kind=8) B(NM,N).  On input, a real(kind=8) upper triangular matrix.  
!    In addition, location B(n,1) contains the tolerance quantity EPSB
!    computed and saved in QZIT.  On output, B is still in upper triangular 
!    form, although its elements have been altered.  B(N,1) is unaltered.
!
!    Output, real(kind=8) ALFR(N), ALFI(N), the real(kind=8) and imaginary parts of the
!    diagonal elements of the triangular matrix that would be obtained 
!    if A were reduced completely to triangular form by unitary 
!    transformations.  Non-zero values of ALFI occur in pairs, the first 
!    member positive and the second negative.
!
!    Output, real(kind=8) BETA(N), the diagonal elements of the corresponding B,
!    normalized to be real(kind=8) and non-negative.  The generalized eigenvalues 
!    are then the ratios (ALFR + I * ALFI) / BETA.
!
!    Input, logical MATZ, should be TRUE if the right hand transformations
!    are to be accumulated for later use in computing eigenvectors, and 
!    to FALSE otherwise.
!
!    Input/output, real(kind=8) Z(NM,N), is only used if MATZ is TRUE.
!    On input, the transformation matrix produced in the reductions by QZHES
!    and QZIT, if performed, or else the identity matrix.  On output, 
!    the product of the right hand transformations for all three steps.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) a1
  real(kind=8) a11
  real(kind=8) a11i
  real(kind=8) a11r
  real(kind=8) a12
  real(kind=8) a12i
  real(kind=8) a12r
  real(kind=8) a1i
  real(kind=8) a2
  real(kind=8) a21
  real(kind=8) a22
  real(kind=8) a22i
  real(kind=8) a22r
  real(kind=8) a2i
  real(kind=8) an
  real(kind=8) alfi(n)
  real(kind=8) alfr(n)
  real(kind=8) b(nm,n)
  real(kind=8) b11
  real(kind=8) b12
  real(kind=8) b22
  real(kind=8) beta(n)
  real(kind=8) bn
  real(kind=8) c
  real(kind=8) cq
  real(kind=8) cz
  real(kind=8) d
  real(kind=8) di
  real(kind=8) dr
  real(kind=8) e
  real(kind=8) ei
  integer en
  real(kind=8) epsb
  integer i
  integer isw
  integer j
  logical matz
  integer na
  integer nn
  real(kind=8) r
  real(kind=8) s
  real(kind=8) sqi
  real(kind=8) sqr
  real(kind=8) ssi
  real(kind=8) ssr
  real(kind=8) szi
  real(kind=8) szr
  real(kind=8) t
  real(kind=8) ti
  real(kind=8) tr
  real(kind=8) u1
  real(kind=8) u2
  real(kind=8) v1
  real(kind=8) v2
  real(kind=8) z(nm,n)
!
  epsb = b(n,1)
  isw = 1
!
!  Find eigenvalues of quasi-triangular matrices.
!
  do nn = 1, n

     en = n + 1 - nn
     na = en - 1

     if ( isw == 2 ) go to 505

     if ( en == 1 ) go to 410

     if ( a(en,na) /= real(0.0E+00,kind=8) ) go to 420
!
!  1-by-1 block, one real(kind=8) root.
!
410  continue

     alfr(en) = a(en,en)
     if ( b(en,en) < real(0.0E+00,kind=8) ) alfr(en) = -alfr(en)
     beta(en) = abs ( b(en,en) )
     alfi(en) = real(0.0E+00,kind=8)
     go to 510
!
!  2-by-2 block.
!
420  continue

     if ( abs ( b(na,na) ) <= epsb ) then
       a1 = a(na,na)
       a2 = a(en,na)
       go to 460
     end if

     if ( abs ( b(en,en) ) <= epsb ) then
       a1 = a(en,en)
       a2 = a(en,na)
       bn = real(0.0E+00,kind=8)
       go to 435
     end if

     an = abs ( a(na,na) ) + abs ( a(na,en) ) + abs ( a(en,na) ) &
       + abs ( a(en,en) )
     bn = abs ( b(na,na) ) + abs ( b(na,en) ) + abs ( b(en,en) )
     a11 = a(na,na) / an
     a12 = a(na,en) / an
     a21 = a(en,na) / an
     a22 = a(en,en) / an
     b11 = b(na,na) / bn
     b12 = b(na,en) / bn
     b22 = b(en,en) / bn
     e = a11 / b11
     ei = a22 / b22
     s = a21 / ( b11 * b22 )
     t = ( a22 - e * b22 ) / b22

     if ( abs ( e ) > abs ( ei ) ) then
       e = ei
       t = ( a11 - e * b11 ) / b11
     end if

     c = real(0.5E+00,kind=8) * ( t - s * b12 )
     d = c**2 + s * ( a12 - e * b12 )
     if ( d < real(0.0E+00,kind=8) ) go to 480
!
!  Two real(kind=8) roots.
!  Zero both A(EN,NA) and B(EN,NA).
!
     e = e + ( c + sign ( sqrt ( d ), c ) )
     a11 = a11 - e * b11
     a12 = a12 - e * b12
     a22 = a22 - e * b22

     if ( abs ( a11 ) + abs ( a12 ) >= abs ( a21 ) + abs ( a22 ) ) then
       a1 = a12
       a2 = a11
     else
       a1 = a22
       a2 = a21
     end if
!
!  Choose and apply real(kind=8) Z.
!
435  continue

     s = abs ( a1 ) + abs ( a2 )
     u1 = a1 / s
     u2 = a2 / s
     r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
     v1 = - ( u1 + r ) / r
     v2 = - u2 / r
     u2 = v2 / v1

     do i = 1, en
       t = a(i,en) + u2 * a(i,na)
       a(i,en) = a(i,en) + t * v1
       a(i,na) = a(i,na) + t * v2
       t = b(i,en) + u2 * b(i,na)
       b(i,en) = b(i,en) + t * v1
       b(i,na) = b(i,na) + t * v2
     end do

     if ( matz ) then

       do i = 1, n
         t = z(i,en) + u2 * z(i,na)
         z(i,en) = z(i,en) + t * v1
         z(i,na) = z(i,na) + t * v2
       end do

     end if

450  continue

     if ( bn == real(0.0E+00,kind=8) ) go to 475

     if ( an >= abs ( e ) * bn ) then
       a1 = b(na,na)
       a2 = b(en,na)
     else
       a1 = a(na,na)
       a2 = a(en,na)
     end if
!
!  Choose and apply real(kind=8) Q.
!
460  continue

     s = abs ( a1 ) + abs ( a2 )
     if ( s == real(0.0E+00,kind=8) ) go to 475
     u1 = a1 / s
     u2 = a2 / s
     r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     u2 = v2 / v1

     do j = na, n
       t = a(na,j) + u2 * a(en,j)
       a(na,j) = a(na,j) + t * v1
       a(en,j) = a(en,j) + t * v2
       t = b(na,j) + u2 * b(en,j)
       b(na,j) = b(na,j) + t * v1
       b(en,j) = b(en,j) + t * v2
     end do

475  continue

     a(en,na) = real(0.0E+00,kind=8)
     b(en,na) = real(0.0E+00,kind=8)
     alfr(na) = a(na,na)
     alfr(en) = a(en,en)
     if ( b(na,na) < real(0.0E+00,kind=8) ) alfr(na) = -alfr(na)
     if ( b(en,en) < real(0.0E+00,kind=8) ) alfr(en) = -alfr(en)
     beta(na) = abs ( b(na,na) )
     beta(en) = abs ( b(en,en) )
     alfi(en) = real(0.0E+00,kind=8)
     alfi(na) = real(0.0E+00,kind=8)
     go to 505
!
!  Two complex roots.
!
480  continue

     e = e + c
     ei = sqrt ( -d )
     a11r = a11 - e * b11
     a11i = ei * b11
     a12r = a12 - e * b12
     a12i = ei * b12
     a22r = a22 - e * b22
     a22i = ei * b22

     if ( abs ( a11r ) + abs ( a11i ) + abs ( a12r ) + abs ( a12i ) >= &
            abs ( a21 ) + abs ( a22r ) + abs ( a22i ) ) then
       a1 = a12r
       a1i = a12i
       a2 = -a11r
       a2i = -a11i
     else
       a1 = a22r
       a1i = a22i
       a2 = -a21
       a2i = real(0.0E+00,kind=8)
     end if
!
!  Choose complex Z.
!
     cz = sqrt ( a1**2 + a1i**2 )

     if ( cz /= real(0.0E+00,kind=8) ) then
       szr = ( a1 * a2 + a1i * a2i) / cz
       szi = ( a1 * a2i - a1i * a2) / cz
       r = sqrt ( cz**2 + szr**2 + szi**2 )
       cz = cz / r
       szr = szr / r
       szi = szi / r
     else
       szr = real(1.0E+00,kind=8)
       szi = real(0.0E+00,kind=8)
     end if

     if ( an >= ( abs ( e ) + ei ) * bn ) then
       a1 = cz * b11 + szr * b12
       a1i = szi * b12
       a2 = szr * b22
       a2i = szi * b22
     else
       a1 = cz * a11 + szr * a12
       a1i = szi * a12
       a2 = cz * a21 + szr * a22
       a2i = szi * a22
     end if
!
!  Choose complex Q.
!
     cq = sqrt ( a1**2 + a1i**2 )

     if ( cq /= real(0.0E+00,kind=8) ) then
       sqr = ( a1 * a2 + a1i * a2i ) / cq
       sqi = ( a1 * a2i - a1i * a2 ) / cq
       r = sqrt ( cq**2 + sqr**2 + sqi**2 )
       cq = cq / r
       sqr = sqr / r
       sqi = sqi / r
     else
       sqr = real(1.0E+00,kind=8)
       sqi = real(0.0E+00,kind=8)
     end if
!
!  Compute diagonal elements that would result if transformations were applied.
!
     ssr = sqr * szr + sqi * szi
     ssi = sqr * szi - sqi * szr
     i = 1
     tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22
     ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
     dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
     di = cq * szi * b12 + ssi * b22
     go to 503

502  continue

     i = 2
     tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22
     ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
     dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
     di = -ssi * b11 - sqi * cz * b12

503  continue

     t = ti * dr - tr * di

     if ( t < real(0.0E+00,kind=8) ) then
       j = en
     else
       j = na
     end if

     r = sqrt ( dr**2 + di**2 )
     beta(j) = bn * r
     alfr(j) = an * (tr * dr + ti * di) / r
     alfi(j) = an * t / r

     if ( i == 1 ) go to 502

505  continue

     isw = 3 - isw

510  continue

  end do

  b(n,1) = epsb

  return
end
subroutine qzvec ( nm, n, a, b, alfr, alfi, beta, z )
!
!*******************************************************************************
!
!! QZVEC computes eigenvectors for a generalized eigenvalue problem.
!
!
!  Discussion:
!
!    This subroutine is the optional fourth step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This subroutine accepts a pair of real(kind=8) matrices, one of them in
!    quasi-triangular form (in which each 2-by-2 block corresponds to
!    a pair of complex eigenvalues) and the other in upper triangular
!    form.  It computes the eigenvectors of the triangular problem and
!    transforms the results back to the original coordinate system.
!    it is usually preceded by QZHES, QZIT, and QZVAL.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices.
!
!    Input, real(kind=8) A(NM,N), contains a real(kind=8) upper quasi-triangular matrix.
!    Its subdiagonal elements provide information about the storage of 
!    the complex eigenvectors.
!
!    Input/output, real(kind=8) B(NM,N).  On input, a real(kind=8) upper triangular matrix.  
!    In addition, location B(N,1) contains the tolerance quantity EPSB
!    computed and saved in QZIT.  On output, B has been destroyed.
!
!    Input, real(kind=8) ALFR(N), ALFI(N), BETA(N), vectors whose ratios 
!      ( ALFR + I * ALFI ) / BETA 
!    are the generalized eigenvalues.  They are usually obtained from QZVAL.
!
!    Input/output, real(kind=8) Z(NM,N).  On input, the transformation matrix produced 
!    in the reductions by QZHES, QZIT, and QZVAL, if performed.  If the 
!    eigenvectors of the triangular problem are desired, Z must contain the 
!    identity matrix.  On output, Z contains the real(kind=8) and imaginary parts of 
!    the eigenvectors:
!    If ALFI(I) == 0.0, the I-th eigenvalue is real(kind=8) and the I-th column of Z 
!    contains its eigenvector.
!    If ALFI(I) > 0.0, the eigenvalue is the first of a complex pair and the
!    I-th and (I+1)-th columns of Z contain its eigenvector.
!    If ALFI(I) < 0.0, the eigenvalue is the second of a complex pair and the 
!    (I-1)-th and I-th columns of Z contain the conjugate of its eigenvector.
!    Each eigenvector is normalized so that the modulus of its largest 
!    component is real(1.0E+00,kind=8) .
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) alfi(n)
  real(kind=8) alfm
  real(kind=8) alfr(n)
  real(kind=8) almi
  real(kind=8) almr
  real(kind=8) b(nm,n)
  real(kind=8) beta(n)
  real(kind=8) betm
  real(kind=8) d
  real(kind=8) di
  real(kind=8) dr
  integer en
  integer enm2
  real(kind=8) epsb
  integer i
  integer ii
  integer isw
  integer j
  integer jj
  integer k
  integer m
  integer na
  integer nn
  real(kind=8) q
  real(kind=8) r
  real(kind=8) ra
  real(kind=8) rr
  real(kind=8) s
  real(kind=8) sa
  real(kind=8) t
  real(kind=8) t1
  real(kind=8) t2
  real(kind=8) ti
  real(kind=8) tr
  real(kind=8) w
  real(kind=8) w1
  real(kind=8) x
  real(kind=8) x1
  real(kind=8) y
  real(kind=8) z(nm,n)
  real(kind=8) z1
  real(kind=8) zz
!
  epsb = b(n,1)
  isw = 1

  do nn = 1, n

     en = n + 1 - nn
     na = en - 1

     if ( isw == 2 ) go to 795

     if ( alfi(en) /= real(0.0E+00,kind=8) ) go to 710
!
!  real(kind=8) vector.
!
     m = en
     b(en,en) = real(1.0E+00,kind=8)

     if ( na == 0 ) go to 800

     alfm = alfr(m)
     betm = beta(m)

     do ii = 1, na

        i = en - ii
        w = betm * a(i,i) - alfm * b(i,i)
        r = real(0.0E+00,kind=8)

        do j = m, en
          r = r + ( betm * a(i,j) - alfm * b(i,j) ) * b(j,en)
        end do

        if ( i == 1 .or. isw == 2 ) go to 630

        if ( betm * a(i,i-1) == real(0.0E+00,kind=8) ) go to 630

        zz = w
        s = r
        go to 690

630     continue

        m = i

        if ( isw == 2 ) go to 640
!
!  real(kind=8) 1-by-1 block.
!
        t = w
        if ( w == real(0.0E+00,kind=8) ) t = epsb
        b(i,en) = - r / t
        go to 700
!
!  real(kind=8) 2-by-2 block.
!
640     continue

        x = betm * a(i,i+1) - alfm * b(i,i+1)
        y = betm * a(i+1,i)
        q = w * zz - x * y
        t = ( x * s - zz * r ) / q
        b(i,en) = t

        if ( abs ( x ) <= abs ( zz ) ) go to 650

        b(i+1,en) = (-r - w * t) / x

        go to 690

650     continue

        b(i+1,en) = (-s - y * t) / zz

690     continue

        isw = 3 - isw

700     continue

     end do
!
!  End real(kind=8) vector.
!
     go to 800
!
!  Complex vector.
!
710  continue

     m = na
     almr = alfr(m)
     almi = alfi(m)
     betm = beta(m)
!
!  Last vector component chosen imaginary so eigenvector matrix is triangular.
!
     y = betm * a(en,na)
     b(na,na) = -almi * b(en,en) / y
     b(na,en) = ( almr * b(en,en) - betm * a(en,en) ) / y
     b(en,na) = real(0.0E+00,kind=8)
     b(en,en) = real(1.0E+00,kind=8)
     enm2 = na - 1

     do ii = 1, enm2

        i = na - ii
        w = betm * a(i,i) - almr * b(i,i)
        w1 = -almi * b(i,i)
        ra = real(0.0E+00,kind=8)
        sa = real(0.0E+00,kind=8)

        do j = m, en
          x = betm * a(i,j) - almr * b(i,j)
          x1 = -almi * b(i,j)
          ra = ra + x * b(j,na) - x1 * b(j,en)
          sa = sa + x * b(j,en) + x1 * b(j,na)
        end do

        if ( i == 1 .or. isw == 2 ) go to 770
        if ( betm * a(i,i-1) == real(0.0E+00,kind=8) ) go to 770

        zz = w
        z1 = w1
        r = ra
        s = sa
        isw = 2
        go to 790
770     continue

        m = i
        if ( isw == 2 ) go to 780
!
!  Complex 1-by-1 block.
!
        tr = -ra
        ti = -sa

773     continue

        dr = w
        di = w1
!
!  Complex divide (t1,t2) = (tr,ti) / (dr,di),
!
775     continue

        if ( abs ( di ) > abs ( dr ) ) go to 777
        rr = di / dr
        d = dr + di * rr
        t1 = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        go to ( 787, 782 ), isw

777     continue

        rr = dr / di
        d = dr * rr + di
        t1 = ( tr * rr + ti ) / d
        t2 = ( ti * rr - tr ) / d
        go to ( 787, 782 ), isw
!
!  Complex 2-by-2 block.
!
780     continue

        x = betm * a(i,i+1) - almr * b(i,i+1)
        x1 = -almi * b(i,i+1)
        y = betm * a(i+1,i)
        tr = y * ra - w * r + w1 * s
        ti = y * sa - w * s - w1 * r
        dr = w * zz - w1 * z1 - x * y
        di = w * z1 + w1 * zz - x1 * y
        if ( dr == real(0.0E+00,kind=8) .and. di == real(0.0E+00,kind=8) ) dr = epsb
        go to 775

782     continue

        b(i+1,na) = t1
        b(i+1,en) = t2
        isw = 1
        if ( abs ( y ) > abs ( w ) + abs ( w1 ) ) go to 785
        tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
        ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
        go to 773

785     continue

        t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en) ) / y
        t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na) ) / y

787     continue

        b(i,na) = t1
        b(i,en) = t2

790     continue

     end do
!
!  End complex vector.
!
795   continue

      isw = 3 - isw

800   continue

  end do
!
!  End back substitution.
!  Transform to original coordinate system.
!
  do jj = 1, n

     j = n + 1 - jj

     do i = 1, n

        zz = real(0.0E+00,kind=8)

        do k = 1, j
          zz = zz + z(i,k) * b(k,j)
        end do

        z(i,j) = zz

      end do

  end do
!
!  Normalize so that modulus of largest component of each vector is 1.
!  (ISW is 1 initially from before).
!
  do j = 1, n

     d = real(0.0E+00,kind=8)
     if ( isw == 2 ) go to 920
     if ( alfi(j) /= real(0.0E+00,kind=8) ) go to 945

     do i = 1, n
       d = max ( d, abs ( z(i,j) ) )
     end do

     z(1:n,j) = z(1:n,j) / d

     go to 950

920  continue

     do i = 1, n
       r = abs ( z(i,j-1) ) + abs ( z(i,j) )
       if ( r /= real(0.0E+00,kind=8) ) then
         r = r * sqrt ( ( z(i,j-1) / r )**2 + ( z(i,j) / r )**2 )
       end if
       if ( r > d ) d = r
     end do

     z(1:n,j-1) = z(1:n,j-1) / d
     z(1:n,j) = z(1:n,j) / d

945  continue

     isw = 3 - isw

950  continue

  end do

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real(kind=8) values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real(kind=8) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real(kind=8) x
  real(kind=8) y
  real(kind=8) z
!
  z = x
  x = y
  y = z

  return
end
subroutine ratqr ( n, eps1, d, e, e2, m, w, ind, bd, type, idef, ierr )
!
!*******************************************************************************
!
!! RATQR computes selected eigenvalues of a real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the algebraically smallest or largest
!    eigenvalues of a symmetric tridiagonal matrix by the
!    rational QR method with Newton corrections.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) EPS1.  On input, a theoretical absolute error tolerance 
!    for the computed eigenvalues.  If the input EPS1 is non-positive,
!    or indeed smaller than its default value, it is reset at each iteration
!    to the respective default value, namely, the product of the relative 
!    machine precision and the magnitude of the current eigenvalue iterate.
!    The theoretical absolute error in the K-th eigenvalue is usually not 
!    greater than K times EPS1.  On output, EPS1 is unaltered unless it has 
!    been reset to its (last) default value.
!
!    Input, real(kind=8) D(N), the diagonal elements of the input matrix.
!
!    Input, real(kind=8) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, real(kind=8) E2(N).  On input, E2(2:N-1) contains the squares of 
!    the corresponding elements of E, and E2(1) is arbitrary.  On output,
!    elements of E2 corresponding to elements of E regarded as negligible 
!    have been replaced by zero, causing the matrix to split into a direct 
!    sum of submatrices.  E2(1) is set to real(0.0E+00,kind=8) if the smallest eigenvalues 
!    have been found, and to real(2.0E+00,kind=8) if the largest eigenvalues have been
!    found.  E2 is otherwise unaltered (unless overwritten by BD).
!
!    Input, integer M, the number of eigenvalues to be found.
!
!    Output, real(kind=8) W(M), the M algebraically smallest eigenvalues in
!    ascending order, or the M largest eigenvalues in descending order.  
!    If an error exit is made because of an incorrect specification of IDEF, 
!    no eigenvalues are found.  If the Newton iterates for a particular
!    eigenvalue are not monotone, the best estimate obtained is returned 
!    and IERR is set.  W may coincide with D.
!
!    Outpt, integer IND(N), contains in its first M positions the submatrix 
!    indices associated with the corresponding eigenvalues in W:
!    1 for eigenvalues belonging to the first submatrix from the top, 2 for 
!    those belonging to the second submatrix, and so on.
!
!    Output, real(kind=8) BD(N), contains refined bounds for the theoretical errors 
!    of the corresponding eigenvalues in W.  These bounds are usually
!    within the tolerance specified by EPS1.  BD may coincide with E2.
!
!    Input, integer IDEF, should be set to 1 if the input matrix is known to be
!    positive definite, to -1 if the input matrix is known to be negative 
!    definite, and to 0 otherwise.
!
!    Input, logical TYPE, should be set to TRUE if the smallest eigenvalues
!    are to be found, and to FALSE if the largest eigenvalues are to be found.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    6*N+1, if IDEF is set to 1 and TYPE to .true. when the matrix is not 
!      positive definite, or if IDEF is set to -1 and TYPE to .false.
!      when the matrix is not negative definite,
!    5*N+K, if successive iterates to the K-th eigenvalue are not monotone 
!      increasing, where K refers to the last such occurrence.
!
  implicit none
!
  integer n
!
  real(kind=8) bd(n)
  real(kind=8) d(n)
  real(kind=8) delta
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) ep
  real(kind=8) eps1
  real(kind=8) err
  real(kind=8) f
  integer i
  integer idef
  integer ierr
  integer ii
  integer ind(n)
  integer j
  integer jdef
  integer jj
  integer k
  integer m
  real(kind=8) p
  real(kind=8) q
  real(kind=8) qp
  real(kind=8) r
  real(kind=8) s
  real(kind=8) tot
  logical type
  real(kind=8) w(n)
!
  ierr = 0
  jdef = idef
  w(1:n) = d(1:n)

  if ( .not. type ) then
    j = 1
    go to 400
  end if

40 continue

  err = real(0.0E+00,kind=8)
  s = real(0.0E+00,kind=8)
!
!  Look for small sub-diagonal entries and define initial shift 
!  from lower Gerschgorin bound.
!
!  Copy E2 array into BD.
!
  tot = w(1)
  q = real(0.0E+00,kind=8)
  j = 0

  do i = 1, n

     p = q

     if ( i == 1 ) go to 60

     if ( p > ( abs ( d(i) ) + abs (  d(i-1) ) ) * epsilon ( p ) ) then
       go to 80
     end if

60   continue

     e2(i) = real(0.0E+00,kind=8)

80   continue

     bd(i) = e2(i)
!
!  Count also if element of E2 has underflowed.
!
     if ( e2(i) == real(0.0E+00,kind=8) ) j = j + 1
     ind(i) = j
     q = real(0.0E+00,kind=8)
     if ( i /= n ) q = abs ( e(i+1) )
     tot = min ( w(i)-p-q, tot )

  end do

  if ( jdef == 1 .and. tot < real(0.0E+00,kind=8) ) then
    go to 140
  end if

  w(1:n) = w(1:n) - tot

  go to 160

140 continue

  tot = real(0.0E+00,kind=8)

160 continue

  do k = 1, m
!
!  Next QR transformation.
!
180  continue

     tot = tot + s
     delta = w(n) - s
     i = n
     f = abs ( tot ) * epsilon ( f )
     if ( eps1 < f ) eps1 = f
     if ( delta > eps1 ) go to 190
     if ( delta < (-eps1) ) go to 1000
     go to 300
!
!  Replace small sub-diagonal squares by zero to reduce the incidence of 
!  underflows.
!
190  continue

     do j = k+1, n
       if ( bd(j) <= ( abs (  w(j) + w(j-1) ) * epsilon ( bd(j) ) ) ** 2 ) then
         bd(j) = real(0.0E+00,kind=8)
       end if
     end do

     f = bd(n) / delta
     qp = delta + f
     p = real(1.0E+00,kind=8)

     do ii = 1, n-k

       i = n - ii
       q = w(i) - s - f
       r = q / qp
       p = p * r + real(1.0E+00,kind=8)
       ep = f * r
       w(i+1) = qp + ep
       delta = q - ep

       if ( delta > eps1 ) go to 220
       if ( delta < (-eps1) ) go to 1000
       go to 300

220    continue

       f = bd(i) / q
       qp = delta + f
       bd(i+1) = qp * ep

     end do

     w(k) = qp
     s = qp / p

     if ( tot + s > tot ) go to 180
!
!  Set error: irregular end of iteration.
!  Deflate minimum diagonal element.
!
     ierr = 5 * n + k
     s = real(0.0E+00,kind=8)
     delta = qp

     do j = k, n
       if ( w(j) <= delta ) then
         i = j
         delta = w(j)
       end if
     end do
!
!  Convergence.
!
300  continue

     if ( i < n ) bd(i+1) = bd(i) * f / qp
     ii = ind(i)

     do jj = 1, i-k
       j = i - jj
       w(j+1) = w(j) - s
       bd(j+1) = bd(j)
       ind(j+1) = ind(j)
     end do

     w(k) = tot
     err = err + abs ( delta)
     bd(k) = err
     ind(k) = ii

  end do

  if ( type ) then
    return
  end if

  f = bd(1)
  e2(1) = real(2.0E+00,kind=8)
  bd(1) = f
  j = 2
!
!  Negate elements of W for largest values.
!
400 continue

  w(1:n) = - w(1:n)
  jdef = -jdef

  if ( j == 1 ) then
    go to 40
  end if

  return
!
!  Set error: IDEF specified incorrectly.
!
 1000 continue

  ierr = 6 * n + 1
  return
end
subroutine rebak ( nm, n, b, dl, m, z )
!
!*******************************************************************************
!
!! REBAK determines eigenvectors by undoing the REDUC transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a generalized
!    symmetric eigensystem by back transforming those of the
!    derived symmetric matrix determined by REDUC.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) B(NM,N), contains information about the similarity 
!    transformation (Cholesky decomposition) used in the reduction by REDUC
!    in its strict lower triangle.
!
!    Input, real(kind=8) DL(N), further information about the transformation.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) Z(NM,M).  On input, the eigenvectors to be back
!    transformed in its first M columns.  On output, the transformed 
!    eigenvectors.
!
  implicit none
!
  integer m
  integer nm
  integer n
!
  real(kind=8) b(nm,n)
  real(kind=8) dl(n)
  integer i
  integer j
  real(kind=8) z(nm,m)
!
  do j = 1, m
    do i = n, 1, -1
      z(i,j) = z(i,j) - dot_product ( b(i+1:n,i), z(i+1:n,j) ) / dl(i)
    end do
  end do

  return
end
subroutine rebakb ( nm, n, b, dl, m, z )
!
!*******************************************************************************
!
!! REBAKB determines eigenvectors by undoing the REDUC2 transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a generalized
!    symmetric eigensystem by back transforming those of the
!    derived symmetric matrix determined by REDUC2.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) B(NM,N), contains information about the similarity 
!    transformation (Cholesky decomposition) used in the reduction by REDUC2
!    in its strict lower triangle.
!
!    Input, real(kind=8) DL(N), further information about the transformation.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) Z(NM,M).  On input, the eigenvectors to be back
!    transformed in its first M columns.  On output, the transformed 
!    eigenvectors.
!
  implicit none
!
  integer m
  integer nm
  integer n
!
  real(kind=8) b(nm,n)
  real(kind=8) dl(n)
  integer i
  integer j
  real(kind=8) z(nm,m)
!
  do j = 1, m

    do i = n, 1, -1

      z(i,j) = dl(i) * z(i,j) + dot_product ( b(i,1:i-1), z(1:i-1,j) )

    end do

  end do

  return
end
subroutine reduc ( nm, n, a, b, dl, ierr )
!
!*******************************************************************************
!
!! REDUC reduces the eigenvalue problem A*x=lambda*B*x to A*x=lambda*x.
!
!
!  Discussion:
!
!    This subroutine reduces the generalized symmetric eigenproblem
!    ax=(lambda)bx, where B is positive definite, to the standard
!    symmetric eigenproblem using the Cholesky factorization of B.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices A and B.  If the Cholesky
!    factor L of B is already available, N should be prefixed with a minus sign.
!
!    Input/output, real(kind=8) A(NM,N).  On input, A contains a real(kind=8) symmetric matrix.
!    Only the full upper triangle of the matrix need be supplied.
!    On output, A contains in its full lower triangle the full lower triangle
!    of the symmetric matrix derived from the reduction to the
!    standard form.  The strict upper triangle of a is unaltered.
!
!    Input/output, real(kind=8) B(NM,N).  On input, the real(kind=8) symmetric input matrix.  
!    Only the full upper triangle of the matrix need be supplied.  If
!    N is negative, the strict lower triangle of B contains, instead, the 
!    strict lower triangle of its Cholesky factor L.  In any case, on output,
!    B contains in its strict lower triangle the strict lower triangle of 
!    its Cholesky factor L.  The full upper triangle of B is unaltered.
!
!    Input/output, real(kind=8) DL(N).  If N is negative, then the DL contains
!    the diagonal elements of L on input.  In any case, DL will contain
!    the diagonal elements of L on output,
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    7*N+1, if B is not positive definite.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) b(nm,n)
  real(kind=8) dl(n)
  integer i
  integer ierr
  integer j
  integer k
  integer nn
  real(kind=8) x
  real(kind=8) y
!
  ierr = 0
  nn = abs ( n )
!
!  Form L in the arrays B and DL.
!
  do i = 1, n

     do j = i, n

        x = b(i,j)

        do k = 1, i - 1
          x = x - b(i,k) * b(j,k)
        end do

        if ( j == i ) then

          if ( x <= real(0.0E+00,kind=8) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'REDUC - Fatal error!'
            write ( *, '(a)' ) '  The matrix is not positive definite.'
            ierr = 7 * n + 1
            return
          end if

          y = sqrt ( x )
          dl(i) = y
        else
          b(j,i) = x / y
        end if

    end do

  end do
!
!  Form the transpose of the upper triangle of INV(L)*A
!  in the lower triangle of the array A.
!
  do i = 1, nn

     y = dl(i)

     do j = i, nn

        x = a(i,j)

        do k = 1, i - 1
          x = x - b(i,k) * a(j,k)
        end do

        a(j,i) = x / y

      end do

  end do
!
!  Pre-multiply by INV(L) and overwrite.
!
  do j = 1, nn

     do i = j, nn

        x = a(i,j)

        do k = j, i-1
          x = x - a(k,j) * b(i,k)
        end do

        do k = 1, j-1
          x = x - a(j,k) * b(i,k)
        end do

        a(i,j) = x / dl(i)

    end do

  end do

  return
end
subroutine reduc2 ( nm, n, a, b, dl, ierr )
!
!*******************************************************************************
!
!! REDUC2 reduces the eigenvalue problem A*B*x=lamdba*x to A*x=lambda*x.
!
!
!  Discussion:
!
!    This subroutine reduces the generalized symmetric eigenproblems
!    abx=(lambda)x or bay=(lambda)y, where B is positive definite,
!    to the standard symmetric eigenproblem using the Cholesky
!    factorization of B.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices A and B.  If the Cholesky
!    factor L of B is already available, N should be prefixed with a minus sign.
!
!    Input/output, real(kind=8) A(NM,N).  On input, A contains a real(kind=8) symmetric matrix.
!    Only the full upper triangle of the matrix need be supplied.
!    On output, A contains in its full lower triangle the full lower triangle
!    of the symmetric matrix derived from the reduction to the
!    standard form.  The strict upper triangle of a is unaltered.
!
!    Input/output, real(kind=8) B(NM,N).  On input, the real(kind=8) symmetric input matrix.  
!    Only the full upper triangle of the matrix need be supplied.  If
!    N is negative, the strict lower triangle of B contains, instead, the 
!    strict lower triangle of its Cholesky factor L.  In any case, on output,
!    B contains in its strict lower triangle the strict lower triangle of 
!    its Cholesky factor L.  The full upper triangle of B is unaltered.
!
!    Input/output, real(kind=8) DL(N).  If N is negative, then the DL contains
!    the diagonal elements of L on input.  In any case, DL will contain
!    the diagonal elements of L on output,
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    7*N+1, if B is not positive definite.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) b(nm,n)
  real(kind=8) dl(n)
  integer i
  integer ierr
  integer j
  integer k
  integer nn
  real(kind=8) x
  real(kind=8) y
!
  ierr = 0
  nn = abs ( n )
!
!  Form L in the arrays B and DL.
!
  do i = 1, n

     do j = i, n

        x = b(i,j)

        do k = 1, i - 1
          x = x - b(i,k) * b(j,k)
        end do

        if ( j == i ) then

          if ( x <= real(0.0E+00,kind=8) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'REDUC2 - Fatal error!'
            write ( *, '(a)' ) '  The matrix is not positive definite.'
            ierr = 7 * n + 1
            return
          end if

          y = sqrt ( x )
          dl(i) = y

        else

          b(j,i) = x / y

        end if

    end do

  end do
!
!  Form the lower triangle of A*L in the lower triangle of A.
!
  do i = 1, nn

     do j = 1, i

        x = a(j,i) * dl(j)

        do k = j+1, i
          x = x + a(k,i) * b(k,j)
        end do

        do k = i+1, nn
          x = x + a(i,k) * b(k,j)
        end do

        a(i,j) = x

     end do

  end do
!
!  Pre-multiply by L' and overwrite.
!
  do i = 1, nn

    y = dl(i)

    do j = 1, i

      x = y * a(i,j)

      do k = i+1, nn
        x = x + a(k,j) * b(k,i)
      end do

      a(i,j) = x

    end do

  end do

  return
end
subroutine rg ( nm, n, a, wr, wi, matz, z, ierr )
!
!*******************************************************************************
!
!! RG computes eigenvalues and eigenvectors of a real(kind=8) general matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real(kind=8) general matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) A(NM,N), the real(kind=8) general matrix.  On output,
!    A has been overwritten.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) WR(N), WI(N), the real(kind=8) and imaginary parts, respectively, 
!    of the eigenvalues.  Complex conjugate pairs of eigenvalues appear 
!    consecutively with the eigenvalue having the positive imaginary part first.
!
!    Output, real(kind=8) Z(NM,N), contains the real(kind=8) and imaginary parts of the 
!    eigenvectors if MATZ is not zero.  If the J-th eigenvalue is real(kind=8), the
!    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is 
!    complex with positive imaginary part, the J-th and (J+1)-th columns of 
!    Z contain the real(kind=8) and imaginary parts of its eigenvector.  The 
!    conjugate of this vector is the eigenvector for the conjugate eigenvalue.
!
!    Output, integer IERR, an error completion code described in the 
!    documentation for HQR and HQR2.  The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) fv1(n)
  integer ierr
  integer is1
  integer is2
  integer iv1(n)
  integer matz
  real(kind=8) wi(n)
  real(kind=8) wr(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call balanc ( nm, n, a, is1, is2, fv1 )

  call elmhes ( nm, n, is1, is2, a, iv1 )

  if ( matz == 0 ) then

    call hqr ( nm, n, is1, is2, a, wr, wi, ierr )

    if ( ierr /= 0 ) then
      return
    end if

  else

    call eltran ( nm, n, is1, is2, a, iv1, z )

    call hqr2 ( nm, n, is1, is2, a, wr, wi, z, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call balbak ( nm, n, is1, is2, fv1, n, z )

  end if

  return
end
subroutine rgg ( nm, n, a, b, alfr, alfi, beta, matz, z, ierr )
!
!*******************************************************************************
!
!! RGG computes eigenvalues/vectors for the generalized problem A*x = lambda*B*x.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real(kind=8) general generalized eigenproblem  
!
!      A * x = lambda * B * x.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices A and B.
!
!    Input/output, real(kind=8) A(NM,N), B(NM,N), the two real(kind=8) general matrices.
!    On output, A and B have been overwritten.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) ALFR(N), ALFI(N), the real(kind=8) and imaginary parts,
!    respectively, of the numerators of the eigenvalues.
!
!    Output, real(kind=8) BETA(N), the denominators of the eigenvalues,
!    which are thus given by the ratios (ALFR + I * ALFI ) / BETA.
!    Complex conjugate pairs of eigenvalues appear consecutively
!    with the eigenvalue having the positive imaginary part first.
!
!    Output, real(kind=8) Z(NM,N), contains the real(kind=8) and imaginary parts of the 
!    eigenvectors if MATZ is not zero.  If the J-th eigenvalue is real(kind=8), the
!    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is 
!    complex with positive imaginary part, the J-th and (J+1)-th columns of 
!    Z contain the real(kind=8) and imaginary parts of its eigenvector.  The 
!    conjugate of this vector is the eigenvector for the conjugate eigenvalue.
!
!    Output, integer IERR, is set equal to an error completion code 
!    described in the documentation for QZIT.  The normal completion 
!    code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) alfi(n)
  real(kind=8) alfr(n)
  real(kind=8) b(nm,n)
  real(kind=8) beta(n)
  real(kind=8) eps1
  integer ierr
  integer matz
  logical tf
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  eps1 = real(0.0E+00,kind=8)

  if ( matz == 0 ) then
    tf = .false.
  else
    tf = .true.
  end if

  call qzhes ( nm, n, a, b, tf, z )

  call qzit ( nm, n, a, b, eps1, tf, z, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  call qzval ( nm, n, a, b, alfr, alfi, beta, tf, z )

  if ( matz /= 0 ) then
    call qzvec ( nm, n, a, b, alfr, alfi, beta, z )
  end if

  return
end
subroutine rmat_ident ( lda, n, a )
!
!*******************************************************************************
!
!! RMAT_IDENT sets the square matrix A to the identity.
!
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of A.
!
!    Output, real(kind=8) A(LDA,N), the matrix which has been
!    set to the identity.
!
  implicit none
!
  integer lda
  integer n
!
  real(kind=8) a(lda,n)
  integer i
  integer j
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = real(1.0E+00,kind=8)
      else
        a(i,j) = real(0.0E+00,kind=8)
      end if
    end do
  end do

  return
end
subroutine rmat_print ( m, n, a, title )
!
!*******************************************************************************
!
!! RMAT_PRINT prints a real(kind=8) matrix, with an optional title.
!
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real(kind=8) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer m
  integer n
!
  real(kind=8) a(m,n)
  integer i
  integer j
  integer jhi
  integer jlo
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 5

    jhi = min ( jlo + 4, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )

    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do

  end do

  return
end
subroutine rs ( nm, n, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RS computes eigenvalues and eigenvectors of real(kind=8) symmetric matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real(kind=8) symmetric matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) A(NM,N), the real(kind=8) symmetric matrix.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set equal to an error
!    completion code described in the documentation for TQLRAT and TQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  integer ierr
  integer matz
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

  end if

  return
end
subroutine rsb ( nm, n, mb, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSB computes eigenvalues and eigenvectors of a real(kind=8) symmetric band matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real(kind=8) symmetric band matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MB, the half band width of the matrix, defined as the
!    number of adjacent diagonals, including the principal diagonal, required
!    to specify the non-zero portion of the lower triangle of the matrix.
!
!    Input, real(kind=8) A(NM,MB), contains the lower triangle of the real(kind=8) symmetric
!    band matrix.  Its lowest subdiagonal is stored in the last N+1-MB 
!    positions of the first column, its next subdiagonal in the last 
!    N+2-MB positions of the second column, further subdiagonals similarly, 
!    and finally its principal diagonal in the N positions of the last 
!    column.  Contents of storages not part of the matrix are arbitrary.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer mb
  integer n
  integer nm
!
  real(kind=8) a(nm,mb)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  integer ierr
  integer matz
  logical tf
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( mb <= 0 ) then
    ierr = 12 * n
    return
  end if

  if ( mb > n ) then
    ierr = 12 * n
    return
  end if

  if ( matz == 0 ) then
 
    tf = .false.

    call bandr ( nm, n, mb, a, w, fv1, fv2, tf, z )

    call tqlrat ( n, w, fv2, ierr )

  else

    tf = .true.

    call bandr ( nm, n, mb, a, w, fv1, fv1, tf, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

  end if

  return
end
subroutine rsg ( nm, n, a, b, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSG computes eigenvalues/vectors, A*x=lambda*B*x, A symmetric, B pos-def.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real(kind=8) symmetric generalized eigenproblem  ax = (lambda)bx.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Modified:
!
!    12 June 2000
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices A and B.
!
!    Input, real(kind=8) A(NM,N), contains a real(kind=8) symmetric matrix.
!
!    Input, real(kind=8) B(NM,N), contains a positive definite real(kind=8) symmetric matrix.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) b(nm,n)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  integer ierr
  integer matz
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RSG - Fatal error!'
    write ( *, '(a)' ) '  N > NM.'
    return
  end if

  call reduc ( nm, n, a, b, fv2, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RSG - Fatal error!'
    write ( *, '(a)' ) '  Error return from REDUC.'
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSG - Warning!'
      write ( *, '(a)' ) '  Error return from TQLRAT!'
      return
    end if

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSG - Fatal error!'
      write ( *, '(a)' ) '  Error return from TQL2!'
      return
    end if

    call rebak ( nm, n, b, fv2, n, z )

  end if

  return
end
subroutine rsgab ( nm, n, a, b, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSGAB computes eigenvalues/vectors, A*B*x=lambda*x, A symmetric, B pos-def.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real(kind=8) symmetric generalized eigenproblem  abx = (lambda)x.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices A and B.
!
!    Input, real(kind=8) A(NM,N), contains a real(kind=8) symmetric matrix.
!
!    Input, real(kind=8) B(NM,N), contains a positive definite real(kind=8) symmetric matrix.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) b(nm,n)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  integer ierr
  integer matz
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call reduc2 ( nm, n, a, b, fv2, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call rebak ( nm, n, b, fv2, n, z )

  end if

  return
end
subroutine rsgba ( nm, n, a, b, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSGBA computes eigenvalues/vectors, B*A*x=lambda*x, A symmetric, B pos-def.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real(kind=8) symmetric generalized eigenproblem:
! 
!      B * A * x = lambda * x
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrices A and B.
!
!    Input, real(kind=8) A(NM,N), a real(kind=8) symmetric matrix.
!
!    Input, real(kind=8) B(NM,N), a positive definite symmetric matrix.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) b(nm,n)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  integer ierr
  integer matz
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call reduc2 ( nm, n, a, b, fv2, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call rebakb ( nm, n, b, fv2, n, z )

  end if

  return
end
subroutine rsm ( nm, n, a, w, m, z, ierr )
!
!*******************************************************************************
!
!! RSM computes eigenvalues, some eigenvectors, real(kind=8) symmetric matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find all of the eigenvalues and some of the eigenvectors
!    of a real(kind=8) symmetric matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) A(NM,N), the symmetric matrix.
!
!    Input, integer M, specifies the number of eigenvectors to compute.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,M), contains the orthonormal eigenvectors associated
!    with the first M eigenvalues.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for TQLRAT, IMTQLV and
!    TINVIT.  The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
  integer m
!
  real(kind=8) a(nm,n)
  real(kind=8) fwork1(n)
  real(kind=8) fwork2(n)
  real(kind=8) fwork3(n)
  integer ierr
  integer iwork(n)
  integer k1
  integer k2
  integer k3
  integer k4
  real(kind=8) w(n)
  real(kind=8) z(nm,m)
!
  if ( n > nm .or. m > nm ) then
    ierr = 10 * n
    return
  end if

  k1 = 1
  k2 = k1 + n
  k3 = k2 + n
  k4 = k3 + n

  if ( m <= 0 ) then

    call tred1 ( nm, n, a, w, fwork1, fwork2 )

    call tqlrat ( n, w, fwork2, ierr )

  else

    call tred1 ( nm, n, a, fwork1, fwork2, fwork3 )

    call imtqlv ( n, fwork1, fwork2, fwork3, w, iwork, ierr )

    call tinvit ( nm, n, fwork1, fwork2, fwork3, m, w, iwork, z, ierr )

    call trbak1 ( nm, n, a, fwork2, m, z )

  end if

  return
end
subroutine rsp ( nm, n, nv, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSP computes eigenvalues and eigenvectors of real(kind=8) symmetric packed matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of
!    subroutines from the eigensystem subroutine package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real(kind=8) symmetric packed matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NV, the dimension of the array A, which must be at
!    least (N*(N+1))/2.
!
!    Input, real(kind=8) A(NV), contains the lower triangle of the real(kind=8) symmetric
!    packed matrix stored row-wise.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
  integer nv
!
  real(kind=8) a(nv)
  real(kind=8) fv1(n)
  real(kind=8) fv2(n)
  integer ierr
  integer matz
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( ( n * ( n + 1 ) ) / 2 > nv ) then
    ierr = 20 * n
    return
  end if

  call tred3 ( n, nv, a, w, fv1, fv2 )

  if ( matz == 0 ) then

    call tqlrat ( n, w, fv2, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSP - Fatal error!'
      write ( *, '(a)' ) '  Error return from TQLRAT.'
      return
    end if

  else

    call rmat_ident ( nm, n, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSP - Fatal error!'
      write ( *, '(a)' ) '  Error return from TQL2.'
      return
    end if

    call trbak3 ( nm, n, nv, a, n, z )

  end if

  return
end
subroutine rspp ( nm, n, nv, a, w, matz, z, ierr, m, type )
!
!*******************************************************************************
!
!! RSPP computes some eigenvalues/vectors, real(kind=8) symmetric packed matrix.
!
!
!  Discussion:
!
!    This routine calls the appropriate routines for the following problem:
!
!    Given a symmetric matrix A, which is stored in a packed mode, find
!    the M smallest or largest eigenvalues, and corresponding eigenvectors.
!
!    The routine RSP returns all eigenvalues and eigenvectors.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the eigenvector array z.
!
!    Input, integer N, the order of A, the number of rows and columns in the
!    original matrix.
!
!    Input, integer NV, is the of the array A as specified in the calling 
!    program.  NV must not be less than N*(N+1)/2.
!
!    Input, real(kind=8) A((N*(N+1))/2), on input the lower triangle of the
!    real(kind=8) symmetric matrix, stored row-wise in the vector,
!    in the order A(1,1), / A(2,1), A(2,2), / A(3,1), A(3,2), A(3,3)/
!    and so on.
!
!    Output, real(kind=8) W(M), the eigenvalues requested.
!
!    Input, integer MATZ, is set to 0 if only eigenvalues are desired.  
!    Otherwise it is set to any non-zero integer for both eigenvalues 
!    and eigenvectors.
!
!    Output, real(kind=8) Z(NM,M), the eigenvectors.
!
!    Output, integer IERR, error flag from RATQR.  IERR=0 on normal return.
!    IERR nonzero, in this case, means that the algorithm broke
!    down while computing an eigenvalue.
!
!    Input, integer M, the number of eigenvalues to be found.
!
!    Input, logical TYPE, set to .true. if the smallest eigenvalues
!    are to be found, or .false. if the largest ones are sought.
!
  implicit none
!
  integer m
  integer n
  integer nm
  integer nv
!
  real(kind=8) a(nv)
  real(kind=8) bd(n)
  real(kind=8) eps1
  integer idef
  integer ierr
  integer iwork(n)
  integer matz
  logical type
  real(kind=8) w(m)
  real(kind=8) work1(n)
  real(kind=8) work2(n)
  real(kind=8) work3(n)
  real(kind=8) z(nm,m)
!
!  IDEF =
!    -1 if the matrix is known to be negative definite, 
!    +1 if the matrix is known to be positive definite, or
!    0 otherwise.
!
  idef = 0
!
!  Reduce to symmetric tridiagonal form.
!
  call tred3 ( n, nv, a, work1, work2, work3 )
!
!  Find the eigenvalues.
!
  eps1 = real(0.0E+00,kind=8)

  call ratqr ( n, eps1, work1, work2, work3, m, w, iwork, &
    bd, type, idef, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RSPP - Fatal error!'
    write ( *, '(a)' ) '  Error return from RATQR.'
    return
  end if
!
!  Find eigenvectors for the first M eigenvalues.
!
  if ( matz /= 0 ) then

    call tinvit ( nm, n, work1, work2, work3, m, w, iwork, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSPP - Fatal error!'
      write ( *, '(a)' ) '  Error return from TINVIT.'
      return
    end if
!
!  Reverse the transformation.
!
    call trbak3 ( nm, n, nv, a, m, z )

  end if

  return
end
subroutine rst ( nm, n, w, e, matz, z, ierr )
!
!*******************************************************************************
!
!! RST computes eigenvalues/vectors, real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of subroutines 
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real(kind=8) symmetric tridiagonal matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) W(N).  On input, the diagonal elements of the real(kind=8)
!    symmetric tridiagonal matrix.  On output, the eigenvalues in ascending 
!    order.
!
!    Input, real(kind=8) E(N), the subdiagonal elements of the matrix in
!    E(2:N).  E(1) is arbitrary.
!
!    Input, integer MATZ, is zero if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for IMTQL1 and IMTQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) e(n)
  integer ierr
  integer matz
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( matz == 0 ) then

    call imtql1 ( n, w, e, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RST - Fatal error!'
      write ( *, '(a)' ) '  Error return from IMTQL1.'
      return
    end if

  else

    call rmat_ident ( nm, n, z )

    call imtql2 ( nm, n, w, e, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RST - Fatal error!'
      write ( *, '(a)' ) '  Error return from IMTQL2.'
      return
    end if

  end if

  return
end
subroutine rt ( nm, n, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RT computes eigenvalues/vectors, real(kind=8) sign-symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of subroutines 
!    to find the eigenvalues and eigenvectors (if desired)
!    of a special real(kind=8) tridiagonal matrix.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) A(NM,N), contains the special real(kind=8) tridiagonal matrix in 
!    its first three columns.  The subdiagonal elements are stored
!    in the last N-1 positions of the first column, the diagonal elements 
!    in the second column, and the superdiagonal elements in the first N-1 
!    positions of the third column.  Elements A(1,1) and A(N,3) are arbitrary.
!
!    Input, integer MATZ, is 0 if only eigenvalues are desired, and nonzero
!    if both eigenvalues and eigenvectors are desired.
!
!    Output, real(kind=8) W(N), the eigenvalues in ascending order.
!
!    Output, real(kind=8) Z(NM,N), contains the eigenvectors, if MATZ is nonzero.
!
!    Output, integer IERR, is set to an error
!    completion code described in the documentation for IMTQL1 and IMTQL2.  
!    The normal completion code is zero.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) fv1(n)
  integer ierr
  integer matz
  real(kind=8) a(nm,3)
  real(kind=8) w(n)
  real(kind=8) z(nm,n)
!
  if ( n > nm ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RT - Fatal error!'
    write ( *, '(a)' ) '  N greater than NM.'
    ierr = 10 * n
    return
  end if

  if ( matz == 0 ) then

    call figi ( nm, n, a, w, fv1, fv1, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error return from FIGI.'
      return
    end if

    call imtql1 ( n, w, fv1, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error return from IMTQL1.'
      return
    end if

  else

    call figi2 ( nm, n, a, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error return from FIGI2.'
      return
    end if

    call imtql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error return from IMTQL2.'
      return
    end if

  end if

  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real(kind=8) vector, with an optional title.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real(kind=8) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  real(kind=8) a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine rvec2_print ( n, a1, a2, title )
!
!*******************************************************************************
!
!! RVEC2_PRINT prints a pair of real(kind=8) vectors, with an optional title.
!
!
!  Modified:
!
!    14 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real(kind=8) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  real(kind=8) a1(n)
  real(kind=8) a2(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,2g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine svd ( m, n, a, w, matu, u, matv, v, ierr )
!
!*******************************************************************************
!
!! SVD computes the singular value decomposition for a real(kind=8) matrix.
!
!
!  Discussion:
!
!    This subroutine determines the singular value decomposition
!
!      A = U * S * V' 
!
!    of a real(kind=8) M by N rectangular matrix.  Householder bidiagonalization 
!    and a variant of the QR algorithm are used.
!
!  Modified:
!
!    22 April 2002
!
!  Reference:
!
!    Golub and Reinsch,
!    Numerische Mathematik,
!    Volume 14, 1970, pages 403-420.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer M, the number of rows of A and U.
!
!    Input, integer N, the number of columns of A and U, and the order of V.
!
!    Input, real(kind=8) A(M,N), the M by N matrix to be decomposed.
!
!    Output, real(kind=8) W(N), the singular values of A.  These are the
!    diagonal elements of S.  They are unordered.  If an error exit is 
!    made, the singular values should be correct for indices 
!    IERR+1, IERR+2,..., N.
!
!    Input, logical MATU, should be set to TRUE if the U matrix in the
!    decomposition is desired, and to FALSE otherwise.
!
!    Output, real(kind=8) U(M,N), contains the matrix U, with orthogonal columns,
!    of the decomposition, if MATU has been set to TRUE.  Otherwise
!    U is used as a temporary array.  U may coincide with A.
!    If an error exit is made, the columns of U corresponding
!    to indices of correct singular values should be correct.
!
!    Input, logical MATV, should be set to TRUE if the V matrix in the
!    decomposition is desired, and to FALSE otherwise.
!
!    Output, real(kind=8) V(N,N), the orthogonal matrix V of the decomposition if
!    MATV has been set to TRUE.  Otherwise V is not referenced.
!    V may also coincide with A if U is not needed.  If an error
!    exit is made, the columns of V corresponding to indices of
!    correct singular values should be correct.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    K, if the K-th singular value has not been determined after 30 iterations.
!
  implicit none
!
  integer m
  integer n
!
  real(kind=8) a(m,n)
  real(kind=8) c
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  integer i
  integer ierr
  integer its
  integer i1
  integer j
  integer k
  integer kk
  integer k1
  integer l
  integer ll
  integer l1
!  integer m
  logical matu
  logical matv
  integer mn
  real(kind=8) pythag
  real(kind=8) rv1(n)
  real(kind=8) s
  real(kind=8) scale
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) u(m,n)
  real(kind=8) v(n,n)
  real(kind=8) w(n)
  real(kind=8) x
  real(kind=8) y
  real(kind=8) z
!
  ierr = 0
  u(1:m,1:n) = a(1:m,1:n)
!
!  Householder reduction to bidiagonal form.
!
  g = real(0.0E+00,kind=8)
  scale = real(0.0E+00,kind=8)
  x = real(0.0E+00,kind=8)

  do i = 1, n

    l = i + 1
    rv1(i) = scale * g
    g = real(0.0E+00,kind=8)
    s = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)

    if ( i <= m ) then

      scale = sum ( abs ( u(i:m,i) ) )

      if ( scale /= real(0.0E+00,kind=8) ) then

        u(i:m,i) = u(i:m,i) / scale

        s = sum ( u(i:m,i)**2 )

        f = u(i,i)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        u(i,i) = f - g

        if ( i /= n ) then

          do j = l, n
            s = dot_product ( u(i:m,i), u(i:m,j) )
            u(i:m,j) = u(i:m,j) + s * u(i:m,i) / h
          end do

        end if

        u(i:m,i) = scale * u(i:m,i)

      end if

    end if

    w(i) = scale * g
    g = real(0.0E+00,kind=8)
    s = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)

    if ( i <= m .and. i /= n ) then

      scale = sum ( abs ( u(i,l:n) ) )

      if ( scale /= real(0.0E+00,kind=8) ) then

        u(i,l:n) = u(i,l:n) / scale
        s = sum ( u(i,l:n)**2 )
        f = u(i,l)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        u(i,l) = f - g
        rv1(l:n) = u(i,l:n) / h

        if ( i /= m ) then

          do j = l, m

            s = dot_product ( u(j,l:n), u(i,l:n) )

            u(j,l:n) = u(j,l:n) + s * rv1(l:n)

          end do

        end if

        u(i,l:n) = scale * u(i,l:n)

      end if

    end if

    x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

  end do
!
!  Accumulation of right-hand transformations.
!
  if ( matv ) then

    do i = n, 1, -1

      if ( i /= n ) then

         if ( g /= real(0.0E+00,kind=8) ) then

          v(l:n,i) = ( u(i,l:n) / u(i,l) ) / g

          do j = l, n

            s = dot_product ( u(i,l:n), v(l:n,j) )

            v(l:n,j) = v(l:n,j) + s * v(l:n,i)

          end do

        end if

        v(i,l:n) = real(0.0E+00,kind=8)
        v(l:n,i) = real(0.0E+00,kind=8)

      end if

      v(i,i) = real(1.0E+00,kind=8)
      g = rv1(i)
      l = i

    end do

  end if
!
!  Accumulation of left-hand transformations.
!
  if ( matu ) then

    mn = min ( m, n )

    do i = min ( m, n ), 1, -1

      l = i + 1
      g = w(i)

      if ( i /= n ) then
        u(i,l:n) = real(0.0E+00,kind=8)
      end if

      if ( g /= real(0.0E+00,kind=8) ) then

        if ( i /= mn ) then

          do j = l, n
            s = dot_product ( u(l:m,i), u(l:m,j) )
            f = ( s / u(i,i) ) / g
            u(i:m,j) = u(i:m,j) + f * u(i:m,i)
          end do

        end if

        u(i:m,i) = u(i:m,i) / g

      else

        u(i:m,i) = real(0.0E+00,kind=8)

      end if

      u(i,i) = u(i,i) + real(1.0E+00,kind=8)

    end do

  end if
!
!  Diagonalization of the bidiagonal form.
!
  tst1 = x

  do kk = 1, n

     k1 = n - kk
     k = k1 + 1
     its = 0
!
!  Test for splitting.
!
520  continue

     do ll = 1, k

       l1 = k - ll
       l = l1 + 1
       tst2 = tst1 + abs ( rv1(l) )

       if ( tst2 == tst1 ) then
         go to 565
       end if

       tst2 = tst1 + abs ( w(l1) )

       if ( tst2 == tst1 ) then
         exit
       end if

     end do
!
!  Cancellation of rv1(l) if L greater than 1.
!
     c = real(0.0E+00,kind=8)
     s = real(1.0E+00,kind=8)

     do i = l, k

       f = s * rv1(i)
       rv1(i) = c * rv1(i)
       tst2 = tst1 + abs ( f )

       if ( tst2 == tst1 ) then
         go to 565
       end if

       g = w(i)
       h = pythag ( f, g )
       w(i) = h
       c = g / h
       s = -f / h

       if ( matu ) then

         do j = 1, m
           y = u(j,l1)
           z = u(j,i)
           u(j,l1) = y * c + z * s
           u(j,i) = -y * s + z * c
         end do

       end if

    end do
!
!  Test for convergence.
!
565 continue
 
    z = w(k)

    if ( l == k ) go to 650
!
!  Shift from bottom 2 by 2 minor.
!
    if ( its >= 30 ) then
      ierr = k
      return
    end if

    its = its + 1
    x = w(l)
    y = w(k1)
    g = rv1(k1)
    h = rv1(k)
    f = real(0.5E+00,kind=8) * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y )
    g = pythag ( f, real(1.0E+00,kind=8) )
    f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + sign ( g, f ) ) - h)
!
!  Next QR transformation.
!
    c = real(1.0E+00,kind=8)
    s = real(1.0E+00,kind=8)

    do i1 = l, k1

      i = i1 + 1
      g = rv1(i)
      y = w(i)
      h = s * g
      g = c * g
      z = pythag ( f, h )
      rv1(i1) = z
      c = f / z
      s = h / z
      f = x * c + g * s
      g = -x * s + g * c
      h = y * s
      y = y * c

      if ( matv ) then

        do j = 1, n
          x = v(j,i1)
          z = v(j,i)
          v(j,i1) = x * c + z * s
          v(j,i) = -x * s + z * c
        end do

      end if

      z = pythag ( f, h )
      w(i1) = z
!
!  Rotation can be arbitrary if Z is zero.
!
      if ( z /= real(0.0E+00,kind=8) ) then
        c = f / z
        s = h / z
      end if

      f = c * g + s * y
      x = -s * g + c * y

      if ( matu ) then

        do j = 1, m
          y = u(j,i1)
          z = u(j,i)
          u(j,i1) = y * c + z * s
          u(j,i) = -y * s + z * c
        end do

      end if

    end do

    rv1(l) = real(0.0E+00,kind=8)
    rv1(k) = f
    w(k) = x
    go to 520
!
!  Convergence.
!
650 continue

    if ( z <= real(0.0E+00,kind=8) ) then

      w(k) = - z

      if ( matv ) then
        v(1:n,k) = - v(1:n,k)
      end if

    end if

  end do

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tinvit ( nm, n, d, e, e2, m, w, ind, z, ierr )
!
!*******************************************************************************
!
!! TINVIT computes eigenvectors from eigenvalues, real(kind=8) tridiagonal symmetric.
!
!
!  Discussion:
!
!    This subroutine finds those eigenvectors of a tridiagonal
!    symmetric matrix corresponding to specified eigenvalues,
!    using inverse iteration.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) D(N), the diagonal elements of the matrix.
!
!    Input, real(kind=8) E(N), contains the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input, real(kind=8) E2(N), contains the squares of the corresponding elements 
!    of E, with zeros corresponding to negligible elements of E.
!    E(I) is considered negligible if it is not larger than the product of 
!    the relative machine precision and the sum of the magnitudes of D(I) 
!    and D(I-1).  E2(1) must contain real(0.0E+00,kind=8) if the eigenvalues are in 
!    ascending order, or real(2.0E+00,kind=8) if the eigenvalues are in descending order.  
!    If BISECT, TRIDIB, or IMTQLV has been used to find the eigenvalues,
!    their output E2 array is exactly what is expected here.
!
!    Input, integer M, the number of specified eigenvalues.
!
!    Input, real(kind=8) W(M), the eigenvalues.
!
!    Input, integer IND(M), the submatrix indices associated with the 
!    corresponding eigenvalues in W: 1 for eigenvalues belonging to the 
!    first submatrix from the top, 2 for those belonging to the second 
!    submatrix, and so on.
!
!    Output, real(kind=8) Z(NM,M), the associated set of orthonormal eigenvectors.
!    Any vector which fails to converge is set to zero.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    -R, if the eigenvector corresponding to the R-th eigenvalue fails to 
!      converge in 5 iterations.
!
  implicit none
!
  integer m
  integer n
  integer nm
!
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) eps2
  real(kind=8) eps3
  real(kind=8) eps4
  integer group
  integer i
  integer ierr
  integer ii
  integer ind(m)
  integer ip
  integer its
  integer j
  integer jj
  real(kind=8) norm
  real(kind=8) order
  integer p
  real(kind=8) pythag
  integer q
  integer r
  real(kind=8) rv1(n)
  real(kind=8) rv2(n)
  real(kind=8) rv3(n)
  real(kind=8) rv4(n)
  real(kind=8) rv6(n)
  integer s
  integer tag
  real(kind=8) u
  real(kind=8) uk
  real(kind=8) v
  real(kind=8) w(m)
  real(kind=8) x0
  real(kind=8) x1
  real(kind=8) xu
  real(kind=8) z(nm,m)
!
  ierr = 0

  if ( m == 0 ) then
    return
  end if

  u = real(0.0E+00,kind=8)
  x0 = real(0.0E+00,kind=8)

  tag = 0
  order = real(1.0E+00,kind=8) - e2(1)
  q = 0
!
!  Establish and process next submatrix.
!
100 continue

  p = q + 1

  do q = p, n
    if ( q == n ) then
      exit
    end if
    if ( e2(q+1) == real(0.0E+00,kind=8) ) then
      exit
    end if
  end do
!
!  Find vectors by inverse iteration.
!
140 continue

  tag = tag + 1
  s = 0

  do r = 1, m

     if ( ind(r) /= tag ) go to 920

     its = 1
     x1 = w(r)

     if ( s /= 0 ) go to 510
!
!  Check for isolated root.
!
     xu = real(1.0E+00,kind=8)

     if ( p == q ) then
       rv6(p) = real(1.0E+00,kind=8)
       go to 870
     end if

     norm = abs ( d(p) )
     ip = p + 1

     do i = p+1, q
       norm = max ( norm, abs ( d(i) ) + abs ( e(i) ) )
     end do
!
!  EPS2 is the criterion for grouping,
!  EPS3 replaces zero pivots and equal roots are modified by EPS3,
!  EPS4 is taken very small to avoid overflow.
!
     eps2 = 0.001E+00 * norm
     eps3 = abs ( norm ) * epsilon ( eps3 )
     uk = q - p + 1
     eps4 = uk * eps3
     uk = eps4 / sqrt ( uk )
     s = p

505 continue

     group = 0
     go to 520
!
!  Look for close or coincident roots.
!
510  continue

     if ( abs ( x1 - x0 ) >= eps2 ) go to 505

     group = group + 1

     if ( order * (x1 - x0) <= real(0.0E+00,kind=8) ) then
       x1 = x0 + order * eps3
     end if
!
!  Elimination with interchanges and initialization of vector.
!
520  continue

     v = real(0.0E+00,kind=8)

     do i = p, q

        rv6(i) = uk

        if ( i == p ) go to 560

        if ( abs ( e(i) ) < abs ( u ) ) go to 540

        xu = u / e(i)
        rv4(i) = xu
        rv1(i-1) = e(i)
        rv2(i-1) = d(i) - x1
        rv3(i-1) = real(0.0E+00,kind=8)
        if ( i /= q ) rv3(i-1) = e(i+1)
        u = v - xu * rv2(i-1)
        v = - xu * rv3(i-1)
        go to 580

540     continue

        xu = e(i) / u
        rv4(i) = xu
        rv1(i-1) = u
        rv2(i-1) = v
        rv3(i-1) = real(0.0E+00,kind=8)

560     continue

        u = d(i) - x1 - xu * v
        if ( i /= q ) v = e(i+1)

580     continue

     end do

     if ( u == real(0.0E+00,kind=8) ) then
       u = eps3
     end if

     rv1(q) = u
     rv2(q) = real(0.0E+00,kind=8)
     rv3(q) = real(0.0E+00,kind=8)
!
!  Back substitution.
!
600   continue

  do ii = p, q
    i = p + q - ii
    rv6(i) = ( rv6(i) - u * rv2(i) - v * rv3(i) ) / rv1(i)
    v = u
    u = rv6(i)
  end do
!
!  Orthogonalize with respect to previous members of group.
!
     j = r

     do jj = 1, group

       do

         j = j - 1

         if ( ind(j) == tag ) then
           exit
         end if

       end do

       xu = dot_product ( rv6(p:q), z(p:q,j) )

       rv6(p:q) = rv6(p:q) - xu * z(p:q,j)

     end do

     norm = sum ( abs ( rv6(p:q) ) )

     if ( norm >= real(1.0E+00,kind=8) ) go to 840
!
!  Forward substitution.
!
     if ( its == 5 ) go to 830

     if ( norm == real(0.0E+00,kind=8) ) then
       rv6(s) = eps4
       s = s + 1
       if ( s > q ) s = p
       go to 780
     end if

740  continue

     xu = eps4 / norm
     rv6(p:q) = rv6(p:q) * xu
!
!  Elimination operations on next vector iterate.
!
780  continue
!
!  If RV1(I-1) == E(I), a row interchange was performed earlier in the
!  triangularization process.
!
     do i = ip, q

       u = rv6(i)

       if ( rv1(i-1) == e(i) ) then
         u = rv6(i-1)
         rv6(i-1) = rv6(i)
       end if

       rv6(i) = u - rv4(i) * rv6(i-1)

     end do

     its = its + 1
     go to 600
!
!  Set error: non-converged eigenvector.
!
830  continue

     ierr = -r
     xu = real(0.0E+00,kind=8)
     go to 870
!
!  Normalize so that sum of squares is 1 and expand to full order.
!
840  continue

     u = real(0.0E+00,kind=8)
     do i = p, q
       u = pythag ( u, rv6(i) )
     end do

     xu = real(1.0E+00,kind=8) / u

870  continue

     z(1:n,r) = real(0.0E+00,kind=8)
     z(p:q,r) = rv6(p:q) * xu

     x0 = x1

920  continue

  end do

  if ( q < n ) go to 100

  return
end
subroutine tql1 ( n, d, e, ierr )
!
!*******************************************************************************
!
!! TQL1 computes all eigenvalues of a real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric tridiagonal 
!    matrix by the QL method.
!
!  Modified:
!
!    20 March 2002
!
!  References:
!
!    Bowdler, Martin, Reinsch, Wilkinson,
!    Numerische Mathematik,
!    Volume 11, 1968, pages 293-306.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, is the order of the matrix.
!
!    Input/output, real(kind=8) D(N).
!    On input, the diagonal elements of the matrix.
!    On output, the eigenvalues in ascending order.
!    If an error exit is made, the eigenvalues are correct and
!    ordered for indices 1, 2,... IERR-1, but may not be
!    the smallest eigenvalues.
!
!    Input/output, real(kind=8) E(N).  On input, E(2:N) contains the subdiagonal 
!    elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Output, integer IERR, error flag.
!    0, normal return,
!    J, if the J-th eigenvalue has not been determined after 
!    30 iterations.
!
  implicit none
!
  integer n
!
  real(kind=8) c
  real(kind=8) c2
  real(kind=8) c3
  real(kind=8) d(n)
  real(kind=8) dl1
  real(kind=8) e(n)
  real(kind=8) el1
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer l1
  integer l2
  integer m
  integer mml
  real(kind=8) p
  real(kind=8) pythag
  real(kind=8) r
  real(kind=8) s
  real(kind=8) s2
  real(kind=8) tst1
  real(kind=8) tst2
!
  ierr = 0
  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do

  f = real(0.0E+00,kind=8)
  tst1 = real(0.0E+00,kind=8)
  e(n) = real(0.0E+00,kind=8)

  do l = 1, n

    j = 0
    h = abs ( d(l) ) + abs ( e(l) )
    tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
    do m = l, n

      tst2 = tst1 + abs ( e(m) )

      if ( tst2 == tst1 ) then
        exit
      end if

    end do

    if ( m == l ) go to 210

130 continue

    if ( j >= 30 ) then
      ierr = l
      return
    end if

    j = j + 1
!
!  Form the shift.
!
    l1 = l + 1
    l2 = l1 + 1
    g = d(l)
    p = ( d(l1) - g ) / ( real(2.0E+00,kind=8) * e(l) )
    r = pythag ( p, real(1.0E+00,kind=8) )
    d(l) = e(l) / ( p + sign ( r, p ) )
    d(l1) = e(l) * ( p + sign ( r, p ) )
    dl1 = d(l1)
    h = g - d(l)

    d(l2:n) = d(l2:n) - h

    f = f + h
!
!  QL transformation.
!
    p = d(m)
    c = real(1.0E+00,kind=8)
    c2 = c
    el1 = e(l1)
    s = real(0.0E+00,kind=8)
    mml = m - l

    do ii = 1, mml
      c3 = c2
      c2 = c
      s2 = s
      i = m - ii
      g = c * e(i)
      h = c * p
      r = pythag ( p, e(i) )
      e(i+1) = s * r
      s = e(i) / r
      c = p / r
      p = c * d(i) - s * g
      d(i+1) = h + s * ( c * g + s * d(i) )
    end do

    p = - s * s2 * c3 * el1 * e(l) / dl1
    e(l) = s * p
    d(l) = c * p
    tst2 = tst1 + abs ( e(l) )
    if ( tst2 > tst1 ) go to 130

210 continue

    p = d(l) + f
!
!  Order the eigenvalues.
!
    do ii = 2, l
      i = l + 2 - ii
      if ( p >= d(i-1) ) then
        go to 270
      end if
      d(i) = d(i-1)
    end do

    i = 1

270 continue

    d(i) = p

  end do

  return
end
subroutine tql2 ( nm, n, d, e, z, ierr )
!
!*******************************************************************************
!
!! TQL2 computes all eigenvalues/vectors, real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors of a symmetric 
!    tridiagonal matrix by the QL method.  The eigenvectors of a full 
!    symmetric matrix can also be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    Bowdler, Martin, Reinsch, Wilkinson,
!    TQL2,
!    Numerische Mathematik,
!    Volume 11, pages 293-306, 1968.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) D(N).  On input, the diagonal elements of the matrix.
!    On output, the eigenvalues in ascending order.  If an error exit is 
!    made, the eigenvalues are correct but unordered for indices 1,2,...,IERR-1.
!
!    Input/output, real(kind=8) E(N).  On input, E(2:N) contains the subdiagonal 
!    elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Input, real(kind=8) Z(NM,N).  On input, the transformation matrix produced in 
!    the reduction by TRED2, if performed.  If the eigenvectors of the 
!    tridiagonal matrix are desired, Z must contain the identity matrix.
!    On output, Z contains the orthonormal eigenvectors of the symmetric
!    tridiagonal (or full) matrix.  If an error exit is made, Z contains 
!    the eigenvectors associated with the stored eigenvalues.
!
!    Output, integer IERR, error flag.
!    0, normal return,
!    J, if the J-th eigenvalue has not been determined after 
!    30 iterations.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) c
  real(kind=8) c2
  real(kind=8) c3
  real(kind=8) d(n)
  real(kind=8) dl1
  real(kind=8) e(n)
  real(kind=8) el1
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  integer i
  integer ierr
  integer ii
  integer j
  integer k
  integer l
  integer l1
  integer l2
  integer m
  integer mml
  real(kind=8) p
  real(kind=8) pythag
  real(kind=8) r
  real(kind=8) s
  real(kind=8) s2
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) z(nm,n)
!
  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do

  f = real(0.0E+00,kind=8)
  tst1 = real(0.0E+00,kind=8)
  e(n) = real(0.0E+00,kind=8)

  do l = 1, n

     j = 0
     h = abs ( d(l) ) + abs ( e(l) )
     tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
     do m = l, n
       tst2 = tst1 + abs ( e(m) )
       if ( tst2 == tst1 ) then
         exit
       end if
     end do

     if ( m == l ) go to 220

 130 continue

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift.
!
     l1 = l + 1
     l2 = l1 + 1
     g = d(l)
     p = ( d(l1) - g ) / ( real(2.0E+00,kind=8) * e(l) )
     r = pythag ( p, real(1.0E+00,kind=8) )
     d(l) = e(l) / ( p + sign ( r, p ) )
     d(l1) = e(l) * ( p + sign ( r, p ) )
     dl1 = d(l1)
     h = g - d(l)
     d(l2:n) = d(l2:n) - h
     f = f + h
!
!  QL transformation.
!
     p = d(m)
     c = real(1.0E+00,kind=8)
     c2 = c
     el1 = e(l1)
     s = real(0.0E+00,kind=8)
     mml = m - l

     do ii = 1, mml

        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = pythag ( p, e(i) )
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * ( c * g + s * d(i) )
!
!  Form vector.
!
        do k = 1, n
          h = z(k,i+1)
          z(k,i+1) = s * z(k,i) + c * h
          z(k,i) = c * z(k,i) - s * h
        end do

     end do

     p = - s * s2 * c3 * el1 * e(l) / dl1
     e(l) = s * p
     d(l) = c * p
     tst2 = tst1 + abs ( e(l) )

     if ( tst2 > tst1 ) then 
       go to 130
     end if

220  continue

     d(l) = d(l) + f

  end do
!
!  Order eigenvalues and eigenvectors.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n

      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if

    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      do j = 1, n
        call r_swap ( z(j,i), z(j,k) )
      end do

    end if

  end do

  return
end
subroutine tqlrat ( n, d, e2, ierr )
!
!*******************************************************************************
!
!! TQLRAT compute all eigenvalues of a real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric
!    tridiagonal matrix by the rational QL method.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    C Reinsch,
!    Algorithm 464, TQLRAT,
!    Communications of the ACM,
!    Volume 16, page 689, 1973.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) D(N).  On input, D contains the diagonal elements
!    of the matrix.  On output, D contains the eigenvalues in ascending
!    order.  If an error exit was made, then the eigenvalues are correct
!    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, real(kind=8) E2(N), contains in positions 2 through N the
!    squares of the subdiagonal elements of the matrix.  E2(1) is
!    arbitrary.  On output, E2 has been overwritten by workspace
!    information.
!
!    Output, integer IERR, error flag.
!    0, for no error,
!    J, if the J-th eigenvalue could not be determined after 30 iterations.
!
  implicit none
!
  integer n
!
  real(kind=8) b
  real(kind=8) c
  real(kind=8) d(n)
  real(kind=8) e2(n)
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer l1
  integer m
  integer mml
  real(kind=8) p
  real(kind=8) pythag
  real(kind=8) r
  real(kind=8) s
  real(kind=8) t
!
  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e2(i-1) = e2(i)
  end do

  f = real(0.0E+00,kind=8)
  t = real(0.0E+00,kind=8)
  e2(n) = real(0.0E+00,kind=8)

  do l = 1, n

     j = 0
     h = abs ( d(l) ) + sqrt ( e2(l) )

     if ( t <= h ) then

       t = h
       b = abs ( t ) * epsilon ( b )
       c = b * b

     end if
!
!  Look for small squared sub-diagonal element.
!
     do m = l, n
       if ( e2(m) <= c ) then
         exit
       end if
     end do

     if ( m == l ) go to 210

130  continue

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift.
!
     l1 = l + 1
     s = sqrt ( e2(l) )
     g = d(l)
     p = ( d(l1) - g ) / ( real(2.0E+00,kind=8) * s )
     r = pythag ( p, real(1.0E+00,kind=8) )
     d(l) = s / ( p + sign ( r, p ) )
     h = g - d(l)
     d(l1:n) = d(l1:n) - h
     f = f + h
!
!  Rational QL transformation.
!
     g = d(m)
     if ( g == real(0.0E+00,kind=8) ) g = b
     h = g
     s = real(0.0E+00,kind=8)
     mml = m - l

     do ii = 1, mml
       i = m - ii
       p = g * h
       r = p + e2(i)
       e2(i+1) = s * r
       s = e2(i) / r
       d(i+1) = h + s * ( h + d(i) )
       g = d(i) - e2(i) / g
       if ( g == real(0.0E+00,kind=8) ) g = b
       h = g * p / r
     end do

     e2(l) = s * g
     d(l) = h
!
!  Guard against underflow in convergence test.
!
     if ( h == real(0.0E+00,kind=8) ) go to 210
     if ( abs ( e2(l) ) <= abs ( c / h ) ) go to 210
     e2(l) = h * e2(l)
     if ( e2(l) /= real(0.0E+00,kind=8) ) go to 130

210  continue

     p = d(l) + f
!
!  Order the eigenvalues.
!
     do ii = 2, l
       i = l + 2 - ii
       if ( p >= d(i-1) ) go to 270
       d(i) = d(i-1)
     end do

250  continue

     i = 1

270  continue

     d(i) = p

290  continue

  end do

  return
end
subroutine trbak1 ( nm, n, a, e, m, z ) 
!
!*******************************************************************************
!
!! TRBAK1 determines eigenvectors by undoing the TRED1 transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a real(kind=8) symmetric
!    matrix by back transforming those of the corresponding
!    symmetric tridiagonal matrix determined by TRED1.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) A(NM,N), contains information about the orthogonal 
!    transformations used in the reduction by TRED1 in its strict lower
!    triangle.
!
!    Input, real(kind=8) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is arbitrary.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) Z(NM,M).  On input, the eigenvectors to be back 
!    transformed.  On output, the transformed eigenvectors.
!
  implicit none
!
  integer m
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) e(n)
  integer i
  integer j
  integer l
  real(kind=8) s
  real(kind=8) z(nm,m)
!
  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if

  do i = 2, n

    l = i - 1

    if ( e(i) /= real(0.0E+00,kind=8) ) then

      do j = 1, m

        s = dot_product ( a(i,1:l), z(1:l,j) )

        s = ( s / a(i,l) ) / e(i)

        z(1:l,j) = z(1:l,j) + s * a(i,1:l)

      end do

    end if

  end do

  continue

  return
end
subroutine trbak3 ( nm, n, nv, a, m, z )
!
!*******************************************************************************
!
!! TRBAK3 determines eigenvectors by undoing the TRED3 transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a real(kind=8) symmetric
!    matrix by back transforming those of the corresponding
!    symmetric tridiagonal matrix determined by TRED3.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NV, the dimension of the array paramater A,
!    which must be at least N*(N+1)/2.
!
!    Input, real(kind=8) A(NV), information about the orthogonal transformations
!    used in the reduction by TRED3.
!
!    Input, integer M, the number of eigenvectors to be back transformed.
!
!    Input/output, real(kind=8) Z(NM,M).  On input, the eigenvectors to be back 
!    transformed.  On output, the transformed eigenvectors.
!
  implicit none
!
  integer m
  integer nm
  integer nv
!
  real(kind=8) a(nv)
  real(kind=8) h
  integer i
  integer ik
  integer iz
  integer j
  integer k
  integer l
  integer n
  real(kind=8) s
  real(kind=8) z(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  do i = 2, n

    l = i - 1
    iz = ( i * l ) / 2
    ik = iz + i
    h = a(ik)

    if ( h /= real(0.0E+00,kind=8) ) then

      do j = 1, m

        s = real(0.0E+00,kind=8)
        ik = iz

        do k = 1, l
          ik = ik + 1
          s = s + a(ik) * z(k,j)
        end do

        s = ( s / h ) / h
        ik = iz

        do k = 1, l
          ik = ik + 1
          z(k,j) = z(k,j) - s * a(ik)
        end do

      end do

    end if

  end do

  return
end
subroutine tred1 ( nm, n, a, d, e, e2 )
!
!*******************************************************************************
!
!! TRED1 transforms a real(kind=8) symmetric matrix to symmetric tridiagonal form.
!
!
!  Discussion:
!
!    The routine reduces a real(kind=8) symmetric matrix to a symmetric 
!    tridiagonal matrix using orthogonal similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    Martin, Reinsch, Wilkinson,
!    TRED1,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the array A.
!    NM must be at least N.
!
!    Input, integer N, the order of the matrix A.
!
!    Input/output, real(kind=8) A(NM,N), on input, contains the real(kind=8) symmetric matrix.  
!    Only the lower triangle of the matrix need be supplied.
!    On output, A contains information about the orthogonal transformations 
!    used in the reduction in its strict lower triangle.  
!    The full upper triangle of A is unaltered.
!
!    Output, real(kind=8) D(N), contains the diagonal elements of the tridiagonal 
!    matrix.
!
!    Output, real(kind=8) E(N), contains the subdiagonal elements of the tridiagonal
!    matrix in its last n-1 positions.  e(1) is set to zero.
!
!    Output, real(kind=8) E2(N), contains the squares of the corresponding 
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  integer i
  integer ii
  integer j
  integer k
  integer l
  real(kind=8) scale
!
  d(1:n) = a(n,1:n)

  do i = 1, n
    a(n,i) = a(i,i)
  end do

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = real(0.0E+00,kind=8)
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    if ( scale == real(0.0E+00,kind=8) ) then

      do j = 1, l
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = real(0.0E+00,kind=8)
      end do

      e(i) = real(0.0E+00,kind=8)
      e2(i) = real(0.0E+00,kind=8)

      cycle

    end if

    d(1:l) = d(1:l) / scale

    do k = 1, l
      h = h + d(k)**2
    end do

    e2(i) = h * scale**2
    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g

    if ( l >= 1 ) then
!
!  Form A * U.
!
      e(1:l) = real(0.0E+00,kind=8)

      do j = 1, l

        f = d(j)
        g = e(j) + a(j,j) * f

        do k = j+1, l
          g = g + a(k,j) * d(k)
          e(k) = e(k) + a(k,j) * f
        end do

        e(j) = g

      end do
!
!  Form P.
!
      f = real(0.0E+00,kind=8)

      do j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
      end do

      h = f / ( h + h )
!
!  Form Q.
!
      e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
      do j = 1, l

        f = d(j)
        g = e(j)

        a(j:l,j) = a(j:l,j) - f * e(j:l) - g * d(j:l)

      end do

    end if

    do j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
    end do


  end do

  return
end
subroutine tred2 ( nm, n, a, d, e, z )
!
!*******************************************************************************
!
!! TRED2 transforms a real(kind=8) symmetric matrix to symmetric tridiagonal form.
!
!
!  Discussion:
!
!    This subroutine reduces a real(kind=8) symmetric matrix to a
!    symmetric tridiagonal matrix using and accumulating
!    orthogonal similarity transformations.
!
!    A and Z may coincide, in which case a single storage area is used
!    for the input of A and the output of Z.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    Martin, Reinsch, Wilkinson,
!    TRED2,
!    Numerische Mathematik, 
!    Volume 11, pages 181-195, 1968.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real(kind=8) A(NM,N), the real(kind=8) symmetric input matrix.  Only the
!    lower triangle of the matrix need be supplied.
!
!    Output, real(kind=8) D(N), the diagonal elements of the tridiagonal matrix.
!
!    Output, real(kind=8) E(N), contains the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is set to zero.
!
!    Output, real(kind=8) Z(NM,N), the orthogonal transformation matrix produced 
!    in the reduction.
!
  implicit none
!
  integer n
  integer nm
!
  real(kind=8) a(nm,n)
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  real(kind=8) hh
  integer i
  integer ii
  integer j
  integer k
  integer l
  real(kind=8) scale
  real(kind=8) z(nm,n)
!
  do i = 1, n
    z(i:n,i) = a(i:n,i)
  end do

  d(1:n) = a(n,1:n)

  do ii = 2, n

    i = n + 2 - ii
    l = i - 1
    h = real(0.0E+00,kind=8)
    scale = real(0.0E+00,kind=8)
!
!  Scale row.
!
    do k = 1, l
      scale = scale + abs ( d(k) )
    end do

    if ( scale == real(0.0E+00,kind=8) ) then

      e(i) = d(l)

      do j = 1, l
        d(j) = z(l,j)
        z(i,j) = real(0.0E+00,kind=8)
        z(j,i) = real(0.0E+00,kind=8)
      end do

      go to 290

    end if

    d(1:l) = d(1:l) / scale

    h = h + dot_product ( d(1:l), d(1:l) )

    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g
!
!  Form A*U.
!
    e(1:l) = real(0.0E+00,kind=8)

    do j = 1, l

      f = d(j)
      z(j,i) = f
      g = e(j) + z(j,j) * f

      do k = j+1, l
        g = g + z(k,j) * d(k)
        e(k) = e(k) + z(k,j) * f
      end do

      e(j) = g

    end do
!
!  Form P.
!
    e(1:l) = e(1:l) / h

    f = dot_product ( e(1:l), d(1:l) )

    hh = real(0.5E+00,kind=8) * f / h
!
!  Form Q.
!
    e(1:l) = e(1:l) - hh * d(1:l)
!
!  Form reduced A.
!
    do j = 1, l

      f = d(j)
      g = e(j)

      z(j:l,j) = z(j:l,j) - f * e(j:l) - g * d(j:l)

      d(j) = z(l,j)
      z(i,j) = real(0.0E+00,kind=8)

    end do

290 continue

    d(i) = h


  end do
!
!  Accumulation of transformation matrices.
!
  do i = 2, n

    l = i - 1
    z(n,l) = z(l,l)
    z(l,l) = real(1.0E+00,kind=8)
    h = d(i)
    if ( h /= real(0.0E+00,kind=8) ) then

      d(1:l) = z(1:l,i) / h

      do j = 1, l

        g = dot_product ( z(1:l,i), z(1:l,j) )

        do k = 1, l
          z(k,j) = z(k,j) - g * d(k)
        end do

      end do

    end if

    z(1:l,i) = real(0.0E+00,kind=8)

  end do

  d(1:n) = z(n,1:n)

  z(n,1:n-1) = real(0.0E+00,kind=8)
  z(n,n) = real(1.0E+00,kind=8)

  e(1) = real(0.0E+00,kind=8)

  return
end
subroutine tred3 ( n, nv, a, d, e, e2 )
!
!*******************************************************************************
!
!! TRED3 transforms a real(kind=8) symmetric packed matrix to symmetric tridiagonal form.
!
!
!  Discussion:
!
!    This subroutine reduces a real(kind=8) symmetric matrix, stored as
!    a one-dimensional array, to a symmetric tridiagonal matrix
!    using orthogonal similarity transformations.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    Martin, Reinsch, Wilkinson,
!    TRED3,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NV, the dimension of A, which must be at least
!    (N*(N+1))/2.
!
!    Input/output, real(kind=8) A(NV).  On input, the lower triangle of the real(kind=8) 
!    symmetric matrix, stored row-wise.  On output, information about the 
!    orthogonal transformations used in the reduction.
!
!    Output, real(kind=8) D(N), the diagonal elements of the tridiagonal matrix.
!
!    Output, real(kind=8) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is set to zero.
!
!    Output, real(kind=8) E2(N),  the squares of the corresponding elements of E.
!    E2 may coincide with E if the squares are not needed.
!
  implicit none
!
  integer n
  integer nv
!
  real(kind=8) a(nv)
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) f
  real(kind=8) g
  real(kind=8) h
  real(kind=8) hh
  integer i
  integer ii
  integer iz
  integer j
  integer jk
  integer k
  integer l
  real(kind=8) scale
!
  do ii = 1, n

     i = n + 1 - ii
     l = i - 1
     iz = ( i * l ) / 2
     h = real(0.0E+00,kind=8)
     scale = real(0.0E+00,kind=8)
!
!  Scale row.
!
     do k = 1, l
       iz = iz + 1
       d(k) = a(iz)
       scale = scale + abs ( d(k) )
     end do

     if ( scale == real(0.0E+00,kind=8) ) then
       e(i) = real(0.0E+00,kind=8)
       e2(i) = real(0.0E+00,kind=8)
       go to 290
     end if

     do k = 1, l
       d(k) = d(k) / scale
       h = h + d(k)**2
     end do

     e2(i) = scale * scale * h
     f = d(l)
     g = - sign ( sqrt ( h ), f )
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
     a(iz) = scale * d(l)

     if ( l == 1 ) go to 290

     jk = 1

     do j = 1, l

        f = d(j)
        g = real(0.0E+00,kind=8)

        do k = 1, j-1
          g = g + a(jk) * d(k)
          e(k) = e(k) + a(jk) * f
          jk = jk + 1
        end do

        e(j) = g + a(jk) * f
        jk = jk + 1

     end do
!
!  Form P.
!
     e(1:l) = e(1:l) / h
     f = dot_product ( e(1:l), d(1:l) )
     hh = f / ( h + h )
!
!  Form Q.
!
     e(1:l) = e(1:l) - hh * d(1:l)
     jk = 1
!
!  Form reduced A.
!
     do j = 1, l
       f = d(j)
       g = e(j)
       do k = 1, j
         a(jk) = a(jk) - f * e(k) - g * d(k)
         jk = jk + 1
       end do
     end do

290  continue

     d(i) = a(iz+1)
     a(iz+1) = scale * sqrt ( h )

300  continue

  end do

  return
end
subroutine tridib ( n, eps1, d, e, e2, lb, ub, m11, m, w, ind, ierr )
!
!*******************************************************************************
!
!! TRIDIB computes some eigenvalues of a real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds those eigenvalues of a tridiagonal
!    symmetric matrix between specified boundary indices,
!    using bisection.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) EPS1.  On input, an absolute error tolerance for 
!    the computed eigenvalues.  It should be chosen commensurate with
!    relative perturbations in the matrix elements of the order of the 
!    relative machine precision.  If the input EPS1 is non-positive, it 
!    is reset for each submatrix to a default value, namely, minus the
!    product of the relative machine precision and the 1-norm of the submatrix.
!
!    Input, real(kind=8) D(N), the diagonal elements of the input matrix.
!
!    Input, real(kind=8) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, real(kind=8) E2(N).  On input, the squares of the corresponding 
!    elements of E.  E2(1) is arbitrary.  On output, elements of E2 
!    corresponding to elements of E regarded as negligible, have been 
!    replaced by zero, causing the matrix to split into a direct sum of 
!    submatrices.  E2(1) is also set to zero.
!
!    Input, integer M11, the lower boundary index for the desired eigenvalues.
!
!    Input, integer M, the number of eigenvalues desired.  The upper
!    boundary index M22 is then obtained as M22 = M11 + M - 1.
!
!    Output, real(kind=8) LB, UB, define an interval containing exactly the desired
!    eigenvalues.
!
!    Output, real(kind=8) W(M), the eigenvalues between indices M11 and M22 
!    in ascending order.
!
!    Output, integer IND(M), the submatrix indices associated with the 
!    corresponding eigenvalues in W: 1 for eigenvalues belonging to the 
!    first submatrix from the top, 2 for those belonging to the second 
!    submatrix, and so on.
!
!    Output, integer IERR, error flag.
!    0, for normal return,
!    3*N+1, if multiple eigenvalues at index M11 make unique selection 
!      impossible,
!    3*N+2, if multiple eigenvalues at index M22 make unique selection
!      impossible.
!
  implicit none
!
  integer m
  integer n
!
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) eps1
  integer i
  integer ierr
  integer ii
  integer ind(m)
  integer isturm
  integer j
  integer k
  integer l
  real(kind=8) lb
  integer m1
  integer m11
  integer m2
  integer m22
  integer p
  integer q
  integer r
  real(kind=8) rv4(n)
  real(kind=8) rv5(n)
  integer s
  real(kind=8) t1
  real(kind=8) t2
  integer tag
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) u
  real(kind=8) ub
  real(kind=8) v
  real(kind=8) w(m)
  real(kind=8) x0
  real(kind=8) x1
  real(kind=8) xu
!
  ierr = 0
  tag = 0
  xu = d(1)
  x0 = d(1)
  s = 0
  u = real(0.0E+00,kind=8)
!
!  Look for small sub-diagonal entries and determine an
!  interval containing all the eigenvalues.
!
  do i = 1, n

     x1 = u

     if ( i == n ) then
       u = real(0.0E+00,kind=8)
     else
       u = abs ( e(i+1) )
     end if

     xu = min ( xu, d(i)-(x1+u) )
     x0 = max ( x0, d(i)+(x1+u) )

     if ( i >= 1 ) then
       tst1 = abs ( d(i) ) + abs ( d(i-1) )
       tst2 = tst1 + abs ( e(i) )
       if ( tst2 <= tst1 ) then
         e2(i) = real(0.0E+00,kind=8)
       end if
     else
       e2(i) = real(0.0E+00,kind=8)
     end if

  end do

  x1 = n
  x1 = x1 * max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )
  xu = xu - x1
  t1 = xu
  x0 = x0 + x1
  t2 = x0
!
!  Determine an interval containing exactly the desired eigenvalues.
!
  p = 1
  q = n
  m1 = m11 - 1
  if ( m1 == 0 ) go to 75
  isturm = 1

50 continue

  v = x1
  x1 = xu + (x0 - xu) * real(0.5E+00,kind=8)
  if ( x1 == v ) go to 980
  go to 320

60 continue

  if ( s - m1 ) 65, 73, 70

65 continue

  xu = x1
  go to 50

70 continue

  x0 = x1
  go to 50

73 continue

  xu = x1
  t1 = x1

75 continue

  m22 = m1 + m
  if ( m22 == n ) go to 90
  x0 = t2
  isturm = 2
  go to 50

80 continue

   if ( s - m22 ) 65, 85, 70

85 continue

   t2 = x1

90 continue

  q = 0
  r = 0
!
!  Establish and process next submatrix, refining interval by the 
!  Gerschgorin bounds.
!
100 continue

  if ( r == m ) then
    go to 1001
  end if

  tag = tag + 1
  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = real(0.0E+00,kind=8)

  do q = p, n

    x1 = u
    u = real(0.0E+00,kind=8)
    v = real(0.0E+00,kind=8)

    if ( q < n ) then
      u = abs ( e(q+1) )
      v = e2(q+1)
    end if

    xu = min ( d(q)-(x1+u), xu )
    x0 = max ( d(q)+(x1+u), x0 )

    if ( v == real(0.0E+00,kind=8) ) then
      exit
    end if

  end do

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )
  if ( eps1 <= real(0.0E+00,kind=8) ) eps1 = -x1
  if ( p /= q ) go to 180
!
!  Check for isolated root within interval.
!
  if ( t1 > d(p) .or. d(p) >= t2 ) go to 940
  m1 = p
  m2 = p
  rv5(p) = d(p)
  go to 900

180 continue

  x1 = x1 * (q - p + 1)
  lb = max ( t1, xu-x1 )
  ub = min ( t2, x0+x1 )
  x1 = lb
  isturm = 3
  go to 320

200 continue

  m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320

220 continue

  m2 = s
  if ( m1 > m2 ) go to 940
!
!  Find roots by bisection.
!
  x0 = ub
  isturm = 5

  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for the K-th eigenvalue.
!
  k = m2

250 continue

  xu = lb

  do ii = m1, k

    i = m1 + k - ii
    if ( xu < rv4(i) ) then
      xu = rv4(i)
      exit
    end if

  end do

  if ( x0 > rv5(k) ) x0 = rv5(k)
!
!  Next bisection step.
!
300  continue

     x1 = ( xu + x0 ) * real(0.5E+00,kind=8)
     if ( ( x0 - xu ) <= abs ( eps1) ) go to 420
     tst1 = real(2.0E+00,kind=8) * ( abs ( xu ) + abs ( x0 ) )
     tst2 = tst1 + (x0 - xu)
     if ( tst2 == tst1 ) go to 420
!
!  Sturm sequence.
!
320  continue

     s = p - 1
     u = real(1.0E+00,kind=8)

     do i = p, q

       if ( u == real(0.0E+00,kind=8) ) then
         v = abs ( e(i) ) / epsilon ( v )
         if ( e2(i) == real(0.0E+00,kind=8) ) v = real(0.0E+00,kind=8)
       else
         v = e2(i) / u
       end if

       u = d(i) - x1 - v

       if ( u < real(0.0E+00,kind=8) ) then
         s = s + 1
       end if

     end do

     go to (60,80,200,220,360), isturm
!
!  Refine intervals.
!
360  continue

     if ( s >= k) go to 400
     xu = x1
     if ( s >= m1) go to 380
     rv4(m1) = x1
     go to 300

380  continue

     rv4(s+1) = x1
     if ( rv5(s) > x1) rv5(s) = x1
     go to 300

400  continue

     x0 = x1
     go to 300
!
!  K-th eigenvalue found.
!
420  continue

  rv5(k) = x1
  k = k - 1
  if ( k >= m1 ) go to 250
!
!  Order eigenvalues tagged with their submatrix associations.
!
900 continue

  s = r
  r = r + m2 - m1 + 1
  j = 1
  k = m1

  do l = 1, r

     if ( j > s ) go to 910
     if ( k > m2 ) go to 940
     if ( rv5(k) >= w(l) ) go to 915

     do ii = j, s
       i = l + s - ii
       w(i+1) = w(i)
       ind(i+1) = ind(i)
     end do

910  continue

     w(l) = rv5(k)
     ind(l) = tag
     k = k + 1
     go to 920

915  continue

     j = j + 1

920  continue

  end do

940 continue

  if ( q < n ) go to 100

  go to 1001
!
!  Set error: interval cannot be found containing exactly the 
!  desired eigenvalues.
!
980 continue

  ierr = 3 * n + isturm

1001 continue

  lb = t1
  ub = t2
  return
end
subroutine tsturm ( nm, n, eps1, d, e, e2, lb, ub, mm, m, w, z, ierr )
!
!*******************************************************************************
!
!! TSTURM computes some eigenvalues/vectors, real(kind=8) symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds those eigenvalues of a tridiagonal
!    symmetric matrix which lie in a specified interval and their
!    associated eigenvectors, using bisection and inverse iteration.
!
!  Modified:
!
!    20 March 2002
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the two-dimensional arrays.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real(kind=8) EPS1.  On input, an absolute error tolerance for 
!    the computed eigenvalues.  It should be chosen commensurate with
!    relative perturbations in the matrix elements of the order of the 
!    relative machine precision.  If the input EPS1 is non-positive, it 
!    is reset for each submatrix to a default value, namely, minus the
!    product of the relative machine precision and the 1-norm of the submatrix.
!
!    Input, real(kind=8) D(N), the diagonal elements of the input matrix.
!
!    Input, real(kind=8) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, real(kind=8) E2(N).  On input, the squares of the corresponding 
!    elements of E.  E2(1) is arbitrary.  On output, elements of E2 
!    corresponding to elements of E regarded as negligible have been 
!    replaced by zero, causing the matrix to split into a direct sum of 
!    submatrices.  E2(1) is also set to zero.
!
!    Input, real(kind=8) LB, UB, define the interval to be searched for eigenvalues.
!    If LB is not less than UB, no eigenvalues will be found.
!
!    Input, integer MM, an upper bound for the number of eigenvalues in 
!    the interval.  If more than MM eigenvalues are determined to lie in 
!    the interval, an error return is made with no values or vectors found.
!
!    Output, integer M, the number of eigenvalues determined to lie 
!    in (LB, UB).
!
!    Output, real(kind=8) W(M), the eigenvalues in ascending order if the matrix
!    does not split.  If the matrix splits, the eigenvalues are in ascending
!    order for each submatrix.  If a vector error exit is made, W contains 
!    those values already found.
!
!    Output, real(kind=8) Z(NM,MM), the associated set of orthonormal eigenvectors.
!    If an error exit is made, Z contains those vectors already found.
!
!    Output, integer IERR, error flag.
!    0, normal return.
!    3*N+1, if M exceeds MM.
!    4*N+R, if the eigenvector corresponding to the R-th
!      eigenvalue fails to converge in 5 iterations.
!
  implicit none
!
  integer mm
  integer n
  integer nm
!
  real(kind=8) d(n)
  real(kind=8) e(n)
  real(kind=8) e2(n)
  real(kind=8) eps1
  real(kind=8) eps2
  real(kind=8) eps3
  real(kind=8) eps4
  integer group
  integer i
  integer ierr
  integer ii
  integer ip
  integer isturm
  integer its
  integer j
  integer jj
  integer k
  real(kind=8) lb
  integer m
  integer m1
  integer m2
  real(kind=8) norm
  integer p
  real(kind=8) pythag
  integer q
  integer r
  real(kind=8) rv1(n)
  real(kind=8) rv2(n)
  real(kind=8) rv3(n)
  real(kind=8) rv4(n)
  real(kind=8) rv5(n)
  real(kind=8) rv6(n)
  integer s
  real(kind=8) t1
  real(kind=8) t2
  real(kind=8) tst1
  real(kind=8) tst2
  real(kind=8) u
  real(kind=8) ub
  real(kind=8) uk
  real(kind=8) v
  real(kind=8) w(mm)
  real(kind=8) x0
  real(kind=8) x1
  real(kind=8) xu
  real(kind=8) z(nm,mm)
!
  ierr = 0
  s = 0
  t1 = lb
  t2 = ub
!
!  Look for small sub-diagonal entries.
!
  e2(1) = real(0.0E+00,kind=8)

  do i = 2, n

    tst1 = abs ( d(i) ) + abs ( d(i-1) )
    tst2 = tst1 + abs ( e(i) )

    if ( tst2 <= tst1 ) then
      e2(i) = real(0.0E+00,kind=8)
    end if

  end do
!
!  Determine the number of eigenvalues in the interval.
!
  p = 1
  q = n
  x1 = ub
  isturm = 1
  go to 320

60 continue

  m = s
  x1 = lb
  isturm = 2
  go to 320

80 continue

  m = m - s

  if ( m > mm ) go to 980

  q = 0
  r = 0
!
!  Establish and process next submatrix, refining interval by the 
!  Gerschgorin bounds.
!
100 continue

  if ( r == m ) go to 1001

  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = real(0.0E+00,kind=8)

  do q = p, n

     x1 = u
     u = real(0.0E+00,kind=8)
     v = real(0.0E+00,kind=8)

     if ( q /= n ) then
       u = abs ( e(q+1) )
       v = e2(q+1)
     end if

     xu = min ( d(q)-(x1+u), xu )
     x0 = max ( d(q)+(x1+u), x0 )

     if ( v == real(0.0E+00,kind=8) ) then
       exit
     end if

  end do

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )

  if ( eps1 <= real(0.0E+00,kind=8) ) then
    eps1 = -x1
  end if

  if ( p /= q ) go to 180
!
!  Check for isolated root within interval.
!
  if ( t1 > d(p) .or. d(p) >= t2 ) go to 940

  r = r + 1

  z(1:n,r) = real(0.0E+00,kind=8)

  w(r) = d(p)
  z(p,r) = real(1.0E+00,kind=8)
  go to 940

180 continue

  u = q - p + 1
  x1 = u * x1
  lb = max ( t1, xu-x1 )
  ub = min ( t2, x0+x1 )
  x1 = lb
  isturm = 3
  go to 320

200 continue

  m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320

220 continue

  m2 = s
  if ( m1 > m2 ) go to 940
!
!  Find roots by bisection.
!
  x0 = ub
  isturm = 5

  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for K-th eigenvalue.
!
  k = m2

250 continue

  xu = lb

  do ii = m1, k

    i = m1 + k - ii

    if ( xu < rv4(i) ) then
      xu = rv4(i)
      exit
    end if

  end do

280 continue

  if ( x0 > rv5(k) ) x0 = rv5(k)
!
!  Next bisection step.
!
300 continue

     x1 = ( xu + x0 ) * real(0.5E+00,kind=8)
     if ( ( x0 - xu ) <= abs ( eps1 ) ) go to 420
     tst1 = real(2.0E+00,kind=8) * ( abs ( xu ) + abs ( x0 ) )
     tst2 = tst1 + (x0 - xu)
     if ( tst2 == tst1 ) go to 420
!
!  Sturm sequence.
!
320  continue

     s = p - 1
     u = real(1.0E+00,kind=8)

     do i = p, q

        if ( u /= real(0.0E+00,kind=8) ) go to 325
        v = abs ( e(i) ) / epsilon ( v )
        if ( e2(i) == real(0.0E+00,kind=8) ) v = real(0.0E+00,kind=8)
        go to 330

325     continue

        v = e2(i) / u
330     continue

        u = d(i) - x1 - v
        if ( u < real(0.0E+00,kind=8) ) s = s + 1

     end do

     go to ( 60,80,200,220,360 ), isturm
!
!  Refine intervals.
!
360  continue

     if ( s >= k ) go to 400
     xu = x1
     if ( s >= m1 ) go to 380
     rv4(m1) = x1
     go to 300

380  continue

     rv4(s+1) = x1
     if ( rv5(s) > x1 ) then
       rv5(s) = x1
     end if
     go to 300

400  continue

     x0 = x1
     go to 300
!
!  K-th eigenvalue found.
!
420  continue

  rv5(k) = x1
  k = k - 1
  if ( k >= m1 ) go to 250
!
!  Find vectors by inverse iteration.
!
  norm = abs ( d(p) )
  ip = p + 1

  do i = ip, q
    norm = max ( norm, abs ( d(i) ) + abs ( e(i) ) )
  end do
!
!  EPS2 is the criterion for grouping,
!  EPS3 replaces zero pivots and equal roots are modified by eps3,
!  EPS4 is taken very small to avoid overflow.
!
  eps2 = 1.0E-03 * norm
  eps3 = abs ( norm ) * epsilon ( eps3 )
  uk = q - p + 1
  eps4 = uk * eps3
  uk = eps4 / sqrt ( uk )
  group = 0
  s = p

  do k = m1, m2

     r = r + 1
     its = 1
     w(r) = rv5(k)
     x1 = rv5(k)
!
!  Look for close or coincident roots.
!
     if ( k /= m1 ) then
       if ( x1 - x0 >= eps2 ) group = -1
       group = group + 1
       if ( x1 <= x0 ) then
         x1 = x0 + eps3
       end if
     end if
!
!  Elimination with interchanges and initialization of vector.
!
520  continue

     v = real(0.0E+00,kind=8)

     do i = p, q

        rv6(i) = uk

        if ( i == p ) go to 560

        if ( abs ( e(i) ) >= abs ( u ) ) then
          xu = u / e(i)
          rv4(i) = xu
          rv1(i-1) = e(i)
          rv2(i-1) = d(i) - x1
          rv3(i-1) = real(0.0E+00,kind=8)
          if ( i /= q ) rv3(i-1) = e(i+1)
          u = v - xu * rv2(i-1)
          v = -xu * rv3(i-1)
          cycle
        end if

540     continue

        xu = e(i) / u
        rv4(i) = xu
        rv1(i-1) = u
        rv2(i-1) = v
        rv3(i-1) = real(0.0E+00,kind=8)

560     continue

        u = d(i) - x1 - xu * v

        if ( i /= q ) then
          v = e(i+1)
        end if

     end do

     if ( u == real(0.0E+00,kind=8) ) u = eps3
     rv1(q) = u
     rv2(q) = real(0.0E+00,kind=8)
     rv3(q) = real(0.0E+00,kind=8)
!
!  Back substitution.
!
600  continue

     do ii = p, q
        i = p + q - ii
        rv6(i) = ( rv6(i) - u * rv2(i) - v * rv3(i) ) / rv1(i)
        v = u
        u = rv6(i)
     end do
!
!  Orthogonalize with respect to previous members of group.
!
     do jj = 1, group
        j = r - group - 1 + jj
        xu = dot_product ( rv6(p:q), z(p:q,j) )
        rv6(p:q) = rv6(p:q) - xu * z(p:q,j)
     end do

700  continue

     norm = sum ( abs ( rv6(p:q) ) )

     if ( norm >= real(1.0E+00,kind=8) ) then
       go to 840
     end if
!
!  Forward substitution.
!
     if ( its == 5 ) then
       ierr = 4 * n + r
       go to 1001
     end if

     if ( norm == real(0.0E+00,kind=8) ) then
       rv6(s) = eps4
       s = s + 1
       if ( s > q ) then
         s = p
       end if
       go to 780
     end if

740  continue

    xu = eps4 / norm

     rv6(p:q) = rv6(p:q) * xu
!
!  Elimination operations on next vector iterate.
!
780    continue
!
!  If rv1(i-1) == e(i), a row interchange was performed earlier in the
!  triangularization process.
!
     do i = p, q

       u = rv6(i)

       if ( rv1(i-1) == e(i) ) then
         u = rv6(i-1)
         rv6(i-1) = rv6(i)
       end if

       rv6(i) = u - rv4(i) * rv6(i-1)

     end do

     its = its + 1
     go to 600
!
!  Normalize so that sum of squares is 1 and expand to full order.
!
840  continue

     u = real(0.0E+00,kind=8)

     do i = p, q
       u = pythag ( u, rv6(i) )
     end do

     xu = real(1.0E+00,kind=8) / u

     z(1:n,r) = real(0.0E+00,kind=8)
     z(p:q,r) = rv6(p:q) * xu

     x0 = x1
  
  end do

940 continue

  if ( q < n ) then
    go to 100
  end if

  go to 1001
!
!  Set error: underestimate of number of eigenvalues in interval.
!
980 continue

  ierr = 3 * n + 1

1001 continue

  lb = t1
  ub = t2

  return
end
