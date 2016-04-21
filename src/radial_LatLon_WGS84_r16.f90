subroutine radial_LatLon_WGS84_r16(Angle, Range, LAT_Radar, LON_Radar, LAT_Ship, LON_Ship)
!	This subroutine !al!ulates Latitude and Longitude
!	from radial !oordinates using WGS84 !oordinates.

!	Input:
real(kind=8)::Angle ! Angle in geographi!al degrees (East = 90 degrees) from radar lo!ation
real(kind=8)::Range ! Distan!e in km from radar lo!ation
real(kind=8)::LAT_Radar ! Latitude of radar lo!ation in degrees
real(kind=8)::LON_Radar ! Longitude of radar lo!ation in degrees

!	Output:
real(kind=8)::LAT_Ship ! Latitude of ship lo!ation in degrees
real(kind=8)::LON_Ship ! Longitude of ship lo!ation in degrees

!	Local variables:
real(kind=8)::PI
real(kind=8)::rad
real(kind=8)::a, f, esq
real(kind=8)::r1, r2, arc
real(kind=8)::blimit, dd_max
real(kind=8)::glat1, glat2, glon1, glon2, faz, baz, edist

logical::err

!----------------------------------------------------------------------------------
!
PI  = acos(0.0)*2.
rad = 180.0 / PI

a = 6378137.0
f = 1.0 / 298.25722210088
esq = f * (2.0 - f)
r1  = 0.0
r2  = PI / 2.0

call gpnarc ( a, f, esq, PI, r1, r2, arc )

blimit = 2.0 * arc - 1.0
dd_max = blimit                  
!
glat1 = LAT_Radar / rad ! LAT_in_radian
glon1 = LON_Radar / rad ! LON_in radian
faz   = Angle / rad ! Azimuth_in_radian
edist = Range * 1000.0 ! Range_in_meters

if(edist.ge.dd_max )then
     err = .true.
     write(*,'(a,f12.3,a)') ' Invalid Distance: ', edist, ' meters'
endif

call dirct1 (glat1,glon1,glat2,glon2,faz,baz,edist)

LAT_Ship = glat2 * rad
LON_Ship = glon2 * rad

return
end subroutine radial_LatLon_WGS84_r16


!**************************************************************************
!
! NAME:        GPNARC
! VERSION:     200005.26
! WRITTEN BY:  ROBERT (Sid) SAFFORD
! PURPOSE:     SUBROUTINE TO COMPUTE THE LENGTH OF A MERIDIONAL AR! 
!              BETWEEN TWO LATITUDES
!
! INPUT PARAMETERS:
! -----------------
! AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
! FLAT         FLATTENING (0.0033528 ... )
! ESQ          E!!ENTRI!ITY SQUARED FOR REFERENCE ELLIPSOID
! PI           3.14159...
! P1           LAT STATION 1
! P2           LAT STATION 2
!
! OUTPUT PARAMETERS:
! ------------------
! ARC          GEODETIC DISTAN!E 
!
!********1*********2*********3*********4*********5*********6*********7*
! 
subroutine GPNARC(AMAX, FLAT, ESQ, PI, P1, P2, ARC)

real(kind=8)::AMAX
real(kind=8)::FLAT
real(kind=8)::ESQ
real(kind=8)::PI
real(kind=8)::P1, P2
real(kind=8)::ARC

real(kind=8)::TT
real(kind=8)::S1, S2, DA
real(kind=8)::E2, E4, E6, E8, EX
real(kind=8)::T1, T2, T3, T4, T5 
real(kind=8)::A
real(kind=8)::B, C, D, E, F
real(kind=8)::DB, DC, DD, DE, DF
logical::FLAG

!-----------------------------------------------------------
!
!       CHECK FOR A 90 DEGREE LOOKUP
!
TT   = 5.0d-15
FLAG = .FALSE.
!
S1 = ABS(P1)
S2 = ABS(P2)
!
IF( (PI/2.0-TT).LT.S2 .AND. S2.LT.(PI/2.0+TT) )THEN
FLAG = .TRUE.
END IF
!
IF( S1.GT.TT )THEN
FLAG = .FALSE.
END IF
!
DA = P2 - P1
S1 = 0.0
S2 = 0.0
!
!       !OMPUTE THE LENGTH OF A MERIDIONAL AR! BETWEEN TWO LATITUDES
!
E2 = ESQ
E4 = E2 * E2
E6 = E4 * E2
E8 = E6 * E2
EX = E8 * E2

T1 = E2 * (003.0 /      4.0)
T2 = E4 * (015.0 /     64.0)
T3 = E6 * (035.0 /    512.0)
T4 = E8 * (315.0 /  16384.0)
T5 = EX * (693.0 / 131072.0)

A  = 1.0 + T1 + 3.0 * T2 + 10.0 * T3 + 35.0 * T4 + 126.0 * T5

if(.not. FLAG) then
     B  = T1 +  4.0 * T2 + 15.0 * T3 +  56.0 * T4 + 210.0 * T5
     C  = T2 +  6.0 * T3 + 28.0 * T4 + 120.0 * T5
     D  = T3 +  8.0 * T4 + 45.0 * T5
     E  = T4 + 10.0 * T5
     F  = T5

     DB = SIN(P2 *  2.0) - SIN(P1 *  2.0)
     DC = SIN(P2 *  4.0) - SIN(P1 *  4.0)
     DD = SIN(P2 *  6.0) - SIN(P1 *  6.0)
     DE = SIN(P2 *  8.0) - SIN(P1 *  8.0)
     DF = SIN(P2 * 10.0) - SIN(P1 * 10.0)

     !	 COMPUTE THE S2 PART OF THE SERIES EXPANSION
     !
     S2 = -DB * B / 2.0 + DC * C / 4.0 - DD*D / 6.0 + DE * E / 8.0 - DF * F / 10.0
end if

!	COMPUTE THE S1 PART OF THE SERIES EXPANSION
!
S1 = DA * A
!
!	COMPUTE THE ARC LENGTH
!
ARC = AMAX * (1.0 - ESQ) * (S1 + S2)
!
RETURN

END subroutine GPNARC


!**************************************************************************
!
SUBROUTINE DIRCT1(glat1,glon1,glat2,glon2,faz,baz,s)
!
! *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY
! *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
! *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
!
! *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
! *** F IS THE FLATTENING OF THE REFERENCE ELLIPSOID
! *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
! *** AZIMUTHS IN RADIANS CLO!KWISE FROM NORTH
! *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A
!
! *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75
! *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608
!
!	Input:
real(kind=8)::glat1 ! LAT location1 in radian
real(kind=8)::glon1 ! LON location1 in radian
real(kind=8)::faz   ! Azimuth forward in radian
real(kind=8)::s     ! Distance in meters

!	Output:
real(kind=8)::glat2 ! LAT location2 in radian
real(kind=8)::glon2 ! LON location2 in radian
real(kind=8)::baz   ! Azimuth back in radian

real(kind=8)::PI, EPS
real(kind=8)::A, F

real(kind=8)::r, tu 
real(kind=8)::sf, cf, cu, su, sa, c2a
real(kind=8)::x, c, d, y
real(kind=8)::sy, cy, cz
real(kind=8)::e

PI  = 3.14159265358979323846264338
EPS = 0.5d-13
A = 6378137.0
F = 1.0 / 298.25722210088

r   = 1.0 - F
tu  = r * sin(glat1) / cos(glat1)
sf  = sin(faz)
cf  = cos(faz)
baz = 0.0
if(cf.ne.0.0) baz = 2.0 * atan2(tu,cf)
cu  = 1.0 / sqrt(tu * tu + 1.0)
su  = tu * cu
sa  = cu * sf
c2a = -1.0 * sa*sa + 1.0
x   = sqrt((1.0 / R/R - 1.0) * c2a + 1.0) + 1.0
x   = (x - 2.0) / x
c   = 1.0 - x
c   = (x*x / 4.0 + 1.0) / c
d   = (0.375 * x*x - 1.0) * x
tu  = s / r / A / c
y   = tu

do while (abs(y-c) .gt. EPS)
     sy  = sin(y)
     cy  = cos(y)
     cz  = cos(baz + y)
     e   = cz*cz * 2.0 - 1.0
     c   = y
     x   = e * cy
     y   = e + e - 1.0
     y   = (((sy*sy * 4.0 - 3.0) * y * cz * d / 6.0 + x) * d / 4.0 - cz) * sy * d + tu
enddo

baz   = cu * cy * cf - su * sy
c     = r * sqrt(sa*sa + baz*baz)
d     = su * cy + cu * sy * cf
glat2 = atan2(d,c)
c     = cu * cy - su * sy * cf
x     = atan2(sy*sf,c)
c     = ((-3.0 * c2a + 4.0) * F + 4.0) * c2a * F / 16.0
d     = ((e * cy * c + cz) * sy * c + y) * sa
glon2 = glon1 + x -(1.0 - c) * d * F
baz   = atan2(sa,baz) + PI

return
end SUBROUTINE DIRCT1
