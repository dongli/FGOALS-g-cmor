subroutine calc_gamil_lat_bnds(num_lat, B, lat, lat_bnds)

  use flogger
  use const_mod

  implicit none

  integer, intent(in) :: num_lat
  real(8), intent(in) :: B ! A CONTROL PARAMETER RELATED TO THE MERIDIONAL RESOLUTION
  real(8), intent(in) :: lat(num_lat)
  real(8), intent(inout) :: lat_bnds(num_lat+1)

  real(8) DY            ! MERIDIONAL STEPSIZE
  real(8) YTHU(num_lat) ! LATITUDE AT THE NORNAL MERIDIONAL GRID Yj
  real(8) YTHV(num_lat) ! LATITUDE AT THE HALF-MOVE MERIDIONAL GRID Yj+1/2
  real(8) WTGU(num_lat) ! AREA-WEIGHTING AT THE NORNAL MERIDIONAL GRID Yj
  real(8) WTGV(num_lat) ! AREA-WEIGHTING AT THE HALF-MOVE MERIDIONAL GRID Yj+1/2
!
!     4) WORKING VARIABLES: U,S,DA,A2,U2,A,B,J,M1
!
  real(8) AS ! AREA SIZE WITH RESPECT TO LATITUDE
!
! The resolution formula:
!
!   DY = (180-10B)/(num_lat-1) (in degree) PI(1-B/18)/(num_lat-1) (in arc)
!
!   when B=2.0, num_lat=41, DY=4 degree; when B=2.0, num_lat=81, DY=2 degree
!
  real(8) A  ! A  = B/(0.5*PI) A CONTROL PARAMETER RELATED TO THE AREA-SIZE COMPUTING
  real(8) A2 ! A2 = A*2
  real(8) DA ! DA = 1/A
  real(8) S  ! S  = AS(PI/2)   THE TOTAL AREA-SIZE OF THE WHOLE REGION
  real(8) S1 ! S1 = AS(-PI/3)	 THE AREA-SIZE AT THE POINT THETA=-PI/3
  real(8) S2 ! S2 = AS(PI/3)	 THE AREA-SIZE AT THE POINT THETA=PI/3
  real(8) U1 ! U1 = 1+2B/3     NEEDED WHEN CALCULATE YTHU, YTHV
  real(8) U2 ! U2 = (1-B/3)^2	 NEEDED WHEN CALCULATE YTHU, YTHV
  integer I, J
  integer M1 ! M1 = num_lat-1
!
! Description of the even-area method developed by Bin Wang in 2000:
!
!   Under the weighting w(theta)=1.0-a(|theta|-PI/3) when |theta|>=PI/3
!                       w(theta)=1.0                 when |theta|< PI/3
!
!   The area size AS is:
!     AS(theta)=[(1+a*theta+a*PI/3)^2-(1-b/3)^2]/(2a)  when -PI/2<=theta<=-PI/3
!     AS(theta)=AS(-PI/3)+(theta+PI/3)                 when -PI/3< theta<= PI/3
!     AS(theta)=AS(PI/3)+[1-(1-a*theta+a*PI/3)^2]/(2a) when  PI/3< theta<= PI/2
!     where a= b/(0.5*Pi), 0<b<1, b=1.0-0.25*sqrt(25/(m1-20))
!
!   Suppose the total area-size of the whole region S is partitioned into num_lat-1
!   equal small area: AS(theta(j+1))-AS(theta(j))=DY=constant, then theta(j)
!   can be calculated according to the formula of AS. Obviously,
!   theta(j+1)-theta(j) will not be a constant when they are not in the interval
!   [-PI/3, PI/3]. Especially, closer to poles theta(j) is, bigger
!   theta(j+1)-theta(j) will become. In this way, the physical stepsizes in the
!   polar regions increase and the computational stability becomes better.
!   Note that: the physical mesh is not even, but the computing mesh is, which
!   makes the meridional discretization easy.
!
  integer k

  M1 = num_lat - 1
  A  = B * 2.0D0 / PI
  A2 = A * 2.0D0
  DA = 1.0D0 / A
  S  = PI * (1.0D0 - B / 18.0D0)
  S1 = PI * (1.0D0 - B / 6.0D0) / 6.0D0
  S2 = PI * (5.0D0 - B / 6.0D0) / 6.0D0
  U1 = 1.0D0 + 2.0D0 * B / 3.0D0
  U2 = (1.0D0 - B / 3.0D0) * (1.0D0 - B / 3.0D0)
  DY = S / DFLOAT(M1)
  DO J = 0, M1
    AS = DY * DFLOAT(J)
    IF (AS <= S1) THEN
      YTHU(J+1) = (DSQRT(AS * A2 + U2) - U1) * DA + PI * 0.5D0
      IF (YTHU(J+1) < 0.0) YTHU(J+1) = 0.0D0
      WTGU(J+1) = 1.0D0 - A * (DABS(YTHU(J+1) - PI * 0.5D0) - PI / 3.0D0)
    ELSE IF (AS <= S2) THEN
      YTHU(J+1) = AS - S1 - PI / 3.0D0 + PI * 0.5D0
      WTGU(J+1) = 1.0D0
    ELSE
      YTHU(J+1) = (U1 - DSQRT(1.0D0 - (AS - S2) * A2)) * DA + PI * 0.5D0
      WTGU(J+1) = 1.0D0 - A * (DABS(YTHU(J+1) - PI * 0.5D0) - PI / 3.0D0)
    END IF
    AS = DY * (DFLOAT(J) + 0.5D0)
    IF (AS <= S1) THEN
      YTHV(J+1) = (DSQRT(AS * A2 + U2) - U1) * DA + PI * 0.5D0
      WTGV(J+1) = 1.0D0 - A * (DABS(YTHV(J+1) - PI * 0.5D0) - PI / 3.0D0)
    ELSE IF (AS <= S2) THEN
      YTHV(J+1) = AS - S1 - PI / 3.0 + PI *  0.5
      WTGV(J+1) = 1.0D0
    ELSE IF (J < M1) THEN
      YTHV(J+1) = (U1 - DSQRT(1.0D0 - (AS - S2) * A2)) * DA + PI * 0.5D0
      WTGV(J+1) = 1.0D0 - A * (DABS(YTHV(J+1) - PI * 0.5D0) - PI / 3.0D0)
    END IF
  END DO

  YTHV(num_lat) = PI
  WTGV(num_lat) = 1.0D0 - A * (PI * 0.5 - PI / 3.0)

  do j = 1, num_lat
    if (abs(lat(j) - YTHU(J) * deg + 90.0d0) > 1.0d-12) then
      print *, j, lat(j), YTHU(J) * deg - 90.0d0
      call log_error('Lat grids are not matched!', __FILE__, __LINE__)
    end if
    lat_bnds(j+1) = YTHV(J) * deg - 90.0d0
  end do
  lat_bnds(1) = -90.0d0
  lat_bnds(num_lat+1) = 90.0d0

end subroutine calc_gamil_lat_bnds
