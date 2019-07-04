! ------------------------------------------------------------------------------
! Description:
!
!   This module provides several useful interpolation routines.
!
! Authors:
!
!   Li Dong - dongli@lasg.iap.ac.cn
! ------------------------------------------------------------------------------

module interp_mod

  use const_mod

  implicit none

contains

  subroutine interp_bilinear(x1, y1, data1, x2, y2, data2, xspan, ierr)

    real(8), intent(in) :: x1(:), y1(:), data1(:,:), x2(:), y2(:)
    real(8), intent(out) :: data2(:,:)
    real(8), intent(in) :: xspan
    integer, intent(out), optional :: ierr

    integer nx1, ny1, nx2, ny2
    integer size(1), i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(8) tmp1, tmp2
    logical is_found

    size = shape(x1); nx1 = size(1)
    size = shape(y1); ny1 = size(1)
    size = shape(x2); nx2 = size(1)
    size = shape(y2); ny2 = size(1)

    allocate(i1(nx2), i2(nx2))
    allocate(j1(ny2), j2(ny2))
    allocate(xwgt1(nx2), xwgt2(nx2))
    allocate(ywgt1(ny2), ywgt2(ny2))

    do i = 1, nx2
      is_found = .false.
      do ii = 1, nx1-1
        if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
          i1(i) = ii
          i2(i) = ii+1
          tmp1 = x1(i1(i))
          tmp2 = x1(i2(i))
          is_found = .true.
          exit
        else if (x2(i) >= x1(nx1)) then
          i1(i) = nx1
          i2(i) = 1
          tmp1 = x1(i1(i))
          tmp2 = x1(i2(i))+xspan
          is_found = .true.
          exit
        end if
      end do
      if (.not. is_found) then
        if (present(ierr)) then
          ierr = -1
          return
        else
          write(6, "('[Error]: interp_bilinear: Longitude mismatch!')")
          stop
        end if
      end if
      xwgt1(i) = (tmp2-x2(i))/(tmp2-tmp1)
      xwgt2(i) = (x2(i)-tmp1)/(tmp2-tmp1)
    end do

    do j = 1, ny2
      is_found = .false.
      do jj = 1, ny1-1
        if ((y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) .or. &
            (jj == ny1-1 .and. (y2(j) <= y1(ny1) .or. abs(y2(j)-y1(ny1)) < 1.0e-6))) then
          j1(j) = jj
          j2(j) = jj+1
          is_found = .true.
          exit
        end if
      end do
      if (.not. is_found) then
        if (present(ierr)) then
          ierr = -1
          return
        else
          write(6, "('[Error]: interp_bilinear: Latitude mismatch!')")
          write(6, "('[Debug]: interp_bilinear: ', F15.10)") y2(j)
          stop
        end if
      end if
      tmp1 = y1(j1(j))
      tmp2 = y1(j2(j))
      ywgt1(j) = (tmp2-y2(j))/(tmp2-tmp1)
      ywgt2(j) = (y2(j)-tmp1)/(tmp2-tmp1)
    end do

    do j = 1, ny2
      do i = 1, nx2
        tmp1 = data1(i1(i),j1(j))*xwgt1(i)+data1(i2(i),j1(j))*xwgt2(i)
        tmp2 = data1(i1(i),j2(j))*xwgt1(i)+data1(i2(i),j2(j))*xwgt2(i)
        data2(i,j) = tmp1*ywgt1(j)+tmp2*ywgt2(j)
      end do
    end do

  end subroutine interp_bilinear

  subroutine interp_linear(x1, data1, x2, data2, allow_extrap, ierr)

    real(8), intent(in) :: x1(:), data1(:), x2(:)
    real(8), intent(out) :: data2(:)
    logical, intent(in) :: allow_extrap
    integer, intent(out), optional :: ierr

    integer nx1, nx2
    integer size(1), i, ii
    integer, allocatable :: i1(:), i2(:)
    real(8), allocatable :: wgt1(:), wgt2(:)
    real(8) tmp1, tmp2

    size = shape(x1); nx1 = size(1)
    size = shape(x2); nx2 = size(1)

    allocate(i1(nx2), i2(nx2))
    allocate(wgt1(nx2), wgt2(nx2))

    do i = 1, nx2
      do ii = 1, nx1-1
        if ((x2(i) >= x1(ii) .and. x2(i) < x1(ii+1))) then
          i1(i) = ii
          i2(i) = ii+1
          exit
        else if (x2(i) < x1(1)) then
          if (allow_extrap) then
            i1(i) = 1
            i2(i) = 2
          else
            i1(i) = 1
            i2(i) = 1
          end if
          exit
        else if (x2(i) >= x1(nx1)) then
          if (allow_extrap) then
            i1(i) = nx1-1
            i2(i) = nx1
          else
            i1(i) = nx1
            i2(i) = nx1
          end if
          exit
        end if
      end do
      if (i1(i) == i2(i)) then
        wgt1(i) = 0.5d0
        wgt2(i) = 0.5d0
      else
        tmp1 = x1(i1(i))
        tmp2 = x1(i2(i))
        wgt1(i) = (tmp2-x2(i))/(tmp2-tmp1)
        wgt2(i) = (x2(i)-tmp1)/(tmp2-tmp1)
      end if
    end do

    do i = 1, nx2
      data2(i) = data1(i1(i))*wgt1(i)+data1(i2(i))*wgt2(i)
    end do

  end subroutine interp_linear

  subroutine interp_log_linear(x1, data1, x2, data2, allow_extrap, ierr)

    real(8), intent(in) :: x1(:), data1(:), x2(:)
    real(8), intent(out) :: data2(:)
    logical, intent(in) :: allow_extrap
    integer, intent(out), optional :: ierr

    integer nx1, nx2
    integer size(1), i, ii
    integer, allocatable :: i1(:), i2(:)
    real(8), allocatable :: wgt1(:), wgt2(:)
    real(8) tmp1, tmp2
    logical out_bound

    if (any(x1 == 0.0) .or. any(x2 == 0.0)) then
      if (present(ierr)) then
        ierr = -1
        return
      else
        write(6, "('[Error]: interp_log_linear: Input coordinate equals zero!')")
        stop
      end if
    end if

    size = shape(x1); nx1 = size(1)
    size = shape(x2); nx2 = size(1)

    allocate(i1(nx2), i2(nx2))
    allocate(wgt1(nx2), wgt2(nx2))

    do i = 1, nx2
      do ii = 1, nx1-1
        if ((x2(i) >= x1(ii) .and. x2(i) < x1(ii+1))) then
          out_bound = .false.
          i1(i) = ii
          i2(i) = ii+1
          exit
        else if (x2(i) < x1(1)) then
          out_bound = .true.
          if (allow_extrap) then
            i1(i) = 1
            i2(i) = 2
          else
            i1(i) = 1
            i2(i) = 1
          end if
          exit
        else if (x2(i) >= x1(nx1)) then
          out_bound = .true.
          if (allow_extrap) then
            i1(i) = nx1-1
            i2(i) = nx1
          else
            i1(i) = nx1
            i2(i) = nx1
          end if
          exit
        end if
      end do
      if (i1(i) == i2(i)) then
        wgt1(i) = 0.5d0
        wgt2(i) = 0.5d0
      else
        tmp1 = x1(i1(i))
        tmp2 = x1(i2(i))
        wgt1(i) = log(tmp2/x2(i))/log(tmp2/tmp1)
        wgt2(i) = log(x2(i)/tmp1)/log(tmp2/tmp1)
      end if
      if (out_bound) then
        wgt1(i) = missing_value
      end if
    end do

    do i = 1, nx2
      if (wgt1(i) == missing_value) then
        data2(i) = missing_value
      else
        data2(i) = data1(i1(i))*wgt1(i)+data1(i2(i))*wgt2(i)
      end if
    end do

  end subroutine interp_log_linear

end module interp_mod
