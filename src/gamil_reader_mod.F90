module gamil_reader_mod

  use fiona
  use flogger
  use interp_mod

  implicit none

  private

  public gamil_reader_open
  public gamil_reader_get_grids
  public gamil_reader_get_var
  public gamil_reader_get_att
  public gamil_reader_close
  public gamil_reader_final

  public time_interval
  public lon
  public lon_bnds
  public lat
  public lat_bnds
  public lev
  public lev_bnds
  public ptop

  character(30), parameter :: time_interval = '10 minutes'

  integer num_lon
  integer num_lat
  integer num_lev

  interface gamil_reader_get_var
    module procedure gamil_reader_get_var_0d
    module procedure gamil_reader_get_var_1d
    module procedure gamil_reader_get_var_2d
    module procedure gamil_reader_get_var_3d
  end interface gamil_reader_get_var

  interface gamil_reader_get_att
    module procedure gamil_reader_get_var_att_s
  end interface gamil_reader_get_att

  real(8), allocatable :: lon(:)
  real(8), allocatable :: lon_bnds(:)
  real(8), allocatable :: lat(:)
  real(8), allocatable :: lat_bnds(:)
  real(8), allocatable :: lev(:)
  real(8), allocatable :: lev_bnds(:)

  real(8), parameter :: ptop = 219.4
  real(8), allocatable, dimension(:    ) :: pfull
  real(8), allocatable, dimension(:    ) :: phalf
  real(8), allocatable, dimension(:,:  ) :: ps
  real(8), allocatable, dimension(:,:  ) :: buf_2d
  real(8), allocatable, dimension(:,:,:) :: buf_3d

contains

  subroutine gamil_reader_open(file_path)

    character(*), intent(in) :: file_path

    call io_create_dataset('gamil', file_path=file_path, mode='input', mute=.true.)
    call io_get_dim('gamil', 'lon', size=num_lon)
    call io_get_dim('gamil', 'lat', size=num_lat)
    call io_get_dim('gamil', 'lev', size=num_lev)
    call io_start_input('gamil')

    if (.not. allocated(lon     )) allocate(lon     (num_lon                  ))
    if (.not. allocated(lon_bnds)) allocate(lon_bnds(num_lon+1                ))
    if (.not. allocated(lat     )) allocate(lat     (        num_lat          ))
    if (.not. allocated(lat_bnds)) allocate(lat_bnds(        num_lat+1        ))
    if (.not. allocated(lev     )) allocate(lev     (                num_lev  ))
    if (.not. allocated(lev_bnds)) allocate(lev_bnds(                num_lev+1))
    if (.not. allocated(pfull   )) allocate(pfull   (                num_lev  ))
    if (.not. allocated(phalf   )) allocate(phalf   (                num_lev+1))
    if (.not. allocated(ps      )) allocate(ps      (num_lon,num_lat          ))
    if (.not. allocated(buf_2d  )) allocate(buf_2d  (num_lon,num_lat          ))
    if (.not. allocated(buf_3d  )) allocate(buf_3d  (num_lon,num_lat,num_lev  ))

  end subroutine gamil_reader_open

  subroutine gamil_reader_get_grids()

    real(8) dlon
    integer i

    call io_input('gamil', 'lon', lon)
    call io_input('gamil', 'lat', lat)
    call io_input('gamil', 'lev', lev)
    call io_input('gamil', 'ilev', lev_bnds)

    dlon = lon(2) - lon(1)
    do i = 1, num_lon
      lon_bnds(i) = lon(i) - 0.5d0 * dlon
    end do
    lon_bnds(num_lon+1) = lon(num_lon) + 0.5d0 * dlon

    call calc_gamil_lat_bnds(num_lat, 2.0d0, lat, lat_bnds)

  end subroutine gamil_reader_get_grids

  subroutine gamil_reader_get_var_0d(var_name, value, time_step)

    character(*), intent(in   ) :: var_name
    real(8)     , intent(  out) :: value
    integer     , intent(in   ) :: time_step

    call io_input('gamil', var_name, value, time_step=time_step)

  end subroutine gamil_reader_get_var_0d

  subroutine gamil_reader_get_var_1d(var_name, array, time_step)

    character(*), intent(in   ) :: var_name
    real(8)     , intent(  out) :: array(:)
    integer     , intent(in   ) :: time_step

    call io_input('gamil', var_name, array, time_step=time_step)

  end subroutine gamil_reader_get_var_1d

  subroutine gamil_reader_get_var_2d(var_name, array, time_step)

    character(*), intent(in   ) :: var_name
    real(8)     , intent(  out) :: array(:,:)
    integer     , intent(in   ) :: time_step

    if (io_has_var('gamil', var_name)) then
      call io_input('gamil', var_name, array, time_step=time_step)
    end if

    select case (var_name)
    case ('PRECT', 'PRECC')
      array = array * 1000
    case ('PRECSC+PRECSL')
      if (io_has_var('gamil', 'PRECSC') .and. io_has_var('gamil', 'PRECSL')) then
        call io_input('gamil', 'PRECSC', array, time_step=time_step)
        call io_input('gamil', 'PRECSL', buf_2d, time_step=time_step)
        array = array + buf_2d
      else
        call log_warning('Variable PRECSC+PRECSL cannot be outputted!')
      end if
    case ('TGCLDLWP+TGCLDIWP')
      if (io_has_var('gamil', 'TGCLDLWP') .and. io_has_var('gamil', 'TGCLDIWP')) then
        call io_input('gamil', 'TGCLDLWP', array, time_step=time_step)
        call io_input('gamil', 'TGCLDIWP', buf_2d, time_step=time_step)
        array = array + buf_2d
      else
        call log_warning('Variable TGCLDLWP+TGCLDIWP cannot be outputted!')
      end if
    case ('FLNS')
      if (io_has_var('gamil', 'FLUS') .and. io_has_var('gamil', 'FLDS')) then
        call io_input('gamil', 'FLUS', array, time_step=time_step)
        call io_input('gamil', 'FLDS', buf_2d, time_step=time_step)
        array = array - buf_2d
      else
        call log_warning('Variable FLNS (rls) cannot be outputted!')
      end if
    case ('SOLLD+SOLSD')
      if (io_has_var('gamil', 'SOLLD') .and. io_has_var('gamil', 'SOLSD')) then
        call io_input('gamil', 'SOLLD', array, time_step=time_step)
        call io_input('gamil', 'SOLSD', buf_2d, time_step=time_step)
        array = array + buf_2d
      else
        call log_warning('Variable SOLLD+SOLSD cannot be outputted!')
      end if
    case ('FSNT-FLNT')
      if (io_has_var('gamil', 'FSNT') .and. io_has_var('gamil', 'FLNT')) then
        call io_input('gamil', 'FSNT', array, time_step=time_step)
        call io_input('gamil', 'FLNT', buf_2d, time_step=time_step)
        array = array - buf_2d
      else
        call log_warning('Variable FSNT-FLNT cannot be outputted!')
      end if
    end select

  end subroutine gamil_reader_get_var_2d

  subroutine gamil_reader_get_var_3d(var_name, array, time_step, plev, use_log_linear)

    character(*), intent(in   )           :: var_name
    real(8)     , intent(  out)           :: array(:,:,:)
    integer     , intent(in   )           :: time_step
    real(8)     , intent(in   ), optional :: plev(:)
    logical     , intent(in   ), optional :: use_log_linear

    integer i, j

    if (present(plev)) then
      call io_input('gamil', 'PS', ps, time_step=time_step)
      call io_input('gamil', var_name, buf_3d, time_step=time_step)
      do j = 1, num_lat
        do i = 1, num_lon
          pfull = lev * (ps(i,j) - ptop) + ptop
          if (present(use_log_linear) .and. use_log_linear) then
            call interp_log_linear(pfull, buf_3d(i,j,:), plev, array(i,j,:), allow_extrap=.false.)
          else
            call interp_linear(pfull, buf_3d(i,j,:), plev, array(i,j,:), allow_extrap=.false.)
          end if
        end do
      end do
    else if (var_name == '<pfull>') then
      call io_input('gamil', 'PS', ps, time_step=time_step)
      do j = 1, num_lat
        do i = 1, num_lon
          pfull = lev * (ps(i,j) - ptop) + ptop
          array(i,j,:) = pfull
        end do
      end do
    else if (var_name == '<phalf>') then
      call io_input('gamil', 'PS', ps, time_step=time_step)
      do j = 1, num_lat
        do i = 1, num_lon
          phalf = lev_bnds * (ps(i,j) - ptop) + ptop
          array(i,j,:) = phalf
        end do
      end do
    else
      call io_input('gamil', var_name, array, time_step=time_step)
    end if

  end subroutine gamil_reader_get_var_3d

  subroutine gamil_reader_get_var_att_s(var_name, att_name, value)

    character(*), intent(in ) :: var_name
    character(*), intent(in ) :: att_name
    character(*), intent(out) :: value

    call io_get_att('gamil', var_name, att_name, value)

  end subroutine gamil_reader_get_var_att_s

  subroutine gamil_reader_close()

    call io_end_input('gamil')

  end subroutine gamil_reader_close

  subroutine gamil_reader_final()

    if (allocated(lon     )) deallocate(lon     )
    if (allocated(lon_bnds)) deallocate(lon_bnds)
    if (allocated(lat     )) deallocate(lat     )
    if (allocated(lat_bnds)) deallocate(lat_bnds)
    if (allocated(lev     )) deallocate(lev     )
    if (allocated(lev_bnds)) deallocate(lev_bnds)
    if (allocated(pfull   )) deallocate(pfull   )
    if (allocated(phalf   )) deallocate(phalf   )
    if (allocated(ps      )) deallocate(ps      )
    if (allocated(buf_2d  )) deallocate(buf_3d  )
    if (allocated(buf_3d  )) deallocate(buf_3d  )

  end subroutine gamil_reader_final

end module gamil_reader_mod
