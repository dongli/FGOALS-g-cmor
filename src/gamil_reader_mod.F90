module gamil_reader_mod

  use fiona
  use interp_mod

  implicit none

  private

  public gamil_reader_open
  public gamil_reader_get_grids
  public gamil_reader_get_var
  public gamil_reader_get_att
  public gamil_reader_close
  public gamil_reader_final

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

  real(8), parameter :: ptop = 219.4
  real(8), allocatable :: sigma(:)
  real(8), allocatable :: p(:)
  real(8), allocatable :: ps(:,:)
  real(8), allocatable :: array_on_gamil_levels(:,:,:)

contains

  subroutine gamil_reader_open(file_path)

    character(*), intent(in) :: file_path

    call io_create_dataset('gamil', file_path=file_path, mode='input')
    call io_get_dim('gamil', 'lon', size=num_lon)
    call io_get_dim('gamil', 'lat', size=num_lat)
    call io_get_dim('gamil', 'lev', size=num_lev)
    call io_start_input('gamil')

    if (.not. allocated(sigma)) allocate(sigma(num_lev))
    if (.not. allocated(p)) allocate(p(num_lev))
    if (.not. allocated(ps)) allocate(ps(num_lon,num_lat))
    if (.not. allocated(array_on_gamil_levels)) allocate(array_on_gamil_levels(num_lon,num_lat,num_lev))

    call io_input('gamil', 'lev', sigma)

  end subroutine gamil_reader_open

  subroutine gamil_reader_get_grids(lon, lon_bnds, lat, lat_bnds, lev, lev_bnds)

    real(8), intent(out), allocatable :: lon(:)
    real(8), intent(out), allocatable :: lon_bnds(:)
    real(8), intent(out), allocatable :: lat(:)
    real(8), intent(out), allocatable :: lat_bnds(:)
    real(8), intent(out), allocatable :: lev(:)
    real(8), intent(out), allocatable :: lev_bnds(:)

    real(8) dlon
    integer i

    if (allocated(lon)) deallocate(lon); allocate(lon(num_lon))
    if (allocated(lon_bnds)) deallocate(lon_bnds); allocate(lon_bnds(num_lon + 1))
    if (allocated(lat)) deallocate(lat); allocate(lat(num_lat))
    if (allocated(lat_bnds)) deallocate(lat_bnds); allocate(lat_bnds(num_lat + 1))
    if (allocated(lev)) deallocate(lev); allocate(lev(num_lev))
    if (allocated(lev_bnds)) deallocate(lev_bnds); allocate(lev_bnds(num_lev + 1))

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

  subroutine gamil_reader_get_var_0d(var_name, value)

    character(*), intent(in) :: var_name
    real(8), intent(out) :: value

    call io_input('gamil', var_name, value)

  end subroutine gamil_reader_get_var_0d

  subroutine gamil_reader_get_var_1d(var_name, array)

    character(*), intent(in) :: var_name
    real(8), intent(out) :: array(:)

    call io_input('gamil', var_name, array)

  end subroutine gamil_reader_get_var_1d

  subroutine gamil_reader_get_var_2d(var_name, array)

    character(*), intent(in) :: var_name
    real(8), intent(out) :: array(:,:)

    call io_input('gamil', var_name, array)

    if (var_name == 'PS') then
      ps = array
    end if

  end subroutine gamil_reader_get_var_2d

  subroutine gamil_reader_get_var_3d(var_name, array, plev, use_log_linear)

    character(*), intent(in) :: var_name
    real(8), intent(out) :: array(:,:,:)
    real(8), intent(in), optional :: plev(:)
    logical, intent(in), optional :: use_log_linear

    integer i, j

    if (present(plev)) then
      call io_input('gamil', var_name, array_on_gamil_levels)
      do j = 1, num_lat
        do i = 1, num_lon
          p = sigma * (ps(i,j) - ptop) + ptop
          if (present(use_log_linear) .and. use_log_linear) then
            call interp_log_linear(p, array_on_gamil_levels(i,j,:), plev, array(i,j,:), allow_extrap=.false.)
          else
            call interp_linear(p, array_on_gamil_levels(i,j,:), plev, array(i,j,:), allow_extrap=.false.)
          end if
        end do
      end do
    else
      call io_input('gamil', var_name, array)
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

    if (allocated(sigma)) deallocate(sigma)
    if (allocated(p)) deallocate(p)
    if (allocated(ps)) deallocate(ps)
    if (allocated(array_on_gamil_levels)) deallocate(array_on_gamil_levels)

  end subroutine gamil_reader_final

end module gamil_reader_mod
