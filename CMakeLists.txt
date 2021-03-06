cmake_minimum_required(VERSION 2.8)

project(FGOALS-g-cmor LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-summary=none -ffree-line-length-none")
endif ()

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

if (DEFINED ENV{UDUNITS} AND (NOT DEFINED ENV{UDUNITS_ROOT}))
  set(ENV{UDUNITS_ROOT} $ENV{UDUNITS})
endif ()
if (DEFINED ENV{UDUNITS_ROOT})
  include_directories("$ENV{UDUNITS_ROOT}/include")
  link_directories("$ENV{UDUNITS_ROOT}/lib")
else ()
  message(FATAL_ERROR "Unable to find UDUNITS library!")
endif ()

if (DEFINED ENV{LIBUUID} AND (NOT DEFINED ENV{LIBUUID_ROOT}))
  set(ENV{LIBUUID_ROOT} $ENV{LIBUUID})
endif ()
if (DEFINED ENV{LIBUUID_ROOT})
  include_directories("$ENV{LIBUUID_ROOT}/include")
  link_directories("$ENV{LIBUUID_ROOT}/lib")
else ()
  message(FATAL_ERROR "Unable to find LIBUUID library!")
endif ()

if (DEFINED ENV{CMOR} AND (NOT DEFINED ENV{CMOR_ROOT}))
  set(ENV{CMOR_ROOT} $ENV{CMOR})
endif ()
if (DEFINED ENV{CMOR_ROOT})
  include_directories("$ENV{CMOR_ROOT}/include")
  link_directories("$ENV{CMOR_ROOT}/lib")
else ()
  message(FATAL_ERROR "Unable to find CMOR library!")
endif ()

if (DEFINED ENV{JSON_C} AND (NOT DEFINED ENV{JSON_C_ROOT}))
  set(ENV{JSON_C_ROOT} $ENV{JSON_C})
endif ()
if (DEFINED ENV{JSON_C_ROOT})
  include_directories("$ENV{JSON_C_ROOT}/include")
  link_directories("$ENV{JSON_C_ROOT}/lib")
else ()
  message(FATAL_ERROR "Unable to find JSON_C library!")
endif ()

if (DEFINED ENV{JSON_FORTRAN} AND (NOT DEFINED ENV{JSON_FORTRAN_ROOT}))
  set(ENV{JSON_FORTRAN_ROOT} $ENV{JSON_FORTRAN})
endif ()
if (DEFINED ENV{JSON_FORTRAN_ROOT})
  include_directories("$ENV{JSON_FORTRAN_ROOT}/include")
  link_directories("$ENV{JSON_FORTRAN_ROOT}/lib")
  link_directories("$ENV{JSON_FORTRAN_ROOT}/lib64")
else ()
  message(FATAL_ERROR "Unable to find JSON_FORTRAN library!")
endif ()

add_subdirectory(lib/container)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)
add_subdirectory(lib/datetime)
add_subdirectory(lib/string)

set(source_files
  src/const_mod.F90
  src/namelist_mod.F90
  src/calc_gamil_lat_bnds.F90
  src/interp_mod.F90
  src/gamil_reader_mod.F90
  src/cmor_fgoals_g_mod.F90
  src/convert_fgoals_g.F90
)

add_executable(convert_fgoals_g.exe ${source_files})
target_link_libraries(convert_fgoals_g.exe fortran_datetime fortran_string flogger fiona jsonfortran netcdf netcdff udunits2 uuid json-c cmor)
