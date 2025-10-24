# CMake build for canopy-app (Fortran)

cmake_minimum_required(VERSION 3.15)
project(canopy-app LANGUAGES Fortran)

# Option for NetCDF
option(USE_NETCDF "Build with NetCDF support" ON)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags for debug levels
set(CANOPY_DEBUG_LEVEL 0 CACHE STRING "Debug level: 0=Release, 1=Debug, 2=Extra Debug")

if(CANOPY_DEBUG_LEVEL EQUAL 0)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
elseif(CANOPY_DEBUG_LEVEL EQUAL 1)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Og -Wall -Wextra -Wconversion -pedantic -frecursive -fcheck=bounds -fall-intrinsics -fmax-errors=5 -std=f2003")
elseif(CANOPY_DEBUG_LEVEL EQUAL 2)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Og -Wall -Wextra -Wconversion -pedantic -frecursive -fcheck=all -fall-intrinsics -fmax-errors=0 -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999 -std=f2003")
else()
  message(FATAL_ERROR "Invalid CANOPY_DEBUG_LEVEL: ${CANOPY_DEBUG_LEVEL}")
endif()

# Find NetCDF if requested
if(USE_NETCDF)
  find_program(NF_CONFIG nf-config)
  if(NF_CONFIG)
    execute_process(COMMAND ${NF_CONFIG} --flibs OUTPUT_VARIABLE NETCDF_FLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${NF_CONFIG} --includedir OUTPUT_VARIABLE NETCDF_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(NETCDF_DEFS "-DNETCDF")
    include_directories(${NETCDF_INCLUDE})
    set(EXTRA_LIBS ${NETCDF_FLIBS})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_DEFS}")
  else()
    message(FATAL_ERROR "nf-config not found but USE_NETCDF=ON")
  endif()
endif()

# Group common source files
set(COMMON_SRC_FILES
  canopy_const_mod.F90
  canopy_date_mod.F90
  canopy_coord_mod.F90
  canopy_canopts_mod.F90
  canopy_canmet_mod.F90
  canopy_canvars_mod.F90
  canopy_utils_mod.F90
  canopy_files_mod.F90
  canopy_readnml.F90
  canopy_alloc.F90
  canopy_init.F90
  canopy_aero_ddep_mod.F90
  canopy_txt_io_mod.F90
  canopy_dxcalc_mod.F90
  canopy_profile_mod.F90
  canopy_var3din_mod.F90
  canopy_phot_mod.F90
  canopy_rad_mod.F90
  canopy_tleaf_mod.F90
  canopy_wind_mod.F90
  canopy_fire_mod.F90
  canopy_eddy_mod.F90
  canopy_bioparm_mod.F90
  canopy_bioemi_mod.F90
  canopy_drydep_mod.F90
  canopy_calcs.F90
  canopy_write_txt.F90
  canopy_dealloc.F90
)

# NetCDF-related sources
set(NC_SRC_FILES
  canopy_ncf_io_mod.F90
  canopy_check_input.F90
)

# Main program sources
set(MAIN_SRC_FILES
  canopy_app.F90
  canopy_read_txt.F90
)

# Remove NetCDF sources if not using NetCDF
if(NOT USE_NETCDF)
  set(NC_SRC_FILES)
endif()

# Create a static library for common code
add_library(canopy_common STATIC ${COMMON_SRC_FILES})

# Create a static library for NetCDF code if enabled
if(NC_SRC_FILES)
  add_library(canopy_netcdf STATIC ${NC_SRC_FILES})
  target_link_libraries(canopy_netcdf PRIVATE canopy_common)
endif()

# Executable output in parent dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Build the main executable, linking to the libraries
add_executable(canopy ${MAIN_SRC_FILES})
target_link_libraries(canopy PRIVATE canopy_common)
if(NC_SRC_FILES)
  target_link_libraries(canopy PRIVATE canopy_netcdf)
endif()
target_link_libraries(canopy PRIVATE ${EXTRA_LIBS})

# Install targets
install(TARGETS canopy RUNTIME DESTINATION bin)
install(TARGETS canopy_common ARCHIVE DESTINATION lib)
if(NC_SRC_FILES)
  install(TARGETS canopy_netcdf ARCHIVE DESTINATION lib)
endif()
