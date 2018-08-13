# Check for HDF5 library - CXX

# The env var HDF5_DIR needs to be present
if(NOT DEFINED ENV{HDF5_DIR})
  message(FATAL_ERROR "ERROR: env var HDF5_DIR is not defined")
else(NOT DEFINED ENV{HDF5_DIR})
  message(STATUS "HDF5_DIR: $ENV{HDF5_DIR}")
endif(NOT DEFINED ENV{HDF5_DIR})

# Bypass cmake's internal FindHDF5 script
find_package(HDF5 NAMES hdf5 COMPONENTS C static)

# Check for C/CXX include directories
if(NOT DEFINED HDF5_INCLUDE_DIR)
  message(FATAL_ERROR "ERROR: var HDF5_INCLUDE_DIR is not defined")
else(NOT DEFINED HDF5_INCLUDE_DIR)
  message(STATUS "HDF5 include dir: ${HDF5_INCLUDE_DIR}")
endif(NOT DEFINED HDF5_INCLUDE_DIR)
include_directories(${HDF5_INCLUDE_DIR})

# Check for C static libraries
if(NOT DEFINED HDF5_C_STATIC_LIBRARY)
  message(FATAL_ERROR "ERROR: var HDF5_C_STATIC_LIBRARY is not defined")
else(NOT DEFINED HDF5_C_STATIC_LIBRARY)
  message(STATUS "HDF5 C static library: ${HDF5_C_STATIC_LIBRARY}")
endif(NOT DEFINED HDF5_C_STATIC_LIBRARY)
