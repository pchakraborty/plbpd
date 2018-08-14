if(NOT DEFINED ENV{TBBROOT})
  message(FATAL_ERROR "ERROR: env var TBBROOT is not defined")
endif(NOT DEFINED ENV{TBBROOT})

include_directories($ENV{TBBROOT}/include)
if(APPLE)
  set(TBBLibName libtbb.dylib)
elseif(UNIX)
  set(TBBLibName libtbb.so)
else()
  message(FATAL_ERROR "Platform not recognized")
endif()  
find_library(
  LIBTBB
  NAMES ${TBBLibName}
  PATHS $ENV{TBBROOT}/lib $ENV{TBBROOT}/lib/intel64_lin/gcc4.7)
get_filename_component(TBB_LIB_DIR ${LIBTBB} DIRECTORY)
message(STATUS "TBB: ${LIBTBB}")
