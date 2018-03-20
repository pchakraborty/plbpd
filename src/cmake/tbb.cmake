if(NOT DEFINED ENV{TBBROOT})
  message(FATAL_ERROR "ERROR: env var TBBROOT is not defined")
endif(NOT DEFINED ENV{TBBROOT})

include_directories($ENV{TBBROOT}/include)
set(TBBLibName libtbb.so)
find_library(LIBTBB NAMES ${TBBLibName} PATHS $ENV{TBBROOT}/lib)
get_filename_component(TBB_LIB_DIR ${LIBTBB} DIRECTORY)
message(STATUS "TBB: ${LIBTBB}")
