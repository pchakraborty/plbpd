# Run plbpd
execute_process(
  COMMAND ${CMAKE_BINARY_DIR}/../plbpd
  RESULT_VARIABLE plbpd_retcode
  )
if(NOT "${plbpd_retcode}" STREQUAL "0")
  message(FATAL_ERROR "plbpd failed")
endif()

# # Compare plbpd output with baseline output
# execute_process(
#   COMMAND $ENV{HDF5_DIR}/../../bin/h5diff FinalState.h5 ${CMAKE_CURRENT_LIST_DIR}/FinalState.h5
#   RESULT_VARIABLE cmp_retcode
#   )
# if(NOT "${cmp_retcode}" STREQUAL "0")
#   message(FATAL_ERROR "comparing FinalState.h5 failed")
# endif()
