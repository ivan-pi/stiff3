add_library(stiff3
  stiff3_linalg.f90
  stiff3_solver.f90
)

target_link_libraries(stiff3 PUBLIC ${LAPACK_LIBRARIES})

target_include_directories(stiff3 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
