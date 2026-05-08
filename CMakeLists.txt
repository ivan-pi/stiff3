cmake_minimum_required(VERSION 3.14)

project(stiff3
  VERSION 0.2.0
  DESCRIPTION "Adaptive solver for stiff systems of ODEs"
  LANGUAGES Fortran
)

# Find LAPACK (which typically pulls in BLAS as well)
find_package(LAPACK REQUIRED)

include(CTest)

add_subdirectory(src)
add_subdirectory(example)
add_subdirectory(test)
