foreach(example robertson vanpol lorenz predator_prey pendant_drop three_equation_system oregonator e5_chemical hires)
  add_executable(${example} ${example}.f90)
  target_link_libraries(${example} PRIVATE stiff3)
  add_test(NAME ${example} COMMAND ${example})
endforeach()

# The test takes too long
add_executable(ring_modulator ring_modulator.f90)
target_link_libraries(ring_modulator PRIVATE stiff3)

add_executable(caroll_driver caroll_driver.f90 caroll_test_problems.f90)
target_link_libraries(caroll_driver PRIVATE stiff3)

add_executable(hosea_driver hosea_driver.f90)
target_link_libraries(hosea_driver PRIVATE stiff3)

add_executable(fluid_bed_reactor fluid_bed_reactor.f90)
target_link_libraries(fluid_bed_reactor PRIVATE stiff3)