foreach(example robertson vanpol lorenz predator_prey pendant_drop three_equation_system oregonator e5_chemical hires)
  add_executable(${example} ${example}.f90)
  target_link_libraries(${example} PRIVATE stiff3)
  add_test(NAME ${example} COMMAND ${example})
endforeach()

# The test takes too long
add_executable(ring_modulator ring_modulator.f90)
target_link_libraries(ring_modulator PRIVATE stiff3)
