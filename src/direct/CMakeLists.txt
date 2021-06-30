set(SUBSRCDIR ./naive/)

SET(NAIVESRC force_grav_submat.f90 part_init_default.f90
    sim_params.f90)
list(TRANSFORM NAIVESRC PREPEND ${SUBSRCDIR})

list(APPEND SOURCES ${NAIVESRC})