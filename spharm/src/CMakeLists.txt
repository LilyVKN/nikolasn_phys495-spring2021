set(ALLSRC dpmlcos_full.f90 dpmlcos.f90 dyml_full.f90 dyml.f90 pmlcos_full.f90
    pmlcos.f90 yml_full.f90 yml.f90)
set(MUPSRC local_expansion.f90 mml.f90 multipole.f90 quadrupole.f90)

set(SOURCES ${ALLSRC} ${MUPSRC})

add_library(${SPHARMLIB} ${SOURCES})
set_target_properties(
    ${SPHARMLIB} PROPERTIES 
    VERSION ${FMM_VERSION}
    SOVERSION ${FMM_MAJOR_VERSION})
fmm_install_library(${SPHARMLIB})