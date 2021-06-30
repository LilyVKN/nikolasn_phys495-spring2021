set(SOURCES)

include(./naive/CMakeLists.txt)
include(./fmm/CMakeLists.txt)

add_library(${FMMLIB} ${SOURCES})
set_target_properties(
    ${FMMLIB} PROPERTIES
    VERSION ${FMM_VERSION}
    SOVERSION ${FMM_MAJOR_VERSION}
)

list(REMOVE_DUPLICATES SOURCES)

target_link_libraries(${FMMLIB} PRIVATE ${SPHARM_LIBRARIES})
fmm_install_library(${FMMLIB})