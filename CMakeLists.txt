cmake_minimum_required(VERSION 3.12)

project(FMM Fortran)

set(FMM_MAJOR_VERSION 0)
set(FMM_MINOR_VERSION 0)
set(FMM_PATCH_VERSION 0)
set(FMM_VERSION ${FMM_MAJOR_VERSION}.${FMM_MINOR_VERSION}.${FMM_PATCH_VERSION})

set(CMAKE_MODULE_PATH "${FMM_SOURCE_DIR}/CMAKE" ${CMAKE_MODULE_PATH})

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

option(BUILD_INDEX64 "Build Index-64 API libraries" OFF)
if(BUILD_INDEX64)
    set(SPHARMLIB "spharm64")
    set(FMMLIB "fmm64")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" -fdefault-integer-8)
else()
    set(SPHARMLIB "spharm")
    set(FMMLIB "fmm")
endif()

include(GNUInstallDirs)

set(FMM_INSTALL_EXPORT_NAME ${FMMLIB}-targets)

macro(fmm_install_library lib)
    install(
        TARGETS ${lib}
        EXPORT ${FMM_INSTALL_EXPORT_NAME}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
        RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries)
endmacro()

add_subdirectory(spharmpack)
add_subdirectory(src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${FMM_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${FMM_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${FMM_BINARY_DIR}/lib)

set(CPACK_PACKAGE_NAME "FMM Astrophysical")
set(CPACK_PACKAGE_VENDOR "Nikolas Nguyen")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
    "Fast Multipole Method for N-body, astrophysical simulations")