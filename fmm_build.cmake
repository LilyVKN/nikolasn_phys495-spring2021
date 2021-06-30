cmake_minimum_required(VERSION 3.2)

if(UNIX)
    if(NOT compiler)
        set(compiler gcc)
    endif()
    if(NOT c_compiler)
        set(compiler gcc)
    endif()
    if(NOT full_compiler)
        set(full_compiler g++)
    endif()
else()
    message(FATAL_ERROR "only UNIX systems are currently supported")
endif()

if(EXISTS "/proc/cpuinfo")
    set(parallel 1)
    file(STRINGS "/proc/cpuinfo" CPUINFO)
    foreach(line ${CPUINFO})
        if("${line}" MATCHES processor)
            math(EXPR parallel "${parallel} + 1")
        endif()
    endforeach()
endif()

find_program(HOSTNAME NAMES hostname)
find_program(UNAME NAMES uname)

execute_process(${HOSTNAME} OUTPUT_VARIABLE hostname)
string(REGEX REPLACE "[/\\\\+<> #]" "-" hostname "${hostname}")
message("HOSTNAME: ${hostname}")
if(NOT DEFINED parallel)
    set(parallel 1)
endif()

find_package(Git REQUIRED)

set(CTEST_GIT_COMMAND       ${GIT_EXECUTABLE})
set(CTeST_UPDATE_COMMAND    ${GIT_EXECUTABLE})
macro(getuname name flag)
    execute_process(COMMAND "${UNAME}" "${flag}" OUTPUT_VARIABLE "${name}")
    string(REGEX REPLACE "[/\\\\+<> #]" "-" "${name}" "${${name}}")
    string(REGEX REPLACE 
        "^(......|.....|....|...|..|.).*" "\\1" "${name}" "${${name}}")
endmacro()

getuname(osname -s)
getuname(osver  -v)
getuname(osrel  -r)
getuname(cpu    -m)

set(BUILDNAME "${osname}${osver}${osrel}${cpu}-${compiler}")
message("BUILDNAME: ${BUILDNAME}")

set(CTEST_DIR_NAME "fmmGIT")

set(CTEST_DASHBOARD_ROOT "$ENV{HOME}/Dashboards/MyTests-${BUILDNAME}")
set(CTEST_SITE "${hostname}")
set(CTEST_BUILD_NAME ${BUILDNAME})
set(CTEST_TEST_TIMEOUT "36000")

if(NOT DEFINED dashboard_git_url)
    set(dashboard_git_url 
        "https://github.com/NikolasVKN/nikolasn_phys495-spring2021")
endif()
if(NOT DEFINED dashboard_git_branch)
    set(dashboard_git_branch master)
endif()

if(NOT EXISTS "${CTEST_DASHBOARD_ROOT}/${CTEST_DIR_NAME}")
    set(CTEST_CHECKOUT_COMMAND "\"${CTeST_UPDATE_COMMAND}\" 
        clone ${dashboard_git_url} ${CTEST_DIR_NAME}")
endif()

if(NOT DEFINED CTEST_CMAKE_GENERATOR)
    set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()
set(CTEST_PROJECT_NAME          "FMM Astrophysical")
set(CTEST_BUILD_CONFIGURATION   "release")

set(ENV{DISPLAY}                "")
if(CTEST_CMAKE_GENERATOR MATCHES Makefiles)
    set(ENV{CC}                 "${c_compiler}")
    set(ENV{FC}                 "${f_compiler}")
    set(ENV{CXX}                "${full_compiler}")
endif()

make_directory("${CTEST_DASHBOARD_ROOT}")
set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_DIR_NAME}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}-${CTEST_BUILD_NAME}")
set(CTEST_NOTES_FILES "${CTEST_NOTES_FILES}" "${CMAKE_CURRENT_LIST_FILE}")

if(parallel GREATER 1)
    if(NOT CTEST_BUILD_COMMAND)
        set(CTEST_BUILD_COMMAND "make -j${parallel} -i")
    endif()

    message("Use parallel build")
    message("CTEST_BUILD_COMMAND: ${CTEST_BUILD_COMMAND}")
    message("CTEST_CONFIGURE_COMMAND: ${CTEST_CONFIGURE_COMMAND}")
endif()

set(CACHE_CONTENTS "
SITE:STRING=${hostname}
BUILDNAME:STRING=${BUILDNAME}
DART_ROOT:PATH=
GITCOMMAND:FILEPATH=${CTEST_UPDATE_COMMAND}
DROP_METHOD:STRING=https
DART_TESTING_TIMEOUT:STRING=${CTEST_TEST_TIMEOUT}
#Set build type to use optimized build
CMAKE_BUILD_TYPE:STRING=Release
")

message("Remove binary directory...")
ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")

message("CTest Directory: ${CTEST_DASHBOARD_ROOT}")
message("Initial checkout: ${CTEST_CVS_CHECKOUT}")
message("Initial cmake: ${CTEST_CMAKE_COMMAND}")
message("CTest command: ${CTEST_COMMAND}")

file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "${CACHE_CONTENTS}")

if(NOT DEFINED dashboard_model)
    set(dashboard_model nightly)
endif()
if(NOT "${dashboard_model}" MATCHES "^(Nightly|Experimental|Continuous")
    message(FATAL_ERROR "dashboard_model must be 
        Nightly, Experimental, or Continuous"
endif()

message("Start dashboard...")
ctest_start(${Dashboard_model})
message("  Update")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
message("  Configure")
ctest_update(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("read custom files after configure")
ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
message("  Build")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("  Test")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("  Submit")
ctest_submit(RETURN_VALUE res)
message("  All done")