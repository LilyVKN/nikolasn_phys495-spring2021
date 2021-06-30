add_subdirectory(src)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/spharm.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${SPHARMLIB}.cmake @ONLY)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${SPHARMLIB}.cmake
    DESTINATION ${PKG_CONFIG_DIR}
    COMPONENT Development)