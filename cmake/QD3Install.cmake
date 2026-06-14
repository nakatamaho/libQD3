include(CMakePackageConfigHelpers)

function(qd3_install_targets)
  install(TARGETS ${QD3_INSTALL_TARGETS}
    EXPORT QD3Targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(FILES ${QD_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qd
  )
  install(FILES "${CMAKE_BINARY_DIR}/include/qd/qd_config.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qd
  )

  install(FILES
    "${CMAKE_SOURCE_DIR}/README"
    "${CMAKE_SOURCE_DIR}/README.md"
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
  )

  if(QD_HAVE_FORTRAN)
    install(DIRECTORY "${QD_FORTRAN_MODULE_DIR}/"
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qd
      FILES_MATCHING
        PATTERN "*.mod"
        PATTERN "*.MOD"
        PATTERN "*.smod"
        PATTERN "*.SMOD"
    )
  endif()

  install(FILES "${CMAKE_BINARY_DIR}/qd.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
  )

  install(PROGRAMS "${CMAKE_BINARY_DIR}/qd-config"
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

  configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/QD3Config.cmake.in"
    "${CMAKE_BINARY_DIR}/QD3Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QD3
  )
  write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/QD3ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  install(FILES
    "${CMAKE_BINARY_DIR}/QD3Config.cmake"
    "${CMAKE_BINARY_DIR}/QD3ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QD3
  )
  install(EXPORT QD3Targets
    NAMESPACE QD3::
    FILE QD3Targets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QD3
  )
endfunction()
