include(CheckLanguage)

macro(qd3_check_fortran)
  set(QD_HAVE_FORTRAN OFF)
  set(QD_FORTRAN_STATUS "disabled")
  set(QD_FORTRAN_MODULE_DIR "${CMAKE_BINARY_DIR}/fortran/modules")
  set(QD_FORTRAN_MODULE_FLAG "")
  set(QD_FORTRAN_ETIME "etime")

  if(NOT QD_BUILD_FORTRAN STREQUAL "OFF")
    check_language(Fortran)
    if(NOT CMAKE_Fortran_COMPILER)
      if(QD_BUILD_FORTRAN STREQUAL "ON")
        message(FATAL_ERROR "QD_BUILD_FORTRAN=ON requested, but no Fortran compiler was found")
      endif()
      message(WARNING "No Fortran compiler found; Fortran interfaces are disabled")
    else()
      enable_language(Fortran)
      set(CMAKE_Fortran_MODULE_DIRECTORY "${QD_FORTRAN_MODULE_DIR}")
      file(MAKE_DIRECTORY "${QD_FORTRAN_MODULE_DIR}")

      include(FortranCInterface)
      FortranCInterface_VERIFY(CXX)
      FortranCInterface_HEADER(
        "${CMAKE_BINARY_DIR}/qd_fortran_mangle.h"
        MACRO_NAMESPACE "FC_"
      )

      if(CMAKE_Fortran_MODDIR_FLAG)
        set(QD_FORTRAN_MODULE_FLAG "${CMAKE_Fortran_MODDIR_FLAG}")
      else()
        set(QD_FORTRAN_MODULE_FLAG "-I")
      endif()

      set(QD_HAVE_FORTRAN ON)
      set(QD_FORTRAN_STATUS "enabled")
      set(QD_FORTRAN_ETIME "etime")
    endif()
  endif()
endmacro()
