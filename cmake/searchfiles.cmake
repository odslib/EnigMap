set_property( GLOBAL PROPERTY INCLUDE_FILES_G "")
set_property( GLOBAL PROPERTY SOURCE_FILES_NM_G "")

macro(searchfiles)
  get_property(INCLUDE_FILES GLOBAL PROPERTY INCLUDE_FILES_G)
  get_property(SOURCE_FILES_NM GLOBAL PROPERTY SOURCE_FILES_NM_G)
  file(GLOB INCLUDE_FILES_A *.hpp)
  list(APPEND INCLUDE_FILES "${INCLUDE_FILES_A}")
  file(GLOB SOURCE_FILES_NM_A *.cpp)
  list(APPEND SOURCE_FILES_NM "${SOURCE_FILES_NM_A}")
  set_property( GLOBAL PROPERTY INCLUDE_FILES_G "${INCLUDE_FILES}")
  set_property( GLOBAL PROPERTY SOURCE_FILES_NM_G "${SOURCE_FILES_NM}")
endmacro()