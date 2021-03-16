# - Try to find Cuba
# Once done this will define
#  Cuba_FOUND - System has Cuba
#  Cuba_INCLUDE_DIRS - The Cuba include directories
#  Cuba_LIBRARIES - The libraries needed to use Cuba


find_path(Cuba_INCLUDE_DIR cuba.h
          HINTS /usr/include
          PATH_SUFFIXES Cuba )

find_library(Cuba_LIBRARY NAMES cuba Cuba
             HINTS /usr/lib64 )

set(Cuba_LIBRARIES ${Cuba_LIBRARY} )
set(Cuba_INCLUDE_DIRS ${Cuba_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set Cuba_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Cuba  DEFAULT_MSG
                                  Cuba_LIBRARY Cuba_INCLUDE_DIR)

mark_as_advanced(Cuba_INCLUDE_DIR Cuba_LIBRARY )
