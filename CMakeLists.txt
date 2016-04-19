cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_LIBRARY_NAME ${PROJECT_NAME})
set(LT_VERSION_CURRENT 14)
set(LT_VERSION_AGE 1)
set(LT_VERSION_REVISION 0)
set(BUILD_SHARED_LIBS ON)

if (MSVC)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

# Prefix 
if (NOT PREFIX)
    set (PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

# Compilers
if (CC)
    set (CMAKE_C_COMPILER ${CC})
endif (CC)
if (CXX)
    set (CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

project(premia C CXX)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/Src/libpnl/src/include)

if (MSVC)
    set (LINK_TYPE STATIC)
else ()
    set (LINK_TYPE SHARED)
endif()


set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(premia_INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(premia_USE_RPATH ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_FILES_DIRECTORY "/CMakeFiles")


set(CMAKE_MODULE_PATH "${premia_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}")
include(CMakeMacroLibtoolFile)

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()

# Compiler Flags
if (UNIX)
    set(CMAKE_C_FLAGS_DEBUG "-Wall -g")
    set(CMAKE_C_FLAGS_RELEASE "-Wall -O2")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -fno-inline")
    set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O2")
endif (UNIX)
if (WIN32)
    add_definitions (-DMSDOS)
endif (WIN32)

if(MSVC)
    # Default is to statically link runtime.
    if(NOT MSVC_RUNTIME)
        set(MSVC_RUNTIME "static")
    endif(NOT MSVC_RUNTIME)
    if(MSVC_RUNTIME MATCHES "static")
        message(STATUS "MSVC -> forcing use of statically linked runtime.")
        foreach(variable ${variables})
            if(${variable} MATCHES "/MD")
                string(REGEX REPLACE "/MD" "/MT" ${variable} "${${variable}}")
            endif()
        endforeach()
    else()
        message(STATUS "MSVC -> forcing use of dynamically linked runtime.")
        foreach(variable ${variables})
            if(${variable} MATCHES "/MT")
                string(REGEX REPLACE "/MT" "/MD" ${variable} "${${variable}}")
            endif()
        endforeach()
    endif(MSVC_RUNTIME MATCHES "static")
endif(MSVC)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(PREMIA_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    set(PREMIA_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif (CMAKE_BUILD_TYPE MATCHES "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
    set(PREMIA_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    set(PREMIA_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif ()

# Add -lm if needed
if (UNIX)
    set(LIBS ${LIBS} m)
endif (UNIX)

# Find PNL
if (PNL_LIBRARY)
    set(LIBS ${LIBS} ${PNL_LIBRARY})
endif (PNL_LIBRARY)

if (PNL_LIBRARY AND (NOT PNL_INCLUDES))
    message (FATAL_ERROR "--> If you specify PNL_LIBRARY you must also specify PNL_INCLUDES")
endif ()

# Premia version
if (PREMIAVERSION)
    add_definitions(-DPremiaCurrentVersion=${PREMIAVERSION})
    set(PREMIA_C_FLAGS "${CMAKE_C_FLAGS} -DPremiaCurrentVersion=${PREMIAVERSION}")
    set(PREMIA_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPremiaCurrentVersion=${PREMIAVERSION}")
endif(PREMIAVERSION)

# Fix paths for Premia
set (PREMIADIR ${CMAKE_CURRENT_SOURCE_DIR})
set (PATH_SRC_PREMIA "${CMAKE_CURRENT_SOURCE_DIR}/Src")
set (PATH_UNIX_PREMIA "${CMAKE_CURRENT_SOURCE_DIR}/Unix")
set (PATH_COMMON "${CMAKE_CURRENT_SOURCE_DIR}/Src/common")
set (PATH_OPT "${CMAKE_CURRENT_SOURCE_DIR}/Src/opt")
set (PATH_MOD "${CMAKE_CURRENT_SOURCE_DIR}/Src/mod")
set (PATH_MAN "${CMAKE_CURRENT_SOURCE_DIR}/man")
set (PATH_MAN_TEX "${CMAKE_CURRENT_SOURCE_DIR}/man/tex")
set (PATH_MAN_PDF "${CMAKE_CURRENT_SOURCE_DIR}/man/pdf_html")
set (PATH_MAN_UTILS "${CMAKE_CURRENT_SOURCE_DIR}/man/utils")
set (PREMIA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Src")
set (PREMIA_HOME_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Unix")
set (PREMIA_MAN_PDF_DIR "${CMAKE_CURRENT_SOURCE_DIR}/man/pdf_html")
set (PREMIA_MAN_TEX_DIR "${CMAKE_CURRENT_SOURCE_DIR}/man/tex")
set (PREMIA_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
if (UNIX)
    if (APPLE)
        set (LAUNCH_USE_OPEN true)
        set (PDF_VIEWER_PROG)
    else (APPLE)
        set (LAUNCH_USE_OPEN false)
        set (PDF_VIEWER_PROG acroread)
    endif (APPLE)
endif (UNIX)


# Prepare for creating libtool .la
if (UNIX)
    LIST( REMOVE_DUPLICATES LIBS )
    foreach( _lib ${LIBS} )
        # Remove path to the library and suffixes. Transformation example: libXXX-2.0.(so|dylib) => XXX-2.0
        STRING(REGEX REPLACE ".*/lib(.*)\\.(dylib|so)" "\\1" _lib_no_path "${_lib}" )
        STRING(REGEX REPLACE "(.*)/lib[^/]*\\.(dylib|so)" "\\1" _lib_path "${_lib}" )
        if (NOT ${_lib_no_path} STREQUAL ${_lib})
            set(premia_LDFLAGS "${premia_LDFLAGS} ${CMAKE_LIBRARY_PATH_FLAG}${_lib_path}" )
        endif ()
        if (NOT ${_lib_no_path} STREQUAL ${_lib})
            set(premia_LDFLAGS "${premia_LDFLAGS} ${CMAKE_LINK_LIBRARY_FLAG}${_lib_no_path}" )
        endif ()
    endforeach( _lib )
    if (APPLE)
        set(premia_LDFLAGS "${premia_LDFLAGS} -lblas -llapack" )
    endif (APPLE)
    set(premia_LDFLAGS "${premia_LDFLAGS} -lm" )
endif (UNIX)

set(PREMIA_DIRS_TO_INCLUDE
	${PROJECT_SOURCE_DIR}/Src
	${PROJECT_SOURCE_DIR}/Src/common
    )
if (NOT PNL_LIBRARY)
    set(PREMIA_DIRS_TO_INCLUDE ${PREMIA_DIRS_TO_INCLUDE}
        ${PROJECT_SOURCE_DIR}/Src/libpnl/src/include
        ${PROJECT_BINARY_DIR}/Src/libpnl/src/include
        )
    set(LINK_PNL "pnl")
else (NOT PNL_LIBRARY)
    set(PREMIA_DIRS_TO_INCLUDE ${PREMIA_DIRS_TO_INCLUDE}
        ${PNL_INCLUDES}
        )
    set(LINK_PNL)
endif (NOT PNL_LIBRARY)


add_subdirectory(Src)

if (MSVC)
    include_directories(${PROJECT_SOURCE_DIR}/Windows/include)
    add_subdirectory (Windows EXCLUDE_FROM_ALL)
endif (MSVC)

add_subdirectory (bindings EXCLUDE_FROM_ALL)

message ("Generate configure file ${CMAKE_CURRENT_BINARY_DIR}/Src/configunix.h")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Src/configunix.h.cmake 
    ${CMAKE_CURRENT_BINARY_DIR}/Src/configunix.h)


if (UNIX AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR/scripts/}")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/tex2pdfone.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/scripts/tex2pdfone @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/tex2htmlone.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/scripts/tex2htmlone @ONLY)
endif ()

# distclean 
add_custom_target (distclean
    COMMAND echo "Cleaning all generated files"
    COMMAND make clean
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/" -name 'CMakeFiles' -type d | xargs rm -rf 
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/" -name 'CTestTestfile.cmake' | xargs rm -f 
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/" -name 'cmake_install.cmake' | xargs rm -f 
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt" 
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt"
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/libpremia.la"
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/Src/libpnl/libnl.la"
    COMMAND rm -f "${CMAKE_CURRENT_BINARY_DIR}/Src/Makefile"
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/Src/common" -name 'Makefile' | xargs rm -f
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/Src/mod" -name 'Makefile' | xargs rm -f
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/Src/opt" -name 'Makefile' | xargs rm -f
    COMMAND find "${CMAKE_CURRENT_BINARY_DIR}/Src/libpnl/src" -name 'Makefile' | xargs rm -f
    )

## Output configuration stuff ##
message("C Compiler :             ${CMAKE_C_COMPILER}")
message("     FLAGS :             ${PREMIA_C_FLAGS}")
message("C++ Compiler :           ${CMAKE_CXX_COMPILER}")
message("     FLAGS   :           ${PREMIA_CXX_FLAGS}")

message("Installation directory : ${CMAKE_INSTALL_PREFIX}") 
