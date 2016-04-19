include_directories(
    ${PROJECT_SOURCE_DIR}/Src/opt
    ${PROJECT_SOURCE_DIR}/Src/mod
    ${PROJECT_BINARY_DIR}/Src
    ${PREMIA_DIRS_TO_INCLUDE}
    )

if (NOT PNL_LIBRARY)
    option(USE_MPI "Compiles the MPI bindings for PNL. Default = on" OFF)
    set (PREFIX ${CMAKE_INSTALL_PREFIX})
    set (USE_MPI OFF)
    set (SUB_CMAKE_CALL true)
    add_subdirectory (libpnl)
endif (NOT PNL_LIBRARY)

set(CMAKE_POSITION_INDEPENDENT_CODE ON) 

add_subdirectory (common)
add_subdirectory (mod)
add_subdirectory (opt)

set(CMAKE_POSITION_INDEPENDENT_CODE ON) 

SET(premia_LIB_SRCS
    premia_obj.c
    )

set (OBJECT_LIBS
    libHawkes_Intensity
    libaffine3d
    libalsabr11d
    libalsabr21d
    libbergomi2d
    libbergomirev2d
    libbharchiarella1d
    libblack_cox_extended
    libblackkarasinski1d
    libbns
    libbs1d
    libbs2d
    libbscir2d
    libbsdisdiv1d
    libbsnd
    libcallable
    libcev1d
    libcgmy1d
    libcir1d
    libcir2d
    libcirpp1d
    libcirpp2d
    libcommon
    libcopula
    libdoublehes1d
    libdoublim
    libdps
    libdup1d
    libdynamic
    libexoi
    libfps1d
    libfps2d
    libgarch1d
    libhes1d
    libhescir1d
    libhhw4d
    libhk1d
    libhullwhite1d
    libhullwhite1dgeneralized
    libhullwhite2d
    libhw1d
    libinflation_lmm_heston1d
    libjarrowyildirim1d
    libjump1d
    libkou1d
    liblibor_affine_cir1d
    liblibor_affine_gou1d
    liblim
    liblimdisc
    liblmm1d
    liblmm1d_cgmy
    liblmm_heston1d
    liblmm_jump1d
    liblmm_stochvol_piterbarg
    liblocal_vol
    liblocvolhw1d
    liblrshjm1d
    libmath
    libmer1d
    libmerhes1d
    libmrc30d
    libnig1d
    libnig1fact1d
    libnonpar1d
    libou1d
    libpad
    libpurejump1d
    libqtsm2d
    librskou1d
    librstemperedstable1d
    libschwartz
    libschwartztrolle
    libscott1d
    libsg1d
    libstatic_merton
    libstd
    libstd2d
    libstd2dg
    libstda
    libstdc
    libstdf
    libstdg
    libstdi
    libstdnd
    libstdndc
    libstdr
    libstdx
    libstein1d
    libtemperedstable1d
    libtimehes1d
    libuvm1d
    libvariancegamma1d
    libvariancegamma2d
    libvarswap3d
    libvasicek1d
    libvol
    libwishart2d
    )

foreach(lib ${OBJECT_LIBS})
    set (premia_LIB_SRCS ${premia_LIB_SRCS} $<TARGET_OBJECTS:${lib}>)
endforeach ()

add_library(premia ${LINK_TYPE} ${premia_LIB_SRCS})
add_executable(PremiaTerm premia.c)

target_link_libraries(premia ${LIBS} ${LINK_PNL})
target_link_libraries(PremiaTerm ${LIBS} premia ${LINK_PNL})

if (MINGW)
    set_target_properties(premia PROPERTIES 
        LINK_FLAGS "-Wl,--output-def,${CMAKE_CURRENT_BINARY_DIR}/libpremia.def")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libpremia.def DESTINATION lib)
endif (MINGW)


if (UNIX)
    set_target_properties(premia PROPERTIES 
        LT_VERSION_CURRENT ${LT_VERSION_CURRENT}
        LT_VERSION_AGE ${LT_VERSION_AGE}
        LT_VERSION_REVISION ${LT_VERSION_REVISION}
        LT_DEPENDENCY_LIBS ${premia_LDFLAGS}
        SOVERSION ${LT_VERSION_CURRENT}
        VERSION ${LT_VERSION_CURRENT}.${LT_VERSION_AGE}.${LT_VERSION_REVISION}
        )
    CREATE_LIBTOOL_FILE(premia /lib)
endif (UNIX)

install(TARGETS premia DESTINATION lib)
install(TARGETS PremiaTerm DESTINATION bin)
install(FILES 
    premia_obj.h config.h
    DESTINATION include
    )
