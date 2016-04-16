set (RAPMAP_MAIN_SRCS
    EMPHFCommon.cpp
    RapMap.cpp
    RapMapIndexer.cpp
    RapMapSAIndexer.cpp
    RapMapUtils.cpp
    RapMapMapper.cpp
    RapMapSAMapper.cpp
    RapMapFileSystem.cpp
    RapMapSAIndex.cpp
    RapMapIndex.cpp
    HitManager.cpp
    #sherwood_map.cpp
    rank9b.cpp
    #intersection.cpp
    #bitmap.cpp
    stringpiece.cc
    sais.c
    xxhash.c
    bit_array.c
)

set (RAPMAP_ALIGN_SRCS
)

set (RAPMAP_LIB_SRCS
)

set (RSDICT_LIB_SRCS
    EnumCoder.cpp
    RSDic.cpp
    RSDicBuilder.cpp
)

include_directories(
${GAT_SOURCE_DIR}/include
#${GAT_SOURCE_DIR}/include/eigen3
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.2.5
#${GAT_SOURCE_DIR}/external/install/include/emphf
${ZLIB_INCLUDE_DIR}
${CEREAL_INCLUDE_DIRS}
#${Boost_INCLUDE_DIRS}
)

if (JELLYFISH_FOUND)
    include_directories(${JELLYFISH_INCLUDE_DIR})
else()
    include_directories(${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.2.5)
endif()

link_directories(
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
#${Boost_LIBRARY_DIRS}
#${TBB_LIBRARY_DIRS}
#${LAPACK_LIBRARY_DIR}
#${BLAS_LIBRARY_DIR}
)

#message("Boost_LIBRARIES = ${Boost_LIBRARIES}")

# Set the RPATH
if (APPLE)
    ## This DOES NOT do what I / any one sane, expects.  Setting the
    ## linker path on OSX is messed up.  Just tell the user to use
    ## DYLD_FALLBACK_LIBRARY_PATH for now
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:$ORIGIN/:$ORIGIN/../../external/install/lib")
endif()

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Build the rsdic library
add_library(rsdic STATIC ${RSDICT_LIB_SRCS} )

# Build the rapmap executable
add_executable(rapmap ${RAPMAP_MAIN_SRCS})

#set_target_properties(salmon_core salmon PROPERTIES LINK_SEARCH_END_STATIC TRUE)

# our suffix array construction libraries
set (SUFFARRAY_LIB ${GAT_SOURCE_DIR}/external/install/lib/libdivsufsort.a)
set (SUFFARRAY64_LIB ${GAT_SOURCE_DIR}/external/install/lib/libdivsufsort64.a)

# Link the executable
target_link_libraries(rapmap
    rsdic
    ${PTHREAD_LIB}
    #${Boost_LIBRARIES}
    ${ZLIB_LIBRARY}
    ${SUFFARRAY_LIB}
    ${SUFFARRAY64_LIB}
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a
    m
    ${LIBLZMA_LIBRARIES}
    #${LIBSALMON_LINKER_FLAGS}
    ${NON_APPLECLANG_LIBS}
    ${FAST_MALLOC_LIB}
)

#add_dependencies(salmon libbwa)

##
#  This ensures that the salmon executable should work with or without `make install`
##
if (APPLE)
	add_custom_command(TARGET rapmap
		POST_BUILD
		COMMAND install_name_tool -add_rpath ${GAT_SOURCE_DIR}/external/install/lib rapmap
		COMMAND install_name_tool -add_rpath @executable_path/../lib rapmap
		)
endif()

##### ======================================

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "${GAT_SOURCE_DIR}" CACHE PATH "Default install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(INSTALL_LIB_DIR lib )
set(INSTALL_BIN_DIR bin )
set(INSTALL_INCLUDE_DIR include )

# install(FILES ${Boost_LIBRARIES}
# 	           DESTINATION ${INSTALL_LIB_DIR})

install(TARGETS rapmap
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
        )

    #set(POST_INSTALL_SCRIPT ${GAT_SOURCE_DIR}/cmake/PostInstall.cmake)

    #install(
    #    CODE
    #    "
    #    execute_process(COMMAND \"${CMAKE_COMMAND}\"
    #                            -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
    #                            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    #                            -P \"${POST_INSTALL_SCRIPT}\")
    #    "
    #)
    #
    #include(InstallRequiredSystemLibraries)
    #add_test( NAME salmon_read_test COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestSalmon.cmake )
