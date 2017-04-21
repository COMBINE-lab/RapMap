set (RAPMAP_MAIN_SRCS
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
    HitManagerSA.cpp
    EditDistance.cpp
    FastxParser.cpp
    rank9b.cpp
    edlib.cpp
    stringpiece.cc
    xxhash.c
    bit_array.c
)

set (RAPMAP_ALIGN_SRCS
)

set (RAPMAP_LIB_SRCS
)

#set (RSDICT_LIB_SRCS
#    EnumCoder.cpp
#    RSDic.cpp
#    RSDicBuilder.cpp
#)

include_directories(
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.2.6
${ZLIB_INCLUDE_DIR}
${CEREAL_INCLUDE_DIRS}
)

if (JELLYFISH_FOUND)
    include_directories(${JELLYFISH_INCLUDE_DIR})
else()
    include_directories(${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.2.6)
endif()

link_directories(
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
#${Boost_LIBRARY_DIRS}
#${TBB_LIBRARY_DIRS}
#${LAPACK_LIBRARY_DIR}
#${BLAS_LIBRARY_DIR}
)

# Set the RPATH
if (NOT APPLE)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:$ORIGIN/:$ORIGIN/../../external/install/lib")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
else()
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

  # the RPATH to be used when installing
  set(CMAKE_INSTALL_RPATH "")

  # don't add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
endif()

# Build the rapmap executable
add_executable(rapmap ${RAPMAP_MAIN_SRCS})

# our suffix array construction libraries
set (SUFFARRAY_LIB ${GAT_SOURCE_DIR}/external/install/lib/libdivsufsort.a)
set (SUFFARRAY64_LIB ${GAT_SOURCE_DIR}/external/install/lib/libdivsufsort64.a)

# Link the executable
target_link_libraries(rapmap
    ${PTHREAD_LIB}
    ${ZLIB_LIBRARY}
    ${SUFFARRAY_LIB}
    ${SUFFARRAY64_LIB}
    ${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a
    m
    ${LIBLZMA_LIBRARIES}
    ${NON_APPLECLANG_LIBS}
    ${FAST_MALLOC_LIB}
)


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

install(FILES ${GAT_SOURCE_DIR}/scripts/RunRapMap.sh 
              PERMISSIONS WORLD_EXECUTE WORLD_READ OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
              DESTINATION bin)
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
