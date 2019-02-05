set (RAPMAP_MAIN_SRCS
    RapMap.cpp
    RapMapSAIndexer.cpp
    RapMapUtils.cpp
    RapMapSAMapper.cpp
    RapMapFileSystem.cpp
    RapMapSAIndex.cpp
    HitManager.cpp
    FastxParser.cpp
    rank9b.cpp
    stringpiece.cc
    xxhash.c
    bit_array.c
    edlib.cpp
    metro/metrohash64.cpp
    #RapMapIndex.cpp
    #RapMapMapper.cpp
    #RapMapIndexer.cpp
)

set (RAPMAP_ALIGN_SRCS
)

set (RAPMAP_LIB_SRCS
)

include_directories(
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/digestpp
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
${ZLIB_INCLUDE_DIR}
${CEREAL_INCLUDE_DIRS}
)

link_directories(
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
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

# build KSW2 library
set (KSW2PP_BASIC_LIB_SRCS
ksw2pp/kalloc.c
ksw2pp/ksw2_extd.c
ksw2pp/ksw2_extz.c
ksw2pp/ksw2_gg.c
ksw2pp/ksw2_gg2.c
ksw2pp/ksw2_gg2_sse.c
ksw2pp/KSW2Aligner.cpp
)

set (KSW2PP_ADVANCED_LIB_SRCS
ksw2pp/ksw2_extd2_sse.c
ksw2pp/ksw2_extf2_sse.c
ksw2pp/ksw2_extz2_sse.c
)

# check if we know how to do IPO
check_ipo_supported(RESULT HAS_IPOHAS_IPO)

add_library(ksw2pp_sse2 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_sse4 OBJECT ${KSW2PP_ADVANCED_LIB_SRCS})
add_library(ksw2pp_basic OBJECT ${KSW2PP_BASIC_LIB_SRCS})

set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_FLAGS "-O3 -msse2 -mno-sse4.1")
set_target_properties(ksw2pp_sse2 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;KSW_SSE2_ONLY;HAVE_KALLOC")
set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_FLAGS "-O3 -msse4.1")
set_target_properties(ksw2pp_sse4 PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
set_target_properties(ksw2pp_basic PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")

set_target_properties(ksw2pp_basic PROPERTIES INCLUDE_DIRECTORIES ${GAT_SOURCE_DIR}/include)
set_target_properties(ksw2pp_sse4 PROPERTIES INCLUDE_DIRECTORIES ${GAT_SOURCE_DIR}/include)

# Build the ksw2pp library
add_library(ksw2pp STATIC $<TARGET_OBJECTS:ksw2pp_sse2> $<TARGET_OBJECTS:ksw2pp_sse4> $<TARGET_OBJECTS:ksw2pp_basic>)
set_target_properties(ksw2pp PROPERTIES COMPILE_DEFINITIONS "KSW_CPU_DISPATCH;HAVE_KALLOC")
if(HAS_IPO)
  set_property(TARGET ksw2pp PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()



# Link the executable
target_link_libraries(rapmap
    # ${PTHREAD_LIB}
    ${ZLIB_LIBRARY}
    ${SUFFARRAY_LIB}
    ${SUFFARRAY64_LIB}
    m
    ${NON_APPLECLANG_LIBS}
    ${FAST_MALLOC_LIB}
    ${CMAKE_THREAD_LIBS_INIT}
    ${LIBRT}
    ksw2pp
    #${GAT_SOURCE_DIR}/external/install/lib/libjellyfish-2.0.a
    #${LIBLZMA_LIBRARIES}
)

##
# External dependencies of salmon_core and salmon
##
if (${FETCHED_JEMALLOC})
   add_dependencies(rapmap libjemalloc)
   add_dependencies(rapmap libjemalloc)
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
    add_test( NAME quasi_map_test COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestQuasiMap.cmake )
    add_test( NAME quasi_map_test_ph COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/TestQuasiMapPerfectHash.cmake )
