ROOT_GENERATE_DICTIONARY(nTracksDict
    nTracks.h
    LINKDEF LinkDef.h
)

add_library(nTracks SHARED
    nTracks.cxx
    nTracksDict.cxx
)


include_directories(
    ${SNDSW_ROOT}/include
    ${ROOTSYS}/include
    ${FAIRSHIP}/include
    ${FAIRROOT_ROOT}/include
    ${FAIRLOGGER_ROOT}/include
    ${VMC_ROOT}/include/vmc
    ${FAIRROOT_ROOT}/include
    ${FMT_ROOT}/include
    ${CMAKE_SOURCE_DIR}/muonFluxUtils
    ${CMAKE_SOURCE_DIR}/nTracks
)

link_directories(
    ${SNDSW_ROOT}/lib
    ${FAIRSHIP}/lib
    ${ROOTSYS}/lib
)

target_link_libraries(nTracks PRIVATE
    muonFluxUtils
    ${ROOT_LIBRARIES}
)

add_executable(getNTracks nTracks.cxx)

target_include_directories(getNTracks PRIVATE ${ROOT_INCLUDE_DIRS})

target_link_libraries(getNTracks PRIVATE
    ${ROOT_LIBRARIES}
    Core
    Physics
    Tree
    ShipData
    shipLHC
    veto
    muonFluxUtils
)
