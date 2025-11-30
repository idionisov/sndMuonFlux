ROOT_GENERATE_DICTIONARY(trkeffDict
    trkeffUtils.h
    trkeffTC.h
    LINKDEF LinkDef.h
)

add_library(trkeffUtils SHARED
    trkeffTC.cxx
    trkeffUtils.cxx
    histograms.cxx
    trkeffDict.cxx
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
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/muonFluxUtils
    ${CMAKE_SOURCE_DIR}/trkeffUtils
)

link_directories(
    ${SNDSW_ROOT}/lib
    ${FAIRSHIP}/lib
)

find_library(SHIPDATA_LIB ShipData PATHS ${FAIRSHIP}/lib)
find_library(SHIPLHC_LIB ShipData PATHS ${SNDSW_ROOT}/lib)
target_link_libraries(trkeffUtils PRIVATE ${SHIPDATA_LIB} ${SHIPLHC_LIB} muonFluxUtils)

add_executable(getTrkeffTC trkeffTC.cxx histograms.cxx)
target_link_libraries(getTrkeffTC PRIVATE
    ${ROOT_LIBRARIES}
    ROOT::ROOTTPython
    ROOT::Core
    ROOT::Hist
    ROOT::MathCore
    ShipData
    shipLHC
    veto
    muonFluxUtils
    trkeffUtils
)
