ROOT_GENERATE_DICTIONARY(muonFluxUtilsDict
    muonFluxUtils.h
    LINKDEF LinkDef.h
)

add_library(muonFluxUtils SHARED muonFluxUtils.cxx)

target_include_directories(muonFluxUtils PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SNDSW_ROOT}/include
    ${ROOTSYS}/include
    ${FAIRROOT_ROOT}/include
    ${FAIRLOGGER_ROOT}/include
    ${FMT_ROOT}/include
    ${VMC_ROOT}/include/vmc
    ${FAIRSHIP}/include
)

link_directories(
    ${SNDSW_ROOT}/lib
    ${FAIRSHIP}/lib
)

find_library(SHIPDATA_LIB ShipData PATHS ${FAIRSHIP}/lib NO_DEFAULT_PATH)
find_library(SHIPLHC_LIB shipLHC PATHS ${FAIRSHIP}/lib NO_DEFAULT_PATH)

if(NOT SHIPDATA_LIB)
    message(FATAL_ERROR "Cannot find ShipData library!")
endif()

if(NOT SHIPLHC_LIB)
    message(FATAL_ERROR "Cannot find shipLHC library!")
endif()

target_link_libraries(muonFluxUtils
    PRIVATE
    ${ROOT_LIBRARIES}
    ${SHIPDATA_LIB}
    ${SHIPLHC_LIB}
)
