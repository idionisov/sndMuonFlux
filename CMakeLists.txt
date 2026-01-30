cmake_minimum_required(VERSION 3.23)
project(muonFlux)

find_package(ROOT REQUIRED COMPONENTS RIO Core Hist MathCore Tree ROOTTPython)
find_package(Python3 COMPONENTS Development REQUIRED)
include(${ROOT_USE_FILE})

set(SNDSW_ROOT $ENV{SNDSW_ROOT})
set(ROOTSYS $ENV{ROOTSYS})
set(FAIRROOT_ROOT $ENV{FAIRROOT_ROOT})
set(FAIRLOGGER_ROOT $ENV{FAIRLOGGER_ROOT})
set(FMT_ROOT $ENV{FMT_ROOT})
set(VMC_ROOT $ENV{VMC_ROOT})
set(FAIRSHIP $ENV{FAIRSHIP})

set(ENV_ROOT_INCLUDE_PATH $ENV{ROOT_INCLUDE_PATH})
separate_arguments(ENV_ROOT_INCLUDE_PATH NATIVE_COMMAND "${ENV_ROOT_INCLUDE_PATH}")
if(UNIX)
    string(REPLACE ":" ";" ENV_ROOT_INCLUDE_PATH "${ENV_ROOT_INCLUDE_PATH}")
endif()


add_subdirectory(muonFluxUtils)
add_subdirectory(nTracks)
add_subdirectory(trkeff)
