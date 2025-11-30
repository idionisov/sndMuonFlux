cmake_minimum_required(VERSION 3.23)
project(muonFlux)

find_package(ROOT REQUIRED COMPONENTS RIO Core Hist MathCore Tree ROOTTPython)
include(${ROOT_USE_FILE})

set(SNDSW_ROOT $ENV{SNDSW_ROOT})
set(ROOTSYS $ENV{ROOTSYS})
set(FAIRROOT_ROOT $ENV{FAIRROOT_ROOT})
set(FAIRLOGGER_ROOT $ENV{FAIRLOGGER_ROOT})
set(FMT_ROOT $ENV{FMT_ROOT})
set(VMC_ROOT $ENV{VMC_ROOT})
set(FAIRSHIP $ENV{FAIRSHIP})



add_subdirectory(muonFluxUtils)
add_subdirectory(nTracks)
add_subdirectory(trkeff)
