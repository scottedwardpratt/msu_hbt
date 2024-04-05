cmake_minimum_required(VERSION 3.5)
set(CMAKE_VERBOSE_MAKEFILE on)
project(msu_hbt_project CXX)
set(GITHOME_MSU_HBT $ENV{GITHOME_MSU_HBT})
if( NOT DEFINED ENV{GITHOME_MSU_HBT} )
 set(GITHOME_MSU_HBT "$ENV{PWD}/../..")
endif()
message("GITHOME_MSU_HBT = ${GITHOME_MSU_HBT}")

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{PWD}/../bin")

add_custom_target(extern_msu_hbt COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_HBT}/software)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU_HBT}/software/include)

link_directories(${GITHOME_MSU_HBT}/software/lib)

find_library(${GITHOME_MSU_HBT}/software/lib msu_commonutils.a)
find_library(${GITHOME_MSU_HBT}/software/lib msu_coral.a)
find_library(${GITHOME_MSU_HBT}/software/lib msu_hbt.a)

#####################

set(SOURCE_FILES ./main_programs/hbt_fromGauss.cc)
add_executable(hbt_fromGauss ${SOURCE_FILES})
SET_SOURCE_FILES_PROPERTIES(hbt_fromGauss.cc PROPERTIES OBJECT_DEPENDS extern_msu_hbt)
add_dependencies(hbt_fromGauss extern_msu_hbt)

target_link_libraries(hbt_fromGauss msu_hbt msu_commonutils msu_coral ${GITHOME_MSU_HBT}/software/lib/libmsu_commonutils.a ${GITHOME_MSU_HBT}/software/lib/libmsu_coral.a
${GITHOME_MSU_HBT}/software/lib/libmsu_hbt.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#####################
#####################

set(SOURCE_FILES ./main_programs/hbt_fromOSCAR.cc)
add_executable(hbt_fromOSCAR ${SOURCE_FILES})
SET_SOURCE_FILES_PROPERTIES(hbt_fromOSCAR.cc PROPERTIES OBJECT_DEPENDS extern_msu_hbt)
add_dependencies(hbt_fromOSCAR extern_msu_hbt)

target_link_libraries(hbt_fromOSCAR msu_hbt msu_commonutils msu_coral ${GITHOME_MSU_HBT}/software/lib/libmsu_commonutils.a ${GITHOME_MSU_HBT}/software/lib/libmsu_coral.a
${GITHOME_MSU_HBT}/software/lib/libmsu_hbt.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#####################
#####################

set(SOURCE_FILES ./main_programs/hbt_CalcWF.cc)
add_executable(hbt_CalcWF ${SOURCE_FILES})
SET_SOURCE_FILES_PROPERTIES(hbt_CalcWF.cc PROPERTIES OBJECT_DEPENDS extern_msu_hbt)
add_dependencies(hbt_CalcWF extern_msu_hbt)

target_link_libraries(hbt_CalcWF msu_hbt msu_commonutils msu_coral ${GITHOME_MSU_HBT}/software/lib/libmsu_commonutils.a ${GITHOME_MSU_HBT}/software/lib/libmsu_coral.a
${GITHOME_MSU_HBT}/software/lib/libmsu_hbt.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#####################