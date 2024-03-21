cmake_minimum_required(VERSION 3.0)
set(GITHOME_MSU $ENV{GITHOME_MSU})
set(GSL_INCLUDE_DIR $ENV{GSL_INCLUDE_DIR})
project(hbt_project CXX)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
include_directories(SYSTEM include include/msu_coral include/msu_commonutils ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}")

file(GLOB msu_hbt_src ./src/hbt/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(msu_hbt STATIC ${msu_hbt_src})
target_link_libraries(msu_hbt ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})

file(GLOB msu_commonutils_src ./src/commonutils/*/*.cc ./src/commonutils/*/*/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(msu_commonutils STATIC ${msu_commonutils_src})
target_link_libraries(msu_commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})

file(GLOB msu_coral_src ./src/coral/*/*.cc ./src/coral/*/*/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(msu_coral STATIC ${msu_coral_src})
target_link_libraries(msu_coral ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})
