cmake_minimum_required(VERSION 3.2)

project(sketch)

option(BUILDCXX "build pure c++ library instead of python bind interface" OFF)
#SET(CMAKE_CXX_COMPILER /usr/bin/g++)
#SET(CMAKE_C_COMPILER /usr/bin/gcc)

INCLUDE_DIRECTORIES(src)
#INCLUDE_DIRECTORIES(pybind11)

AUX_SOURCE_DIRECTORY(src DIR_SRCS)

#linux only 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ggdb -std=c++14 -march=native -fpermissive ")



if (BUILDCXX)
	add_compile_options(-DNOPYTHON)
	ADD_LIBRARY("${PROJECT_NAME}_dynamic" SHARED ${DIR_SRCS})
	ADD_LIBRARY("${PROJECT_NAME}_static" STATIC ${DIR_SRCS})
	
	set_target_properties("${PROJECT_NAME}_dynamic" PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
	set_target_properties("${PROJECT_NAME}_static" PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
else ()
	add_subdirectory(pybind11)
	pybind11_add_module(sketch ${DIR_SRCS})
endif ()
