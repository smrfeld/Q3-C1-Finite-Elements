cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(q3c1)

include_directories(/usr/local/include)

AUX_SOURCE_DIRECTORY(src SRC_FILES) 

file (GLOB HEADER_FILES "include/q3c1_bits/*.hpp")
file (GLOB HEADER_FILES_FWDS "include/q3c1_bits/fwds/*.hpp")

add_library(q3c1 SHARED 
	src/basis_func.cpp
	src/cell.cpp
	src/dimension_1d.cpp
	src/general.cpp
	src/grid.cpp
	src/idx_set.cpp
	src/vertex.cpp
	include/q3c1_bits/basis_func.hpp
	include/q3c1_bits/cell.hpp
	include/q3c1_bits/dimension_1d.hpp
	include/q3c1_bits/general.hpp
	include/q3c1_bits/grid.hpp
	include/q3c1_bits/idx_set.hpp
	include/q3c1_bits/vertex.hpp
	${HEADER_FILES_FWDS}
)

target_include_directories(q3c1 PRIVATE include/q3c1_bits)

find_library(ARMADILLO_LIB armadillo HINTS /usr/local/lib/)
target_link_libraries(q3c1 PUBLIC ${ARMADILLO_LIB})

install(TARGETS q3c1 DESTINATION lib)

install(FILES 
	include/q3c1 
	DESTINATION include)

install(DIRECTORY 
	include/q3c1_bits 
	DESTINATION include)