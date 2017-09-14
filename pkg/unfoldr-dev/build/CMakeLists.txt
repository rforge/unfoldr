PROJECT( unfoldr )

cmake_minimum_required(VERSION 2.6)

SET(CMAKE_CXX_FLAGS "-fopenmp")
SET(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")

# Load compiler-specific information.
IF(CMAKE_CXX_COMPILER_ID)
  INCLUDE(Compiler/${CMAKE_CXX_COMPILER_ID}-CXX OPTIONAL)
ENDIF(CMAKE_CXX_COMPILER_ID)

# Some interface includes
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} /usr/include/ /usr/local/include/ /usr/local/lib/R/include/ /usr/share/R/include/)
LINK_DIRECTORIES( /usr/lib /usr/lib/R/lib /usr/local/lib )

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib)
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

set(thirdparty
   R
)

set( SRC_GL
 ../src/unfold.cpp
 ../src/ellipt.cpp
 ../src/GeometricPrimitives.cpp
 ../src/Vector.h 
 ../src/Intersector.cpp 
 ../src/SimSphere.cpp
 ../src/SimEllipsoid.cpp
 ../src/SimCylinder.cpp 
 ../src/directions.cpp
 ../src/Rheaders.h
 ../src/Utils.cpp 
)

add_library( unfoldr SHARED ${SRC_GL} )
#add_executable(libmain ../src/main.cpp)
#target_link_libraries(libmain STLib ${thirdparty} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libFBASpheroid.so)

