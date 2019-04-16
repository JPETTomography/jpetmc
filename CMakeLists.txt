# CMAKE SETUPS
#######################################################
cmake_minimum_required (VERSION 2.6)
project(jpet_mc)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
file(GLOB sources   ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers   ${PROJECT_SOURCE_DIR}/include/*.h)
set(CMAKE_LIBRARY_PATH ${PROJECT_SOURCE_DIR}/lib)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wunused-parameter ")

# THE PROJECT VERSION NUMBER
#######################################################
set(VERSION_MAJOR   1   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)


# Provide scripts in a folder with a binary 
#######################################################
set(SCRIPT_FILES
    vis.mac
    run3.mac
    run5.mac
    run6.mac
    nema.mac
    init_vis.mac
    singleBeam.mac
    extendedSource.mac
    )


# Add a target to generate API documentation with Doxygen
#######################################################
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc)
endif()


# Include  GEANT4 WITH VIS DRIVERS
######################################################
# remember to compile geant4 with following flags -DGEANT4_USE_QT=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_GDML=ON -DGEANT4_BUILD_MULTITHREADED=ON
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})


# Include CADMesh 
######################################################
find_package(cadmesh REQUIRED)
include_directories(SYSTEM ${CADMESH_INCLUDE_DIRS})


# Include ROOT
######################################################
find_package(ROOT REQUIRED)
message(STATUS "ROOT version: ${ROOT_VERSION}")
if(ROOT_FOUND AND ROOT_VERSION VERSION_LESS "6.0")
    message(FATAL_ERROR "ROOT 6.0 is required ")
    # tested version above 6.10 can not display graphics while /run/beamOn due to llvm problem 
endif()

set(ROOT_CXX_FLAGS "-std=c++11 -Wunused-parameter")
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include) 
include_directories(SYSTEM  ${Geant4_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})


# Create ROOT dictionaries 
######################################################
ROOT_GENERATE_DICTIONARY(JPetMCClasses   include/JPetGeantDecayTree.h include/JPetGeantScinHits.h include/JPetGeantEventPack.h include/JPetGeantEventInformation.h LINKDEF LinkDef.h)
#---Create a shared library with geneated dictionary
add_library(JPetMCClassesDict SHARED  src/JPetGeantDecayTree.cpp src/JPetGeantScinHits.cpp  src/JPetGeantEventPack.cpp  src/JPetGeantEventInformation.cpp  JPetMCClasses.cxx )
target_link_libraries(JPetMCClassesDict ${ROOT_LIBRARIES})


# Create  a main program using the library
######################################################
add_executable(jpet_mc ${sources} ${headers})
target_link_libraries(jpet_mc ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} JPetMCClassesDict ${cadmesh_LIBRARIES})

foreach(file_i ${SCRIPT_FILES})
     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/${file_i} ${CMAKE_CURRENT_BINARY_DIR}/bin/${file_i} COPYONLY)
endforeach( file_i )

add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ln -sf  ${CMAKE_CURRENT_SOURCE_DIR}/stl_geometry ${CMAKE_CURRENT_BINARY_DIR}/bin
)
