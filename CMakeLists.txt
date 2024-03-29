# Declare minimum required CMake version before anything to avoid problems.
SET(HW2Testing_TRIBITS_CMAKE_MINIMUM_REQUIRED 2.8.11)
CMAKE_MINIMUM_REQUIRED(VERSION ${HW2Testing_TRIBITS_CMAKE_MINIMUM_REQUIRED})

#
# Define project name and set up major options
#
INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)

#Indicate to TriBITS that we might use Python
SET(${PROJECT_NAME}_USES_PYTHON TRUE)

#Automatically read in external repositories
SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE "Continuous" CACHE STRING
  "Enabled in top level CMakeLists.txt")

#Generate repo versions text file
SET(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE ON CACHE BOOL
  "Enabled in top level CMakeLists.txt")

#CMake still requires the PROJECT be defined in the top level CMakeLists.txt file
PROJECT(${PROJECT_NAME} NONE)

#
#Pull in TriBITS for the project
#
SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits"
  CACHE PATH "" )
MESSAGE(STATUS "${PROJECT_NAME}_TRIBITS_DIR = '${${PROJECT_NAME}_TRIBITS_DIR}'")
INCLUDE(${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake)

#
#Enable Tests always
#
SET(${PROJECT_NAME}_ENABLE_HW2Tests ON CACHE BOOL "")

SET(TPL_ENABLE_MPI OFF CACHE BOOL "")

#
#Process the TriBITS Project
#
TRIBITS_PROJECT()
