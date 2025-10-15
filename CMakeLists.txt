cmake_minimum_required(VERSION 3.10)

# Project name and language
project(DipoleCalculation LANGUAGES CXX)

find_package(ROOT REQUIRED COMPONENTS Core Hist MathCore RIO Tree Physics Graf Graf3d Gpad Postscript)
include(${ROOT_USE_FILE})

if(NOT DEFINED ROOT_CXX_STANDARD)
  string(REPLACE " " ";" ROOT_CXX_FLAGS_LIST ${ROOT_CXX_FLAGS})
  foreach(fl ${ROOT_CXX_FLAGS_LIST})
    if ("${fl}" MATCHES ".*\\+\\+11$" OR "${fl}" MATCHES ".*\\+\\+0x$")
      set(ROOT_CXX_STANDARD 11)
    endif()
    if ("${fl}" MATCHES ".*\\+\\+14$" OR "${fl}" MATCHES ".*\\+\\+1y$")
      set(ROOT_CXX_STANDARD 14)
    endif()
    if ("${fl}" MATCHES ".*\\+\\+17$" OR "${fl}" MATCHES ".*\\+\\+1z$")
      set(ROOT_CXX_STANDARD 17)
    endif()
    if ("${fl}" MATCHES ".*\\+\\+20$" OR "${fl}" MATCHES ".*\\+\\+2a$")
      set(ROOT_CXX_STANDARD 20)
    endif()
    if ("${fl}" MATCHES ".*\\+\\+23$" OR "${fl}" MATCHES ".*\\+\\+2b$")
      set(ROOT_CXX_STANDARD 23)
    endif()
    if ("${fl}" MATCHES ".*\\+\\+26$" OR "${fl}" MATCHES ".*\\+\\+2c$")
      set(ROOT_CXX_STANDARD 26)
    endif()
  endforeach(fl ${ROOT_CXX_FLAGS_LIST})
endif()

if(DEFINED ROOT_CXX_STANDARD)
  if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
  elseif(NOT "${ROOT_CXX_STANDARD}" EQUAL "${CMAKE_CXX_STANDARD}")
    message(WARNING "ROOT_CXX_STANDARD=${ROOT_CXX_STANDARD} != CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}")
  endif()
else()
  message(WARNING "ROOT_CXX_STANDARD is undefined.")
endif()

# Ensure C++ standard is defined
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(WITH_HEPMC3
  "Path to HepMC3 installation (root of include/ and lib/)"
  "")
if (WITH_HEPMC3)
  list(APPEND CMAKE_PREFIX_PATH "${WITH_HEPMC3}")
endif()
find_package(HepMC3)
if (HepMC3_FOUND)
  message(STATUS "Found HepMC3 xxx: ${HEPMC3_INCLUDE_DIR} : ${HEPMC3_LIB}")
  if(NOT DEFINED HEPMC3_CXX_STANDARD)
    message(WARNING "HEPMC3_CXX_STANDARD is undefined.")
  elseif(NOT "${HEPMC3_CXX_STANDARD}" EQUAL "${CMAKE_CXX_STANDARD}")
    message(WARNING "HEPMC3_CXX_STANDARD=${HEPMC3_CXX_STANDARD} != CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}")
  endif()
  include_directories(${HEPMC3_INCLUDE_DIR})
  set(HEPMC_LIBRARIES ${HEPMC3_LIB})
else()
  message(FATAL_ERROR "HepMC3 not found. Define WITH_HEPMC3 manually.")
endif()

# Add optimization flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Source files
include_directories(${CMAKE_SOURCE_DIR}/include)
file(GLOB PROJECT_HEADERS "${CMAKE_SOURCE_DIR}/include/*.h")
file(GLOB PROJECT_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cxx")
set(MAIN_SOURCE "${CMAKE_SOURCE_DIR}/main.cxx")


# Define executable
add_executable(Run ${MAIN_SOURCE} ${PROJECT_SOURCES} ${PROJECT_HEADERS})

# Include headers from the current directory
target_include_directories(Run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# Link against ROOT libraries
target_link_libraries(Run
     PRIVATE
         HepMC3::HepMC3
         ROOT::Core
         ROOT::Hist
         ROOT::MathCore
         ROOT::RIO
         ROOT::Tree
         ROOT::Physics
         ROOT::Graf
         ROOT::Graf3d
         ROOT::Gpad
         ROOT::Postscript
)

# Glob all files starting with "hepmc" in your source directory
file(GLOB HEPMC_FILES "${CMAKE_SOURCE_DIR}/hepmc*")

# Copy to build directory (CMAKE_BINARY_DIR is the build directory)
file(COPY ${HEPMC_FILES} DESTINATION ${CMAKE_BINARY_DIR})
