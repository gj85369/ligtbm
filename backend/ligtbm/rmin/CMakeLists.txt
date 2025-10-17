cmake_minimum_required(VERSION 3.0)

project(rmin
  LANGUAGES C)

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always HEAD
    COMMAND sed -e "s/-/./g"
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(GIT_VERSION "unknown")
endif()

# Options
option( PROFILE				"Generate profiling information"		OFF )


# Build flags
add_definitions(
  -D _GIT_VERSION_="\\"${GIT_VERSION}\\""
)

set(CMAKE_C_FLAGS
  "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_C_FLAGS_DEBUG
  "${CMAKE_C_FLAGS_DEBUG} -Wshadow -Wpointer-arith -Wcast-qual -Winline")

if(PROFILE)
  set(CMAKE_C_FLAGS "-pg ${CMAKE_C_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
endif()


# Locate libraries
find_package(mol2 REQUIRED)
find_package(sampling REQUIRED)

# Define executable target
set(SOURCES
  src/main.c)

add_executable(rmin ${SOURCES})
target_link_libraries(rmin
  mol2
  sampling
  m)


install(TARGETS rmin
  RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
