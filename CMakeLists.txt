cmake_minimum_required(VERSION 3.16)

project(pathplanner LANGUAGES CXX)

set(CMAKE_AUTOUIC OFF)
set(CMAKE_AUTOMOC OFF)
set(CMAKE_AUTORCC OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core Gui)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core Gui)

set(EIGEN5_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/eigen-5.0.0)
include_directories(${EIGEN5_INCLUDE_DIR})

add_executable(pathplanner
  main.cpp
  pathplanner.cpp
  pathplanner.h
)
target_link_libraries(pathplanner
  PRIVATE Qt${QT_VERSION_MAJOR}::Core
)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

include(GNUInstallDirs)
install(TARGETS pathplanner
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)


