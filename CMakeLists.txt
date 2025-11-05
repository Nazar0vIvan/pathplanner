cmake_minimum_required(VERSION 3.16)

project(linalg LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core Gui)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core Gui)

set(EIGEN5_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/eigen-5.0.0)
include_directories(${EIGEN5_INCLUDE_DIR})

add_executable(linalg
  main.cpp
  linalg.h
)
target_link_libraries(linalg
  PRIVATE Qt${QT_VERSION_MAJOR}::Core
  PRIVATE Qt${QT_VERSION_MAJOR}::Gui
)

include(GNUInstallDirs)
install(TARGETS linalg
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
