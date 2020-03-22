if(BUILD_VIEW)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif(BUILD_VIEW)

find_package(PkgConfig REQUIRED)
find_package(Boost 1.58 REQUIRED COMPONENTS program_options)
find_package(jsoncpp REQUIRED)

if(BUILD_VIEW)
    find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets OpenGL Svg)
    find_package(CGAL REQUIRED COMPONENTS CGAL Core Qt5)
endif(BUILD_VIEW)

set(SOURCES decomp.cpp file_io.cpp gmobj.cpp main.cpp mesh.cpp view.cpp
    ${PROJECT_BINARY_DIR}/gitversion.h
)
set(UIS mcd-ref-merge.ui)
set(RESOURCES mcd-ref-merge.qrc)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(mcd-ref-merge ${SOURCES} ${UIS} ${RESOURCES})

set_target_properties(mcd-ref-merge PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

target_link_libraries(mcd-ref-merge Boost::program_options)
target_link_libraries(mcd-ref-merge jsoncpp)
target_link_libraries(mcd-ref-merge triangle)

if(BUILD_VIEW)
    target_link_libraries(mcd-ref-merge CGAL::CGAL CGAL::CGAL_Core CGAL::CGAL_Qt5)
    target_link_libraries(mcd-ref-merge Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL Qt5::Svg)
endif(BUILD_VIEW)

add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/gitversion.h"
  DEPENDS "${CMAKE_SOURCE_DIR}/.git/logs/HEAD"
  COMMAND ${CMAKE_COMMAND}
    -Dlocal_dir="${CMAKE_SOURCE_DIR}"
    -Doutput_dir="${PROJECT_BINARY_DIR}"
    -P "${CMAKE_SOURCE_DIR}/cmake/gitversion.cmake"
)
