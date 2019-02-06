# Install script for directory: /Users/sodtaj/surface_projects/MC_bd/qhull-2011.1

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libqhullcpp.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib" TYPE STATIC_LIBRARY FILES "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/libqhullcpp.a")
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/libqhullcpp.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libqhullcpp.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}/usr/local/lib/libqhullcpp.a")
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libqhull6.6.2.0.1385.dylib;/usr/local/lib/libqhull6.dylib")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib" TYPE SHARED_LIBRARY FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/libqhull6.6.2.0.1385.dylib"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/libqhull6.dylib"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/lib/libqhull6.6.2.0.1385.dylib"
      "$ENV{DESTDIR}/usr/local/lib/libqhull6.dylib"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      execute_process(COMMAND "/usr/bin/install_name_tool"
        -id "/usr/local/lib/libqhull6.6.2.0.1385.dylib"
        "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libqhullstatic.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib" TYPE STATIC_LIBRARY FILES "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/libqhullstatic.a")
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/libqhullstatic.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libqhullstatic.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}/usr/local/lib/libqhullstatic.a")
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libqhullstatic_p.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib" TYPE STATIC_LIBRARY FILES "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/libqhullstatic_p.a")
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/libqhullstatic_p.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libqhullstatic_p.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}/usr/local/lib/libqhullstatic_p.a")
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/qhull-6.2.0.1385;/usr/local/bin/qhull")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qhull-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qhull"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/qhull-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/qhull"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/rbox-6.2.0.1385;/usr/local/bin/rbox")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/rbox-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/rbox"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/rbox-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/rbox"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/qconvex-6.2.0.1385;/usr/local/bin/qconvex")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qconvex-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qconvex"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/qconvex-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/qconvex"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/qdelaunay-6.2.0.1385;/usr/local/bin/qdelaunay")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qdelaunay-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qdelaunay"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/qdelaunay-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/qdelaunay"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/qvoronoi-6.2.0.1385;/usr/local/bin/qvoronoi")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qvoronoi-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qvoronoi"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/qvoronoi-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/qvoronoi"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/qhalf-6.2.0.1385;/usr/local/bin/qhalf")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qhalf-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/qhalf"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/qhalf-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/qhalf"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/user_eg-6.2.0.1385;/usr/local/bin/user_eg")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/user_eg-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/user_eg"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/user_eg-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/user_eg"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      execute_process(COMMAND "/usr/bin/install_name_tool"
        -change "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/libqhull6.6.2.0.1385.dylib" "/usr/local/lib/libqhull6.6.2.0.1385.dylib"
        "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/user_eg2-6.2.0.1385;/usr/local/bin/user_eg2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/user_eg2-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/user_eg2"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/user_eg2-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/user_eg2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/user_eg3-6.2.0.1385;/usr/local/bin/user_eg3")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/user_eg3-6.2.0.1385"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/user_eg3"
    )
  foreach(file
      "$ENV{DESTDIR}/usr/local/bin/user_eg3-6.2.0.1385"
      "$ENV{DESTDIR}/usr/local/bin/user_eg3"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/libqhull/libqhull.h;/usr/local/include/libqhull/geom.h;/usr/local/include/libqhull/io.h;/usr/local/include/libqhull/mem.h;/usr/local/include/libqhull/merge.h;/usr/local/include/libqhull/poly.h;/usr/local/include/libqhull/qhull_a.h;/usr/local/include/libqhull/qset.h;/usr/local/include/libqhull/random.h;/usr/local/include/libqhull/stat.h;/usr/local/include/libqhull/user.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/include/libqhull" TYPE FILE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/libqhull.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/geom.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/io.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/mem.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/merge.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/poly.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/qhull_a.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/qset.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/random.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/stat.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhull/user.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/libqhullcpp/Coordinates.h;/usr/local/include/libqhullcpp/functionObjects.h;/usr/local/include/libqhullcpp/PointCoordinates.h;/usr/local/include/libqhullcpp/Qhull.h;/usr/local/include/libqhullcpp/QhullError.h;/usr/local/include/libqhullcpp/QhullFacet.h;/usr/local/include/libqhullcpp/QhullFacetList.h;/usr/local/include/libqhullcpp/QhullFacetSet.h;/usr/local/include/libqhullcpp/QhullHyperplane.h;/usr/local/include/libqhullcpp/QhullIterator.h;/usr/local/include/libqhullcpp/QhullLinkedList.h;/usr/local/include/libqhullcpp/QhullPoint.h;/usr/local/include/libqhullcpp/QhullPoints.h;/usr/local/include/libqhullcpp/QhullPointSet.h;/usr/local/include/libqhullcpp/QhullQh.h;/usr/local/include/libqhullcpp/QhullRidge.h;/usr/local/include/libqhullcpp/QhullSet.h;/usr/local/include/libqhullcpp/QhullSets.h;/usr/local/include/libqhullcpp/QhullStat.h;/usr/local/include/libqhullcpp/QhullVertex.h;/usr/local/include/libqhullcpp/QhullVertexSet.h;/usr/local/include/libqhullcpp/RboxPoints.h;/usr/local/include/libqhullcpp/UsingLibQhull.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/include/libqhullcpp" TYPE FILE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/Coordinates.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/functionObjects.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/PointCoordinates.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/Qhull.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullError.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullFacet.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullFacetList.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullFacetSet.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullHyperplane.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullIterator.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullLinkedList.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullPoint.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullPoints.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullPointSet.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullQh.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullRidge.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullSet.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullSets.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullStat.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullVertex.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/QhullVertexSet.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/RboxPoints.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/libqhullcpp/UsingLibQhull.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/road/RoadError.h;/usr/local/include/road/RoadLogEvent.h;/usr/local/include/road/RoadTest.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/include/road" TYPE FILE FILES
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/road/RoadError.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/road/RoadLogEvent.h"
    "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/src/road/RoadTest.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/man/man1/qhull.1")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/man/man1" TYPE FILE RENAME "qhull.1" FILES "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/html/qhull.man")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/man/man1/rbox.1")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/man/man1" TYPE FILE RENAME "rbox.1" FILES "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/html/rbox.man")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/doc/packages/qhull/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/share/doc/packages/qhull" TYPE DIRECTORY FILES "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/html/")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/sodtaj/surface_projects/MC_bd/qhull-2011.1/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
