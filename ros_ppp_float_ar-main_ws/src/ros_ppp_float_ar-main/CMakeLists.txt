cmake_minimum_required(VERSION 3.0.2)
project(rtkrcv)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  message_generation
  gnss_comm
)

## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)


## Uncomment this if the package has a setup.py. This macro ensures
## modules and global scripts declared therein get installed
## See http://ros.org/doc/api/catkin/html/user_guide/setup_dot_py.html
# catkin_python_setup()

################################################
## Declare ROS messages, services and actions ##
################################################

## To declare and build messages, services or actions from within this
## package, follow these steps:
## * Let MSG_DEP_SET be the set of packages whose message types you use in
##   your messages/services/actions (e.g. std_msgs, actionlib_msgs, ...).
## * In the file package.xml:
##   * add a build_depend tag for "message_generation"
##   * add a build_depend and a exec_depend tag for each package in MSG_DEP_SET
##   * If MSG_DEP_SET isn't empty the following dependency has been pulled in
##     but can be declared for certainty nonetheless:
##     * add a exec_depend tag for "message_runtime"
## * In this file (CMakeLists.txt):
##   * add "message_generation" and every package in MSG_DEP_SET to
##     find_package(catkin REQUIRED COMPONENTS ...)
##   * add "message_runtime" and every package in MSG_DEP_SET to
##     catkin_package(CATKIN_DEPENDS ...)
##   * uncomment the add_*_files sections below as needed
##     and list every .msg/.srv/.action file to be processed
##   * uncomment the generate_messages entry below
##   * add every package in MSG_DEP_SET to generate_messages(DEPENDENCIES ...)

## Generate messages in the 'msg' folder
add_message_files(
   FILES
   epoch_time_msg.msg
   ssr_msg.msg
   corr_msg.msg
)

## Generate services in the 'srv' folder
# add_service_files(
#   FILES
#   Service1.srv
#   Service2.srv
# )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
generate_messages(
   DEPENDENCIES
   std_msgs
)

################################################
## Declare ROS dynamic reconfigure parameters ##
################################################

## To declare and build dynamic reconfigure parameters within this
## package, follow these steps:
## * In the file package.xml:
##   * add a build_depend and a exec_depend tag for "dynamic_reconfigure"
## * In this file (CMakeLists.txt):
##   * add "dynamic_reconfigure" to
##     find_package(catkin REQUIRED COMPONENTS ...)
##   * uncomment the "generate_dynamic_reconfigure_options" section below
##     and list every .cfg file to be processed

## Generate dynamic reconfigure parameters in the 'cfg' folder
# generate_dynamic_reconfigure_options(
#   cfg/DynReconf1.cfg
#   cfg/DynReconf2.cfg
# )

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if your package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES rtkrcv
   CATKIN_DEPENDS roscpp std_msgs message_runtime
#  DEPENDS system_lib
)

###########
## Build ##
###########

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

add_compile_definitions(ENAGLO ENAGAL ENACMP NFREQ=5 SVR_REUSEADDR)

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
# include
  ${catkin_INCLUDE_DIRS}
  src/rtklib
)

## Declare a C++ library
add_library(rtklib SHARED
   src/rtklib/vt.c
   src/rtklib/rtkcmn.c
   src/rtklib/rtksvr.c
   src/rtklib/rtkpos.c
   src/rtklib/geoid.c
   src/rtklib/solution.c
   src/rtklib/lambda.c
   src/rtklib/sbas.c
   src/rtklib/stream.c
   src/rtklib/rcvraw.c
   src/rtklib/rtcm.c
   src/rtklib/rtcm2.c
   src/rtklib/rtcm3.c
   src/rtklib/rtcm3e.c
   src/rtklib/preceph.c
   src/rtklib/options.c
   src/rtklib/pntpos.c
   src/rtklib/ppp.c
   src/rtklib/ppp_ar.c
   src/rtklib/ppp_corr.c
   src/rtklib/ephemeris.c
   src/rtklib/qzslex.c
   src/rtklib/rinex.c
   src/rtklib/ionex.c
   src/rtklib/tides.c
   src/rtklib/rcv/novatel.c
   src/rtklib/rcv/ublox.c
   src/rtklib/rcv/ss2.c
   src/rtklib/rcv/crescent.c
   src/rtklib/rcv/skytraq.c
   src/rtklib/rcv/javad.c
   src/rtklib/rcv/nvs.c
   src/rtklib/rcv/binex.c
   src/rtklib/rcv/rt17.c
   src/rtklib/rcv/septentrio.c
   src/rtklib/rcv/swiftnav.c
   src/rtklib/rcv/gw10.c
   src/rtklib/rcv/cmr.c
   src/rtklib/rcv/comnav.c
   src/rtklib/rcv/tersus.c
)

add_executable(rtkrcv-node src/ros_rtkrcv.cpp)
target_link_libraries(rtkrcv-node ${catkin_LIBRARIES} rtklib m rt pthread)
add_dependencies(rtkrcv-node rtkrcv_generate_messages_cpp)


