project(MPC)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11 -O3)

# set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources src/kinematic_model.cpp src/MPC.cpp src/main.cpp)

set(UWS_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/uWebSockets)

# compile uWebSockets
execute_process(COMMAND make WORKING_DIRECTORY ${UWS_LIB_PATH})

set(BUILD_INCLUDE_PATH ${CMAKE_BINARY_DIR}/include)
set(UWS_INCLUDE_PATH ${BUILD_INCLUDE_PATH}/uWS)
file(MAKE_DIRECTORY ${UWS_INCLUDE_PATH})

# Install uWebSockets header files
file (COPY
  ${UWS_LIB_PATH}/src/Extensions.h
  ${UWS_LIB_PATH}/src/WebSocketProtocol.h
  ${UWS_LIB_PATH}/src/Networking.h
  ${UWS_LIB_PATH}/src/WebSocket.h
  ${UWS_LIB_PATH}/src/Hub.h
  ${UWS_LIB_PATH}/src/Group.h
  ${UWS_LIB_PATH}/src/Node.h
  ${UWS_LIB_PATH}/src/Socket.h
  ${UWS_LIB_PATH}/src/HTTPSocket.h
  ${UWS_LIB_PATH}/src/uWS.h
  ${UWS_LIB_PATH}/src/uUV.h
  DESTINATION ${UWS_INCLUDE_PATH}
  )

include_directories(/usr/local/include)
link_directories(/usr/local/lib)
include_directories(src/Eigen-3.3)
include_directories(${BUILD_INCLUDE_PATH})

find_package(OpenSSL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_library(UV_LIBRARY NAMES uv)
find_library(UWS_LIBRARY NAMES uWS PATHS ${UWS_LIB_PATH} NO_DEFAULT_PATH)
find_library(IPOPT_LIBRARY NAMES ipopt)
find_library(M_LIBRARY NAMES m)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1*/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(mpc ${sources})
target_link_libraries(mpc ${IPOPT_LIBRARY} ${UV_LIBRARY} ${UWS_LIBRARY} ${ZLIB_LIBRARIES} ${OPENSSL_LIBRARIES})
