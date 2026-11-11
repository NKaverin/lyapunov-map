cmake_minimum_required(VERSION 3.24)
project(lyapunov-map LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(lyapmap main.cpp)

include(FetchContent)
FetchContent_Declare(
	stb
	GIT_REPOSITORY https://github.com/nothings/stb.git
	GIT_TAG 4e9b368bfc071e97206ea5c6adba68dd0cbf8897
)
FetchContent_MakeAvailable(stb)

target_include_directories(lyapmap PRIVATE ${stb_SOURCE_DIR})
