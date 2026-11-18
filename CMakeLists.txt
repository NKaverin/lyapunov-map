cmake_minimum_required(VERSION 3.24)
project(lyapunov-map LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(lyapmap main.cpp stb_image_write.cpp)

include(FetchContent)
FetchContent_Declare(
	stb
	GIT_REPOSITORY https://github.com/nothings/stb.git
	GIT_TAG 4e9b368bfc071e97206ea5c6adba68dd0cbf8897
)
FetchContent_MakeAvailable(stb)

target_include_directories(lyapmap SYSTEM PRIVATE ${stb_SOURCE_DIR})
find_package(OpenMP REQUIRED)
target_link_libraries(lyapmap PRIVATE OpenMP::OpenMP_CXX)

if (MSVC)
	target_compile_options(lyapmap PRIVATE /W4)
else()
	target_compile_options(lyapmap PRIVATE -Wall -Wextra)
endif()
