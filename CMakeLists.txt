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
#Tests
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG 08092139210881f01a34e60c5f2fd6afa0ea4024 #v3.16.1
)
FetchContent_MakeAvailable(Catch2)

add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

enable_testing()
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(tests)
# GPU
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
	enable_language(CUDA)
	set(CMAKE_CUDA_STANDARD 20)
	set(CMAKE_CUDA_STANDARD_REQUIRED ON)

	add_executable(lyapgpu gpu.cu)
	# RTX 3050
	set_target_properties(lyapgpu PROPERTIES CUDA_ARCHITECTURES 86)
	target_include_directories(lyapgpu PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
	target_compile_options(lyapgpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Werror=all-warnings>)
	target_include_directories(lyapgpu SYSTEM PRIVATE ${stb_SOURCE_DIR})
endif()
