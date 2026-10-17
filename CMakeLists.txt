cmake_minimum_required(VERSION 3.20)
project(options_pricing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_library(pricer
    src/BlackScholes.cpp
    src/Payoff.cpp
    src/PDESolver.cpp
    src/Tridiagonal.cpp
)
target_include_directories(pricer PUBLIC include)

add_executable(price apps/price.cpp)
target_link_libraries(price PRIVATE pricer)

add_executable(benchmark apps/benchmark.cpp)
target_link_libraries(benchmark PRIVATE pricer)

enable_testing()
add_executable(pricer_tests tests/test_pricer.cpp)
target_link_libraries(pricer_tests PRIVATE pricer)
add_test(NAME pricer_tests COMMAND pricer_tests)
