cmake_minimum_required(VERSION 3.15)
project(final_test)

set(CMAKE_CXX_STANDARD 14)

add_executable(final_test main.cpp)

find_package(ALGLIB REQUIRED)
target_link_libraries(${PROJECT_NAME}
        alglib)