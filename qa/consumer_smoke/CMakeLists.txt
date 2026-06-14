cmake_minimum_required(VERSION 3.20)
project(qd3_consumer_smoke LANGUAGES CXX)

find_package(QD3 REQUIRED)
add_executable(app main.cpp)
target_link_libraries(app PRIVATE QD3::qd)
