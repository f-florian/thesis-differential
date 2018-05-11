# -*- mode: cmake; eval: (hl-line-mode 0);-*-
cmake_minimum_required(VERSION 3.0.0)
project(differential VERSION 1.0 LANGUAGES CXX)
set(sources "differential.cpp")

# add_definitions(-Dcmake_MAJOR=${PROJECT_VERSION_MAJOR} -DHAVE_INLINE -Dcmake_MINOR=${PROJECT_VERSION_MINOR} -Dcmake_NAME="${progname}" -D___DEBUG_PRINT="1")
add_definitions(-DHAVE_INLINE -Dcmake_MAJOR=${PROJECT_VERSION_MAJOR} -Dcmake_MINOR=${PROJECT_VERSION_MINOR} -Dcmake_NAME="${progname}")
add_compile_options(-pipe -Werror -Wall -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common --std=gnu++17 -O1)
set(progname ${PROJECT_NAME})

add_library(${progname} SHARED ${sources})
install(TARGETS ${progname} LIBRARY DESTINATION "lib" PERMISSIONS WORLD_READ)
install(FILES "differential.h" DESTINATION "include" PERMISSIONS WORLD_READ)
install(FILES "LICENSE" DESTINATION "var/lib/${progname}" PERMISSIONS WORLD_READ)
