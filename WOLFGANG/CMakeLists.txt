cmake_minimum_required(VERSION 3.18)
project(wolfgang LANGUAGES CXX)

# Build DO_LIB first (IPACK depends on it)
add_subdirectory(DO_LIB)
add_subdirectory(IPACK)
