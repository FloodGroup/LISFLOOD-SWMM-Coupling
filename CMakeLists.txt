cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
# Thin wrapper: the full LISFLOOD-FP + SWMM build lives in src/
project(lisflood-swmm-coupling LANGUAGES NONE)
add_subdirectory(src)
