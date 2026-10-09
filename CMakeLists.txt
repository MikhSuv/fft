cmake_minimum_required(VERSION 3.20)

project(fft
        DESCRIPTION "Fast Fourier Transform"
        LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DERECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
  # set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "\
    -g \
    -O0 \
    -Wall \
    -Wextra \
    -fcheck=all \
    -fbacktrace \
    -fsanitize=address,undefined \
    -fno-omit-frame-pointer \
    -ffpe-trap=invalid,zero,overflow \
    -finit-real=snan \
    -fstack-protector-strong \
    -DDEBUG")


set(CMAKE_Fortran_FLAGS_RELEASE "\
    -O3 \
    -march=native \
    -funroll-loops \
    -ftree-vectorize \
    -ffast-math \
    -fno-trapping-math \
    -DNDEBUG")

add_library(fast_Fourier STATIC
            src/fft.f90
)

add_executable(fft
              app/main.f90
)
target_link_libraries(fft PRIVATE fast_Fourier)

add_executable(fft-test
              test/check.f90
)
target_link_libraries(fft-test PRIVATE fast_Fourier)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_link_options(${PROJECT_NAME} PRIVATE
        -fsanitize=address,undefined
    )
endif()
