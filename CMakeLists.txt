cmake_minimum_required(VERSION 3.20)

project(fft
        DESCRIPTION "Fast Fourier Transform"
        LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

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
    -std=f2023 \
    -DDEBUG")


set(CMAKE_Fortran_FLAGS_RELEASE "\
    -O3 \
    -march=native \
    -funroll-loops \
    -ftree-vectorize \
    -ffast-math \
    -fno-trapping-math \
    -std=f2023 \
    -DNDEBUG")

add_library(fast_Fourier STATIC
            src/precision_mod.f90
            src/fft.f90
)

add_executable(fft
              app/main.f90
)
target_link_libraries(fft PRIVATE fast_Fourier)

add_executable(test_make_twiddles
              test/test_make_twiddles.f90
)
target_link_libraries(test_make_twiddles PRIVATE fast_Fourier)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # PUBLIC so that every executable linking the library also links the
  # sanitizer runtime it was compiled with
  target_link_options(fast_Fourier PUBLIC
        -fsanitize=address,undefined
    )
endif()

enable_testing()
add_test(NAME make_twiddles COMMAND test_make_twiddles)
