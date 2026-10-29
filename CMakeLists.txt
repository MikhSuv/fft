cmake_minimum_required(VERSION 3.20)

project(fft
        DESCRIPTION "Fast Fourier Transform"
        LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

find_package(OpenMP REQUIRED COMPONENTS Fortran)

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
            src/iofft_mod.f90
            src/fft.f90
)

target_link_libraries(fast_Fourier PUBLIC OpenMP::OpenMP_Fortran)

add_executable(fft
              app/main.f90
)
target_link_libraries(fft PRIVATE fast_Fourier)

foreach(t IN ITEMS test_make_twiddles test_fft_dft test_fft_identity test_bit_reverse)
  add_executable(${t}
              test/${t}.f90
  )
  target_link_libraries(${t} PRIVATE fast_Fourier)
endforeach()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # PUBLIC so that every executable linking the library also links the
  # sanitizer runtime it was compiled with
  target_link_options(fast_Fourier PUBLIC
        -fsanitize=address,undefined
    )
endif()

enable_testing()
foreach(t IN ITEMS test_make_twiddles test_fft_dft test_fft_identity test_bit_reverse)
  add_test(NAME ${t} COMMAND ${t})
endforeach()
