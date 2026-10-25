program main
  use precision_mod
  use fft
  implicit none

  complex(dp), allocatable :: x(:)
  x = read_data("data.dat")
  print *, x
  call write_result(x)
  deallocate(x)
  x = make_twiddles(10, .false.)
  print *, x
  deallocate(x)
end program main
