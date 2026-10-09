module fft
  use precision_mod
  implicit none
  ! private

  ! public :: say_hello
contains
  function read_data(filename) result(x)
    character(len=*), intent(in) :: filename
    complex(dp), allocatable :: x(:)
    integer :: n, size
    real(dp) :: re, im

    integer :: iostatus, iunit, it
    character(len=256) :: line
    open(newunit=iunit, file=filename, status='old', &
    action = 'read', iostat=iostatus)

    if (iostatus /= 0) then
      error stop 'Error occured while opening file'
    end if
    read(iunit, '(a)', iostat=iostatus) line
    if (iostatus /= 0) then
      error stop 'Error occured while reading line'
    end if
    read(line(2:), *) n

    size = shiftl(1, BIT_SIZE(n) - leadz(n-1))

    allocate(x(size), source = complex(0.0_dp, 0.0_dp))

    do it = 1, n
      read(iunit, *) re, im
      x(it) = complex(re, im)
    end do

    close(iunit)
  end function read_data
end module fft
