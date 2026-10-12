module fft
  use precision_mod
  implicit none
  ! private
  real(dp), parameter ::pi = 4.0_dp * atan(1.0_dp)
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)

  ! public :: say_hello
contains
  function read_data(filename) result(x)
    character(len=*), intent(in) :: filename
    complex(dp), allocatable :: x(:)
    integer :: n, sz
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

    sz = shiftl(1, BIT_SIZE(n) - leadz(n-1))

    allocate(x(sz), source = cmplx(0.0_dp, 0.0_dp, kind=dp))

    do it = 1, n
      read(iunit, *) re, im
      x(it) = cmplx(re, im, kind=dp)
    end do

    close(iunit)
  end function read_data

  subroutine write_result(x)
    complex(dp), intent(in) :: x(:)
    integer :: n
    integer :: it
    integer :: ounit, iostatus

    open(newunit=ounit,file="result.dat", &
     action='write', iostat=iostatus)

    if (iostatus /= 0) then
      error stop 'Error occured while opening file'
    end if
    
    n = size(x)

    do it = 1, n
      write(ounit, '(*(e23.15, 1x))') x(it)%re, x(it)%im
    end do
    close(ounit)

    open(newunit=ounit,file="abs.dat", &
     action='write', iostat=iostatus)

    if (iostatus /= 0) then
      error stop 'Error occured while opening file'
    end if
    do it = 1, n
      write(ounit, '(*(e23.15, 1x))') abs(x(it))
    end do
    close(ounit)

  end subroutine write_result

  function make_twiddles(n, inverse) result(w)
    integer, intent(in) :: n
    logical, intent(in) :: inverse
    complex(dp), allocatable :: w(:)
    real(dp) :: sgn
    integer :: k

    if (inverse) then
      sgn = 1.0_dp
    else 
      sgn = -1.0_dp
    endif

    allocate(w(n/2))
    
    do k = 1, n/2
      w(k) = exp(sgn * 2.0_dp * pi * i * real(k,dp)/real(n,dp))
    end do

  end function make_twiddles
end module fft
