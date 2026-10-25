module fft
  use precision_mod
  implicit none
  ! private
  real(dp), parameter ::pi = 4.0_dp * atan(1.0_dp)
  ! Imaginary unit
  complex(dp), parameter :: im_u = (0.0_dp, 1.0_dp)

  ! public :: say_hello
contains
  ! Reads a complex signal from a text file.
  !
  ! File layout:
  !   line 1      : "# <n>" -- header, the count of samples follows the '#'
  !   lines 2..n+1: "<re> <im>" -- one sample per line
  !
  ! The result is padded up to the next power of two with zeros
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
    ! Skip the leading '#' of the header and parse the sample count
    read(line(2:), *) n

    ! Round n up to the next power of two:
    ! BIT_SIZE(n) - leadz(n-1) is the index of the highest set bit of n-1,
    ! so 1 shifted there is the smallest power of two >= n
    sz = shiftl(1, BIT_SIZE(n) - leadz(n-1))

    ! Allocate a zero-padded buffer of power-of-two length
    allocate(x(sz), source = cmplx(0.0_dp, 0.0_dp, kind=dp))

    do it = 1, n
      read(iunit, *) re, im
      x(it) = cmplx(re, im, kind=dp)
    end do

    close(iunit)
  end function read_data

  ! Writes a complex array to two plain-text files in the current directory:
  !   result.dat -- "<re> <im>" per line
  !   abs.dat    -- "|x|" per line
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

  ! Builds the table of twiddle factors w(k) = exp(sgn * 2*pi*i * k/n)
  ! for k = 1, 2, ..., n/2.
  !
  !   inverse = .false. -> sgn = -1: forward transform
  !   inverse = .true.  -> sgn = +1: inverse transform
  function make_twiddles(n, inverse) result(w)
    integer, intent(in) :: n
    logical, intent(in) :: inverse
    complex(dp), allocatable :: w(:)
    real(dp) :: sgn
    integer :: k

    ! Sign of the exponent: +1 for the inverse transform, -1 for the forward one
    if (inverse) then
      sgn = 1.0_dp
    else 
      sgn = -1.0_dp
    endif

    allocate(w(n/2))
    
    do k = 1, n/2
      ! k/n cast to real to avoid integer division
      w(k) = exp(sgn * 2.0_dp * pi * im_u * real(k,dp)/real(n,dp))
    end do

  end function make_twiddles


end module fft
