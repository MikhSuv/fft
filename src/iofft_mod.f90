! Text-file I/O in the format:
!
!   line 1        : "# <n>" -- sample count after '#' and one or more spaces
!   lines 2..n+1  : "<re> <im>" -- real and imaginary parts per sample
!
module iofft_mod
  use precision_mod, only: dp
  use iso_fortran_env, only: error_unit
  implicit none
  private
  public :: read_data, write_result

contains

  ! Reports a fatal I/O problem on stderr and stops with a non-zero status
  subroutine report_error(message)
    character(len=*), intent(in) :: message

    write(error_unit, '(a)') 'io_mod: ' // message
    error stop 1
  end subroutine report_error

  ! Reads a complex signal from a plain-text file.
  ! The buffer is zero padded up to the next power of two.
  function read_data(filename) result(x)
    character(len=*), intent(in) :: filename
    complex(dp), allocatable :: x(:)
    integer :: n, sz
    real(dp) :: re, im
    integer :: iostat, iunit, it
    character(len=256) :: line

    open(newunit=iunit, file=filename, status='old', action='read', iostat=iostat)
    if (iostat /= 0) call report_error('cannot open <' // filename // '> for reading')

    read(iunit, '(a)', iostat=iostat) line
    if (iostat /= 0) call report_error('cannot read the header line of <' // filename // '>')
    if (line(1:1) /= '#') call report_error('the header line of <' // filename // '> must start with #')
    read(line(2:), *, iostat=iostat) n
    if (iostat /= 0 .or. n < 1) call report_error('bad sample count in the header line of <' // filename // '>')

    ! Smallest power of two sz >= n: bit_size(n) - leadz(n-1) is the index
    ! of the highest set bit of n-1, so 1 shifted there goes just past it
    sz = shiftl(1, bit_size(n) - leadz(n - 1))

    allocate(x(sz), source=cmplx(0.0_dp, 0.0_dp, kind=dp))

    do it = 1, n
      read(iunit, *, iostat=iostat) re, im
      x(it) = cmplx(re, im, kind=dp)
    end do

    close(iunit)
  end function read_data

  ! Writes a complex array to two plain-text files:
  !   result_file : "# <n>" header + one "<re> <im>" pair per line (like the input)
  !   abs_file    : one |x(k)| per line, without the header
  subroutine write_result(x, result_file, abs_file)
    complex(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: result_file
    character(len=*), intent(in) :: abs_file
    integer :: n, iunit, iostat, it

    n = size(x)

    open(newunit=iunit, file=result_file, status='replace', action='write', iostat=iostat)
    if (iostat /= 0) call report_error('cannot open <' // result_file // '> for writing')
    write(iunit, '("# ", i0)', iostat=iostat) n
    if (iostat /= 0) call report_error('cannot write the header of <' // result_file // '>')
    do it = 1, n
      write(iunit, '(*(e23.15, 1x))', iostat=iostat) x(it)%re, x(it)%im
    end do
    close(iunit)

    open(newunit=iunit, file=abs_file, status='replace', action='write', iostat=iostat)
    if (iostat /= 0) call report_error('cannot open <' // abs_file // '> for writing')
    do it = 1, n
      write(iunit, '(*(e23.15, 1x))', iostat=iostat) abs(x(it))
    end do
    close(iunit)
  end subroutine write_result

end module iofft_mod 
