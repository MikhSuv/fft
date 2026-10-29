! Reads a complex signal, computes the transform (forward by default) and
! writes the result in the assignment format:
!
!   fft [<data.dat>] [<result.dat>] [<abs.dat>] [inverse|-i]
!
program main
  use precision_mod, only: dp
  use iofft_mod, only: read_data, write_result
  use fft, only: fft_transform
  implicit none

  complex(dp), allocatable :: x(:)
  character(len=256) :: arg
  character(len=:), allocatable :: in_file, result_file, abs_file
  integer :: i, pos
  logical :: inverse

  in_file = 'data.dat'
  result_file = 'result.dat'
  abs_file = 'abs.dat'
  inverse = .false.
  pos = 0
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    if (len_trim(arg) == 0) cycle
    select case (trim(arg))
    case ('-i', '--inverse', 'inverse')
      inverse = .true.
    case default
      pos = pos + 1
      select case (pos)
      case (1)
        in_file = trim(arg)
      case (2)
        result_file = trim(arg)
      case (3)
        abs_file = trim(arg)
      case default
        error stop 'usage: fft [<data.dat>] [<result.dat>] [<abs.dat>] [inverse|-i]'
      end select
    end select
  end do

  x = read_data(in_file)

  call fft_transform(x, inverse)

  call write_result(x, result_file, abs_file)
  deallocate(x, in_file, result_file, abs_file)
end program main
