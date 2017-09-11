      subroutine time_elapsed(s)
      implicit none
! This routine gets the elapsed CPU time in seconds
      real(8) :: s
      call cpu_time(s)
      end subroutine time_elapsed
