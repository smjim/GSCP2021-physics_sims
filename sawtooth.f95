program sawtooth
      !calculates values of sawtooth wave
      implicit none

      real(8) :: x, y, pi
      x = 0.0

      open(unit=100, file='sawtooth.dat')

      pi = 4.0*atan(1.0)
      write (100,*) '#x-values                  y-values'

      do while (x<= 2.0*pi)
        if (x <= pi) then
                y = x
        elseif (x > pi) then
                y = x - 2.0*pi
        else
                exit
        endif
        !y = sin(x)
        !x_in_degrees = x * 180.0/pi
        write(100,*) x, y
        x = x + 0.1
      end do

      close(100)

      end program sawtooth
       
