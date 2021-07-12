program fouriersawtooth
      !calculates values of sawtooth wave using fourier series approximation
      implicit none

      real(8) :: x, y, pi
      integer(4) :: n, iterations
      x = 0.0
      iterations = 20 

      open(unit=100, file='fouriersawtooth.dat')

      pi = 4.0*atan(1.0)
      write (100,*) '#x-values                  y-values'

      do while (x <= 4.0*pi)
        y = 0

        do n = 1, iterations 
                !summation of fourier series formula for sawtooth
                y = y + (-1.0)**(n+1.0)*(2.0/n) * (sin(n*x))
                !random function
                !y = y + sin((2*n+1)*x)
        end do

        write(100,*) x, y
        x = x + 0.01
      end do

      close(100)

      end program fouriersawtooth
