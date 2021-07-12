program eulermethod
        !calculates stuff with euler (chromer(?)) method
        implicit none

	real(8) :: y, v, a, t, tau, v_temp

        !timestep
        tau = 0.01

        !initial values in SI units
        y = 0.015 !position @ first data point
        v = 0.167 !velocity @ first data point 
        t = 0.492 !time     @ first data point
        a = 0.661 !calculate using formula

        open(unit=100, file='y_vs_t.dat')
        open(unit=200, file='v_vs_t.dat')
        open(unit=300, file='a_vs_t.dat')

        write (100,*) '#t-values		y-values'
        write (200,*) '#t-values		v-values'
        write (300,*) '#t-values		a-values'

        do while (t <= 2.4)
                v_temp = v
                v = v + tau*a

                y = y + tau*(v+v_temp)/2
                t = t + tau

                write(100,*) t, y
                write(200,*) t, v
                write(300,*) t, a
        end do

        close(100)
        close(200)
        close(300)

        end program eulermethod
