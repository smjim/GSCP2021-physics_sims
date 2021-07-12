program freefallInArray 
        !solve 3d freefall movement with array math

        implicit none

        integer(4) :: i
        real(8) :: pi, tau, time, maxHeight, v_avg, minSpeed, maxSpeed
	real(8), allocatable, dimension(:) :: r, v, v_tmp, a

        tau = 0.01
        pi = 4.0*atan(1.0)

        allocate(r(3), v(3), a(3))

        r = (/0.0, 0.0, 1.172/)
        v = (/0.0, real(6.618*cos(24.0*pi/180.0)), real(6.618*sin(24.0*pi/180.0))/)
        a = (/0.0, 0.0, -9.8/) 

        open(unit=100, file='position.dat') 
        open(unit=200, file='velocity.dat') 
        open(unit=300, file='acceleration.dat') 

        maxHeight = 0.0
        minSpeed = 100000.0
        maxSpeed = 0.0
        v_avg = 0.0
        time = 0.0
        do while (r(3) > 0)
                v_tmp = v
                a = (/0.0, 0.0, -9.8/) 
                v = v + tau*a
                r = r + tau*(v+v_tmp/2)
                
                time = time + tau !apparently splot only takes the first three lines, so time has to be at the end
                write(100,*) r, time
                write(200,*) v, time
                write(300,*) a, time

                v_avg = v_avg + v(2)
                if (r(3) > maxHeight) then
                        maxHeight = r(3)
                end if
                if (v(2) > maxSpeed) then
                        maxSpeed = v(2)
                end if
                if (v(2) < minSpeed) then
                        minSpeed = v(2)
                end if
        end do
        v_avg = v_avg/(time/tau)

        write(*,*) 'max height:       ', maxHeight
        write(*,*) 'average velocity: ', v_avg
        write(*,*) 'maximum speed:    ', maxSpeed
        write(*,*) 'min speed:        ', minSpeed
        write(*,*) 'final distance:   ', r(2) 

        close(100)
        close(200)
        close(300)

        end program freefallInArray
