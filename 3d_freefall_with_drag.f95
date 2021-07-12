program freefallInArray 
        !solve 3d freefall movement with array math

        implicit none

        integer(4) :: i
        real(8) :: pi, tau, time, maxHeight, drag_co, area, density, mass, v_avg, v_max, v_min
	real(8), allocatable, dimension(:) :: r, v, v_tmp, a, v_wind, v_real, g

        tau = 0.01
        pi = 4.0*atan(1.0)

        drag_co = 0.5               !drag coefficient for sphere 
        area = pi*((0.02541/2)**2)  !cross sectional area for baseball (m)
        density = 1.292             !density of air (kg/m^3)
        mass = 0.06801              !mass of baseball (kg)

        allocate(r(3), v(3), a(3), v_wind(3), v_real(3), g(3))

        r = (/0.0, 0.0, 1.0/)
        v = (/0.0, real(5.362*cos(24.0*pi/180.0)), real(5.362*sin(24.0*pi/180.0))/)
        g = (/0.0, 0.0, -9.8/) 
        v_wind = (/0.0, 0.0, 0.0/)

        open(unit=100, file='position.dat') 
        open(unit=200, file='velocity.dat') 
        open(unit=300, file='acceleration.dat') 

        maxHeight = 0.0
        v_avg = 0.0
        v_max = 0.0
        v_min = 10000.0
        time = 0.0
        do while (r(3) > 0)
                v_tmp = v
                v_real = v - v_wind
                a = g + (-1.0/2.0) * (drag_co*area*density) * norm2(v_real)*v/ mass
                v = v + tau*a
                r = r + tau*(v+v_tmp/2)
                
                time = time + tau !apparently splot only takes the first three lines, so time has to be at the end
                write(100,*) r, time
                write(200,*) v, time
                write(300,*) a, time

                if (maxHeight < r(3)) then
                        maxHeight = r(3)
                end if
                if (v_max < v(2)) then 
                        v_max = v(2)
                end if
                if (v_min > v(2)) then
                        v_min = v(2)
                end if
                v_avg = v_avg + v(2)
        end do
        write(*,*) 'maxHeight: ', maxHeight
        write(*,*) 'max v:     ', v_max
        write(*,*) 'min v:     ', v_min
        write(*,*) 'avg v:     ', v_avg/(time/tau) 
        write(*,*) 'time:      ', time

        close(100)
        close(200)
        close(300)

        end program freefallInArray
