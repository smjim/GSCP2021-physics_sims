program freefall_airdrag
        !calculates mechanical quantities of an object in freefall with air resistance using the euler chromer method
        implicit none

	real(8) :: y, v, a, t, tau, v_temp, drag_co, area, density, mass

        !timestep
        tau = 0.01

        !constants for air drag
        drag_co = 0.5   !drag coefficient for sphere
        !area = 3.1415   !cross sectional area for sphere of r=1
        area = 0.00419  !cross sectional area for a baseball (r=3.65cm) 
        !density = 1000.0!density of water
        density = 1.292 !density of air (kg/m^3)
        !mass = 1        !mass of sphere of r=1
        mass = 0.145    !mass of baseball (kg)

        !initial values in SI units
        y = 5700.0 !position @ t=0 example problem
        v = -3.8 !velocity @ t=0 example problem 
        t =  0.0 !time
        a = -9.8 !acceleration

        open(unit=100, file='y_vs_t.dat')
        open(unit=200, file='v_vs_t.dat')
        open(unit=300, file='a_vs_t.dat')

        write (100,*) '#t-values		y-values'
        write (200,*) '#t-values		v-values'
        write (300,*) '#t-values		a-values'

        do while (y >= 0)
                write(100,*) t, y
                write(200,*) t, v
                write(300,*) t, a
                
                v_temp = v
                v = v + tau*a

                y = y + tau*(v+v_temp)/2
                t = t + tau

                !air drag equation
                a = -9.8 + (drag_co * area * density * (v**2.0)/2.0)
        end do

        close(100)
        close(200)
        close(300)

        end program freefall_airdrag
