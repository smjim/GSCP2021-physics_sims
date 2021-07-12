program heat_equation 
        !heat equation simulation using rk4 and 
        implicit none

        integer(4) :: i, n
        real(8) :: time, tau, kinetic, potential, total, mass, GM, r_mag
        real(8), allocatable, dimension(:) :: r, v, a, state, tmp, k1, k2, k3, k4

        tau = 0.001 !0.1 years ~= 36 days
        time = 0.0  !time (earth years)
        n = 10010    !iterations

        mass = 1.0
        GM = -4.0*((4.0*atan(1.0))**2)

        allocate(r(3), v(3), a(3), state(6), tmp(6), k1(6), k2(6), k3(6), k4(6))

        open(unit=100, file='p.dat')
        open(unit=200, file='v.dat')
        open(unit=300, file='a.dat')

        open(unit=400, file='kinetic.dat')
        open(unit=500, file='potential.dat')
        open(unit=600, file='total.dat')

        r = (/1.0, 0.0, 0.0/)            !distance (AU)
        v = (/0.0, 4.0*atan(1.0), 1.0/)  !velocity (AU/year)

        do i=1, n-1 
                state = (/r, v/)

                r_mag = norm2(r)
                a = GM/(r_mag**3.0)*r
                k1 = (/v, a/)
                tmp = state + (1.0/2.0)*tau*k1

                r = tmp(1:3)
                v = tmp(4:6)
                r_mag = norm2(r)
                a = GM/(r_mag**3.0)*r
                k2 = (/v, a/)
                tmp = state + (1.0/2.0)*tau*k2

                r = tmp(1:3)
                v = tmp(4:6)
                r_mag = norm2(r)
                a = GM/(r_mag**3.0)*r
                k3 = (/v, a/)
                tmp = state + (1.0/2.0)*tau*k3

                r = tmp(1:3)
                v = tmp(4:6)
                r_mag = norm2(r)
                a = GM/(r_mag**3.0)*r
                k4 = (/v, a/)

                state = state + (tau/6.0)*(k1 + (2.0*k2) + (2.0*k3) + k4)
                r = state(1:3)
                v = state(4:6)


                kinetic = 0.5*mass*norm2(v)**2
                potential = GM*mass/norm2(r)
                total = kinetic + potential


                time = time + tau 
                write(100,*) r, time
                write(200,*) v, time
                write(300,*) a, time

                write(400,*) time, kinetic
                write(500,*) time, potential
                write(600,*) time, total
        end do             

        close(100)
        close(200)
        close(300)
        close(400)
        close(500)
        close(600)

        end program heat_equation 
