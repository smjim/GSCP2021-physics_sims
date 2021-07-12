program charged_particle
        !simulation of trajectory of a charged particle in a magnetic field
        implicit none

        integer(4) :: i, n
        real(8) :: time, tau, mass, q, k
        real(8), allocatable, dimension(:) :: r, v, a, state, tmp, k1, k2, k3, k4, E, B, force

        tau = 0.01 !timestep
        time = 0.0 !time (seconds)
        n = 100000 !iterations

        mass = 1.0        !mass of charged particle (kg)
        k = 9.0*(10**(-9))!coulomb constant
        q = 1.0           !charge of particle (C)

        allocate(r(3), v(3), a(3), state(6), tmp(6), k1(6), k2(6), k3(6), k4(6), E(3), B(3), force(3))

        open(unit=100, file='p.dat')
        open(unit=200, file='v.dat')
        open(unit=300, file='a.dat')

        r = (/0.0, 0.0, 0.0/)  !distance (m)
        v = (/1.0, 0.0, 0.0/)  !velocity (m/s)
        E = (/0.0, 0.01, 0.1/) !energy field(?) (N/C)
        B = (/0.0, 0.0, 0.1/)  !magnetic field (T)

        do i=1, n-1 
                state = (/r, v/)

                !formula for calculating force on charged particle from magnetic field
                force = E + q*(/v(2)*B(3) - v(3)*B(2), -v(1)*B(3) - v(3)*B(1), v(1)*B(2) - v(2)*B(1)/)
                !force = E + q*cross_product(v, B)
                a = force/mass
                k1 = (/v, a/)
                tmp = state + (1.0/2.0)*tau*k1

                r = tmp(1:3)
                v = tmp(4:6)
                force = E + q*(/v(2)*B(3) - v(3)*B(2), -v(1)*B(3) - v(3)*B(1), v(1)*B(2) - v(2)*B(1)/)
                a = force/mass
                k2 = (/v, a/)
                tmp = state + (1.0/2.0)*tau*k2

                r = tmp(1:3)
                v = tmp(4:6)
                force = E + q*(/v(2)*B(3) - v(3)*B(2), -v(1)*B(3) - v(3)*B(1), v(1)*B(2) - v(2)*B(1)/)
                a = force/mass
                k3 = (/v, a/)
                tmp = state + (1.0/2.0)*tau*k3

                r = tmp(1:3)
                v = tmp(4:6)
                force = E + q*(/v(2)*B(3) - v(3)*B(2), -v(1)*B(3) - v(3)*B(1), v(1)*B(2) - v(2)*B(1)/)
                a = force/mass
                k4 = (/v, a/)

                state = state + (tau/6.0)*(k1 + (2.0*k2) + (2.0*k3) + k4)
                r = state(1:3)
                v = state(4:6)


                time = time + tau 
                write(100,*) r, time
                write(200,*) v, time
                write(300,*) a, time
        end do             

        close(100)
        close(200)
        close(300)

        end program charged_particle 
