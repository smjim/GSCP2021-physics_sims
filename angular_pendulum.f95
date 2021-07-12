program angular_pendulum
        !solve non-linear oscillation with verlet method 
        !plot Kinetic energy, Potential energy, and Energy in addition to theta, omega, alpha

        implicit none

        integer(4) :: i, n
        real(8) time, gravity, length, object_m, tau, period, inertia, pi
	real(8), allocatable, dimension(:) :: theta, omega, alpha, kinetic, potential, total

        gravity = 9.8
        length = 1.0
        object_m = 1.0
        tau = 0.01

        inertia = object_m * (length**2)
        pi = 4.0*atan(1.0)
        n = 3000 !n=100 -> 1s long
        allocate(theta(n), omega(n), alpha(n), kinetic(n), potential(n), total(n))

        theta(1) = 130.0*(pi/180.0)
        omega(1) = 0.0
        alpha(1) = -(gravity/length) * sin(theta(1))

        theta(2) = theta(1) + tau*theta(1) + (1.0/2.0)*(tau**2.0)*alpha(1) !+error(t^4)

        open(unit=100, file='theta.dat') 
        open(unit=200, file='omega.dat') 
        open(unit=300, file='alpha.dat') 
        open(unit=400, file='kinetic.dat') 
        open(unit=500, file='potential.dat') 
        open(unit=600, file='total.dat') 

        do i=2, n-1
                !calculate acceleration to use as x''(i)
                alpha(i) = -(gravity/length) * sin(theta(i))

                !calculate the vector position using verlet method
                theta(i+1) = 2.0*theta(i) - theta(i-1) + (tau**2.0)*alpha(i) !+error(tau^4)

                !calculate omega based off of change in theta over t
                omega(i) = (theta(i+1)-theta(i-1))/(2.0*tau) !+error(tau^2)

                !calculate energies based off of omega quantities
                kinetic(i) = (1.0/2.0) * inertia * (omega(i)**2.0)
                potential(i) = object_m * gravity * length * (1.0-cos(theta(i)))
                total(i) = kinetic(i) + potential(i)
        end do

        do i=1, n-1
                time = i*tau + 0.0
                !write values to data files
                write(100,*) time, theta(i)
                write(200,*) time, omega(i)
                write(300,*) time, alpha(i)
                write(400,*) time, kinetic(i)
                write(500,*) time, potential(i)
                write(600,*) time, total(i)
        end do

        close(100)
        close(200)
        close(300)
        close(400)
        close(500)
        close(600)

        end program angular_pendulum
