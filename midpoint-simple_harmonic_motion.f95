program simple_harmonic_motion
        !solve simple harmonic motion with midpoint method and arrays
        !plot Kinetic energy, Potential energy, and Energy plus x, v, a

        implicit none

        integer(4) :: i, n
        real(8) spring_k, object_m, tau, period
	real(8), allocatable, dimension(:) :: x, v, a, kinetic, potential, total

        spring_k = 0.5
        object_m = 1.0
        tau = 0.01

        period = 2.0*(4.0*atan(1.0)) * (object_m/spring_k)**(1.0/2.0)

        n = 5*floor(period/tau)
        !n = 1000
        allocate(x(n), v(n), a(n), kinetic(n), potential(n), total(n))
        x(0) = 10
        v(0) = 0

        open(unit=100, file='position.dat') 
        open(unit=200, file='velocity.dat') 
        open(unit=300, file='acceleration.dat') 
        open(unit=400, file='kinetic.dat') 
        open(unit=500, file='potential.dat') 
        open(unit=600, file='total.dat') 

        do i=1, n-1
                !calculate the 1d quantities
                a(i) = -spring_k*x(i-1)/object_m
                v(i) = v(i-1) + tau*(a(i-1)+a(i))/2.0
                x(i) = x(i-1) + tau*(v(i-1)+v(i))/2.0
                
                !calculate energies
                kinetic(i) = 1.0/2.0 * object_m * v(i)**2
                potential(i) = 1.0/2.0 * x(i)**2 
                total(i) = kinetic(i) + potential(i)

                !write values to data files
                write(100,*) i*tau, x(i)
                write(200,*) i*tau, v(i)
                write(300,*) i*tau, a(i)
                write(400,*) i*tau, kinetic(i)
                write(500,*) i*tau, potential(i)
                write(600,*) i*tau, total(i)
        end do

        close(100)
        close(200)
        close(300)
        close(400)
        close(500)
        close(600)

        end program simple_harmonic_motion
