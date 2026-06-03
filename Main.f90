PROGRAM Main
    USE time_dependent_solver

    IMPLICIT NONE
    INTEGER ::  l_idx, nu
    INTEGER :: file, ii
    DOUBLE PRECISION :: sigma0, E_c, m, c_speed, q, Lambda, PI, K_constant, n_E, steady, D, P_in, Heat_tot, Ion_tot
    DOUBLE PRECISION :: E_min, E_max, dE, start_time, time_total, dt, daytoseconds
    DOUBLE PRECISION, ALLOCATABLE :: E(:), z_init(:), z_new(:), P_heat(:), P_ion(:), c2_array(:)
    DOUBLE PRECISION, DIMENSION(1) :: n_i, I, E_bar

    ! Defining constants
    daytoseconds = 86400.0

    sigma0 = 1.5e-15 !cm^2
    E_c = 1.6e-8 !erg

    m = 9.11e-28 !g
    c_speed = 2.9e10 !cm/s
    q = 4.8e-10 !statC
    Lambda = 20.0 
    PI = 3.1415
    K_constant = 4.0*PI*q**4/m*Lambda !cm^6 g s^{-4}

    I = [3.4e-11] !erg
    E_bar = 0.6*I !erg
    
    !Defininf the energy array
    E_min = 1.6e-13 !erg
    E_max = 1.6e-12*1e6 !erg
    dE = 1.6e-12*10.0 !erg
    nu = int((E_max - E_min)/dE)

    ALLOCATE(E(nu))
    E = [(E_min + (l_idx-1)*dE, l_idx=1, nu)]

    ALLOCATE(z_init(size(E)), z_new(size(E)))
    z_init = 0.0
    z_new = z_init
    
    ALLOCATE(P_heat(size(E)), P_ion(size(E)))
    P_heat = 0.0
    P_ion = 0.0

    ALLOCATE(c2_array(size(E)))
    c2_array = 0.0

    !time in seconds
    start_time = 1.0*daytoseconds 
    time_total = start_time
    
    !If steady = 0.0, the time-dependent terms are turned off (steady-state)
    !If steady = 1.0, the time-dependent terms are turned on (time-dependent)
    steady = 1.0

    !time step
    dt = 1.0*daytoseconds
    
    open(newunit=file, file="output_file.txt", status="replace", action="write")
   
    do while (time_total <= 1000.0*daytoseconds)

        !Homologous expansion of the number densities
        n_E = 2e8*(time_total/daytoseconds)**(-3)
        n_i = [1e8*(time_total/daytoseconds)**(-3)]
        
        D = 0.0
        P_in = 0.0

        call Numerical(E, E_bar, I, n_i, n_E, z_new, z_init, sigma0, E_c, m, &
             c_speed, K_constant, dt, E_max, steady, time_total, D, P_in, P_heat, P_ion, Heat_tot,&
             Ion_tot, c2_array)

        z_init = z_new

        !Adaptive time step
        dt = 0.2*time_total
        time_total = time_total + dt

        do ii=1, nu
            write(file, *) E(ii), z_new(ii), D, P_in, Heat_tot, Ion_tot, time_total
        end do

        print *, 'Finished an iteration'
        print *, time_total/daytoseconds
    end do

    close(file)

end program