MODULE time_dependent_solver
    IMPLICIT NONE

    PRIVATE
    PUBLIC S_E, f, g, beta, gamma_func, Numerical

    contains

    !Time-dependent source function
    !E_j    energy of the electron at bin j in erg
    !E_max  maximum energy in erg
    !time   time in seconds
    FUNCTION S_E(E_j, E_max, time) result(s_result)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: E_j, E_max, time
        DOUBLE PRECISION :: s_result
        
        if (E_max*0.967< E_j .and. E_j < E_max) then
            s_result = 1.0 / (E_max*0.033)*(time/86400.0)**(-1.5)
        else
            s_result = 0.0
        end if
    END FUNCTION

    !Cross section primaries
    !E_prime    primary electron energy (E') in erg
    !E_j        energy of the electron at bin j in erg
    !I_i        ionisation potential in erg
    !E_bar      energy scale term in erg
    FUNCTION f(E_prime, E_j, I_i, E_bar) result(f_res)
        IMPLICIT NONE 
        DOUBLE PRECISION, INTENT(IN) :: E_prime, E_j, I_i, E_bar
        DOUBLE PRECISION :: f_res, num, den

        num = 1.0
        den = 1.0 + ((E_prime - E_j - I_i)/E_bar)**2

        f_res = num/den
    END FUNCTION

    !Cross section secondaries
    !E_j    energy of the electron at bin j in erg
    !E_bar  energy scale term in erg
    FUNCTION g(E_j, E_bar) result(g_res)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: E_j, E_bar
        DOUBLE PRECISION :: g_res, num, den
        
        num = 1.0
        den = (1.0 + (E_j/E_bar)**2)

        g_res = num/den
    END FUNCTION

    !Normalisation constant of the cross-section
    !sigma0     scale of the normalisation constant in cm^2
    !E_prime    energy of the primary (E') in erg
    !I_i        ionisation potential in erg
    !E_bar      energy scale term in erg
    FUNCTION C(sigma0, E_prime, I_i, E_bar) result(C_res)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: sigma0, E_prime, I_i, E_bar
        DOUBLE PRECISION :: C_res

        C_res = sigma0/(E_bar*atan((E_prime-I_i)/(E_bar)))
    END FUNCTION

    !High energy correction term for the cross-section
    !E_prime    energy of the primary electron (E') in erg
    !E_c        characteristic cross-over energy in erg
    FUNCTION highE_correction(E_prime, E_c) result(corr_res)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: E_prime, E_c
        DOUBLE PRECISION :: corr_res

        corr_res = 1/(1+E_prime/E_c)
    END FUNCTION

    !beta coefficient for relativistic treatment
    !E_j        energy of the electron at bin j in erg
    !mass       electron mass in g
    !c_speed    speed of light in cm/s
    FUNCTION beta(E_j, mass, c_speed) result(beta_res)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: E_j, mass, c_speed
        DOUBLE PRECISION :: beta_res

        beta_res = sqrt(1 - (1 /(1+E_j/(mass*c_speed**2))**2))
    END FUNCTION

    !gamma coefficient for relativistic treatment
    !E_j        energy of the electron at bin j in erg
    !mass       electron mass in g
    !c_speed    speed of light in cm/s
    FUNCTION gamma_func(E_j, mass, c_speed) result(gamma_res)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: E_j, mass, c_speed
        DOUBLE PRECISION :: gamma_res

        gamma_res = 1 + E_j/(mass*c_speed**2)
    END FUNCTION

    !The numerical calculation of z(E,t) using energy and time discretisation
    !Energy         array of energies in erg
    !E_scale        array of scaling energies (\overline{E}) in erg
    !Ionisation     array of ionisation potentials in erg
    !number_dens_i  atom number density in cm^{-3}
    !number_dens_e  electron number density in cm^{-3}
    !K_constant     heat loss constant (4 \pi q^4/m_e*\Lambda) in cm^6 g s^{-4}
    !z_new          flux z at energy j and time step m+1 (z_j^{m+1}) in erg^{-1} cm^{-2} s^{-1}
    !z_initial      flux at energy j and time step m (z_j^m) in erg^{-1} cm^{-2} s^{-1}
    !sigma0         cross-section normalisation constant cm^2
    !E_c            characteristic cross-over energy for high energy corrections in erg
    !mass           electron mass in g
    !c_speed        speed of light in cm/s
    !dtime          time step in s
    !E_max          maximum energy in erg
    !steady         defines if one calculates steady-state or not 
    !time           time in s
    !D              deposition function (D(t)) in erg s^{-1} cm^{-3}
    !P_in           total injected power (P_in(t)) in erg s^{-1} cm^{-3}
    !P_heat         heating (L_heat/(beta^2*c^2)) in erg cm^{-1}
    !P_ion          ionisation (c_2*I_i) in erg cm^{-1}
    !Heat_tot       total heating rate erg s^{-1} cm^{-3}
    !Ion_tot        total ionisation rate erg s^{-1} cm^{-3}
    !c2_array       array for ionisation loss (c_2) 
    SUBROUTINE Numerical(Energy, E_scale, Ionisation, number_dens_i, number_dens_e, z_new, z_initial, sigma0, &
                         E_c, mass, c_speed, K_constant, dtime, E_max, steady, time, D, P_in, P_heat, P_ion, Heat_tot, Ion_tot, &
                        c2_array)
        IMPLICIT NONE
        INTEGER :: j_idx, NE, k_idx, i_idx, i_idx2
        DOUBLE PRECISION :: Ii, E_bari
        DOUBLE PRECISION :: c1, c2, c3, c4, LHS, Ej, Ep, integral1, integral2, dEnergy, n_idx
        DOUBLE PRECISION :: E_low, E_high, E_start, dE_eff, frac, Ep_mid, z_mid, z_start, Ej_plus
        DOUBLE PRECISION, INTENT(IN) :: Energy(:), E_scale(:), Ionisation(:), number_dens_i(:), &
                                        number_dens_e, K_constant, sigma0, E_c, mass, c_speed, dtime, E_max, steady, time
        DOUBLE PRECISION, INTENT(OUT) :: D, P_in, Heat_tot, Ion_tot
        DOUBLE PRECISION, INTENT(INOUT):: P_heat(:), P_ion(:), c2_array(:)
        DOUBLE PRECISION, INTENT(INOUT) :: z_new(:), z_initial(:)

        !Delta E
        dEnergy = Energy(2) - Energy(1)

        !Number of bins
        NE = size(Energy)

        !Finite difference from highest to lowest energy
        do j_idx = NE-1, 1, -1
            Ej = Energy(j_idx)
            Ej_plus = Energy(j_idx+1)

            !c1 time discretisation term
            c1 = 1/(beta(Ej, mass, c_speed)*c_speed)

            !Integral of the downscattered primaries
            integral1 = 0.0

            !Integral of the creation of secondaries
            integral2 = 0.0

            !calculation of the integrals using the midpoint method on a clipped domain
            do k_idx = j_idx+1, NE-1 
                !Energy of the integration variable E'
                Ep = Energy(k_idx)

                !summing over the ith atoms and atom number densities
                do i_idx = 1, size(number_dens_i)
                    n_idx = number_dens_i(i_idx)

                    Ii = Ionisation(i_idx)
                    E_bari = E_scale(i_idx)

                    E_low  = Energy(k_idx)
                    E_high = Energy(k_idx+1) 

                    if (E_high > (Ii + Ej)) then

                        !Effective integration width
                        E_start = max(E_low, Ii + Ej)
                        dE_eff = E_high - E_start

                        !Midpoint method
                        Ep_mid = 0.5*(E_start + E_high)
                        z_start = z_new(k_idx) + (z_new(k_idx+1)-z_new(k_idx)) * (E_start - E_low)/(E_high - E_low)
                        z_mid   = 0.5*(z_start + z_new(k_idx+1)) 

                        if (dE_eff > 0.0) then

                            !Integral of the downscattered primaries
                            integral1 = integral1 + n_idx * z_mid * f(Ep_mid, Ej, Ii, E_bari)&
                                        *C(sigma0, Ep_mid, Ii, E_bari)*highE_correction(Ep_mid, E_c)* dE_eff

                            !Integral of the creation of secondaries
                            integral2 = integral2 + n_idx* z_mid* g(Ej, E_bari)&
                                        *C(sigma0, Ep_mid, Ii, E_bari)*highE_correction(Ep_mid, E_c)* dE_eff
                        end if
                    end if

                end do
            
            end do

            !ionisation loss term
            c2 = 0.0
            
            !summing over the ith atoms and atom number densities
            do i_idx2 = 1, size(number_dens_i)
                Ii = Ionisation(i_idx2)
                n_idx = number_dens_i(i_idx2)
                E_bari = E_scale(i_idx2)

                E_low  = Ej
                E_high = Energy(j_idx+1)
                
                if (E_high <= Ii) then
                    frac = 0.0
                elseif (E_low >= Ii) then
                    frac = 1.0
                else
                    frac = (E_high - Ii) / (E_high - E_low)
                end if

                c2 = c2 + n_idx * sigma0*highE_correction(Ej, E_c)*frac
            end do

            c2_array(j_idx) = c2

            !One part of the heat loss
            c3 = K_constant*number_dens_e/(beta(Ej, mass, c_speed)**2*c_speed**2)&
                *(1.0 / dEnergy + 2.0 /(beta(Ej, mass, c_speed)**2 *mass*c_speed**2 *gamma_func(Ej, mass, c_speed)**3)) 
            
            !LHS of the equation
            LHS = c1/dtime*steady + c2 + c3

            !Other part of the heat loss
            c4 = K_constant*number_dens_e/(beta(Ej, mass, c_speed)**2*c_speed**2)* z_new(j_idx+1)/ dEnergy 

            !Calculation of z
            z_new(j_idx) = (S_E(Ej,E_max, time) + integral1 + integral2 + c4 + z_initial(j_idx)*c1/dtime*steady) / LHS

            !Calculating the Heating
            P_heat(j_idx) = K_constant*number_dens_e/(beta(Ej, mass, c_speed)**2*c_speed**2)
                
            !Calculating th Ionisation
            do i_idx2 = 1, size(number_dens_i)
                Ii = Ionisation(i_idx2)
            
                P_ion(j_idx) = c2_array(j_idx)*Ii
            end do
            
        end do

        !Initialising the total heating and ionisation rate, the deposition function and the total injected power
        Heat_tot = 0.0
        Ion_tot = 0.0
        D = 0.0
        P_in = 0.0

        do j_idx = 1, NE, 1

            D = D + z_new(j_idx)*(P_heat(j_idx) + P_ion(j_idx))* dEnergy

            P_in = P_in + S_E(Energy(j_idx), E_max, time)*Energy(j_idx)*dEnergy

            Heat_tot = Heat_tot + z_new(j_idx)*P_heat(j_idx) * dEnergy

            Ion_tot = Ion_tot + z_new(j_idx)*P_ion(j_idx) * dEnergy
        end do

        return

    END SUBROUTINE Numerical

END MODULE

