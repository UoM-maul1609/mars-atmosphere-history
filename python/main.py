import numpy as np
import copy
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import crater_chronology
import impactor_model
import impactor_velocity
import impactor_atmos_loss_gain_eq4_5
import svp_ice
import volcano_outgassing
import sputtering_co2
import MarsSputtering_Calc6_7
import isotopic_data
import photo_chemical
import IDPs
import mars_temp_forget
import mars_temperature

"""
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
"""
def test(y,t,omega):
    y1, z = y # this is y and dy/dx
    dydt = [z,-omega**2*y1 ]
    return dydt



# if __name__ == "__main__":
def run_model(runNo, obliquity_flag, sputtering_flag=True, pce_flag=True,\
        X_gas=0.01,f_comet=0.001, C_vol=5.0, crater_model=2,C_Ne_IDP = 10.,\
	dynamo_time = 4.1, output_file=None, comet_time=4.1, D_max=None,
	impact_log=None, rate_log=None, f_sp=1.0, f_ph_N=1.0, X_gas_comet=None,
	P_upper_cap=None, T_atm=None):
    """Run one Monte Carlo history (seeded by runNo).

    Impacts are discrete events (Kurokawa §2.2, Fig. 12 caption): N =
    2e21 kg / <m> impactors are drawn from the truncated size distribution,
    each assigned a time from the crater chronology, a velocity, and (after
    comet_time) comet status with number ratio f_comet. Each is applied in
    its 1 Myr step using the atmosphere at that time. Pass a list as
    impact_log to receive (time_yr, D_m, v_m_s, is_comet, eroded_fraction)
    for every impact. Pass a list as rate_log to receive, per step,
    (time_yr, volcanic N2, sputtered N2, photochemical N2 [molecules/s],
    pN2 [bar], pCO2 [bar]) -- comparable with Kurokawa Fig. 12a.

    f_sp and f_ph_N are optional multipliers on all sputtering rates and on
    photochemical N escape (as f_sp and f_pN in Hu & Thomas 2023, who
    sample 0.5-2 and 0.01-20). The defaults of 1 are the nominal Kurokawa
    parameterisation. The nominal N lifetimes come out ~3x shorter than
    Kurokawa Fig. 14; see README.
    """
    """
        initial conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    
    Patm = 1.5 # initial pressure in most calculations - section 2.9, pp449
    tinit=-4.5e9 
    Rmars=3389.5e3 # m
    Mmars=6.39e23 # kg
    Ggrav=6.67e-11 # gravitational constant
    grav_mars=3.721 # m/s^2
    rho_t=4000. # density of martian surface
    # Atmospheric temperature for the impact scale height and surface
    # density. Kurokawa §2.1 uses one climate assumption: 273.16 K before
    # collapse and 210 K after; §4.3 lists "the assumed temperature" among the
    # impact-model inputs. The previous code fixed 210 K throughout.
    # Pass T_atm to force a constant value (e.g. 210. for the old behaviour).
    def _t_atm(is_collapsed):
        if T_atm is not None:
            return T_atm
        return 210. if is_collapsed else 273.16
    tmars = _t_atm(False)
    Rgas = 8.314 # ideal gas constant
    #MolW_atm = 44.01e-3 # molecular weight of martian atmosphere - may change???
    MolW_atm = np.array([44.01e-3, 28.0134e-3, 18.01528e-3])
    
    rho_pr=2600. # density of impactor
    Pcollapse = 0.5  # bar; overwritten per step in the obliquity modes
    collapsed = False  # CO2 collapsed into surface reservoirs?
    printOutput = False
    
    total_impactor_mass = 2.e21 # kg
    tfinal=0.
    tstep=1e6
    omega=2.*np.pi*1.e-9
    output_interval = 1e7
    last_output = tinit-output_interval
    
    # IDP rates, table 2
    IDP_rate_moles = IDPs.get_mole(['Ne','Ar','Kr','Xe'],Ne=C_Ne_IDP)
    
    isotope_comps_sim = np.array(isotopic_data.isotope_comps['Volcanic degassing']).copy()
    ind,=np.where(['Xe' in iso for iso in isotopic_data.isotopes])
    isotope_comps_sim[ind]=np.array(isotopic_data.isotope_comps['Mars'])[ind].copy()
    isotopes = [None]*5
    isotopes[0] = [isotopic_data.isotopes[0]]
    isotopes[1] = [isotopic_data.isotopes[1]]
    isotopes[2] = [isotopic_data.isotopes[2]]
    isotopes[3] = isotopic_data.isotopes[3:8]
    isotopes[4] = isotopic_data.isotopes[8:]
    isotopes_sim = [None]*5
    isotopes_sim[0] = [(isotope_comps_sim[0]/1000.+1)*3676.5e-6]
    isotopes_sim[1] = [isotope_comps_sim[1]]
    isotopes_sim[2] = [isotope_comps_sim[2]]
    isotopes_sim[3] = isotope_comps_sim[3:8]
    isotopes_sim[4] = isotope_comps_sim[8:]
    isotope_sources = dict()
    isotope_sources['Mars'] = [None]*5
    isotope_sources['Volcanic degassing'] = [None]*5
    isotope_sources['Asteroids'] = [None]*5
    isotope_sources['Comets'] = [None]*5
    isotope_sources['IDPs'] = [None]*5
    isotope_sources['Mars'][0] = [(isotopic_data.isotope_comps['Mars'][0]/1000.+1)*3676.5e-6]
    isotope_sources['Mars'][1] = [isotopic_data.isotope_comps['Mars'][1]]
    isotope_sources['Mars'][2] = [isotopic_data.isotope_comps['Mars'][2]]
    isotope_sources['Mars'][3] = isotopic_data.isotope_comps['Mars'][3:8]
    isotope_sources['Mars'][4] = isotopic_data.isotope_comps['Mars'][8:]
    
    isotope_sources['Volcanic degassing'][0] = [(isotopic_data.isotope_comps['Volcanic degassing'][0]/1000.+1)*3676.5e-6]
    isotope_sources['Volcanic degassing'][1] = [isotopic_data.isotope_comps['Volcanic degassing'][1]]
    isotope_sources['Volcanic degassing'][2] = [isotopic_data.isotope_comps['Volcanic degassing'][2]]
    isotope_sources['Volcanic degassing'][3] = isotopic_data.isotope_comps['Volcanic degassing'][3:8]
    isotope_sources['Volcanic degassing'][4] = isotopic_data.isotope_comps['Volcanic degassing'][8:]
    
    isotope_sources['Asteroids'][0] = [(isotopic_data.isotope_comps['Asteroids'][0]/1000.+1)*3676.5e-6]
    isotope_sources['Asteroids'][1] = [isotopic_data.isotope_comps['Asteroids'][1]]
    isotope_sources['Asteroids'][2] = [isotopic_data.isotope_comps['Asteroids'][2]]
    isotope_sources['Asteroids'][3] = isotopic_data.isotope_comps['Asteroids'][3:8]
    isotope_sources['Asteroids'][4] = isotopic_data.isotope_comps['Asteroids'][8:]

    isotope_sources['Comets'][0] = [(isotopic_data.isotope_comps['Comets'][0]/1000.+1)*3676.5e-6]
    isotope_sources['Comets'][1] = [isotopic_data.isotope_comps['Comets'][1]]
    isotope_sources['Comets'][2] = [isotopic_data.isotope_comps['Comets'][2]]
    isotope_sources['Comets'][3] = isotopic_data.isotope_comps['Comets'][3:8]
    isotope_sources['Comets'][4] = isotopic_data.isotope_comps['Comets'][8:]

    isotope_sources['IDPs'][0] = [(isotopic_data.isotope_comps['IDPs'][0]/1000.+1)*3676.5e-6]
    isotope_sources['IDPs'][1] = [isotopic_data.isotope_comps['IDPs'][1]]
    isotope_sources['IDPs'][2] = [isotopic_data.isotope_comps['IDPs'][2]]
    isotope_sources['IDPs'][3] = isotopic_data.isotope_comps['IDPs'][3:8]
    isotope_sources['IDPs'][4] = isotopic_data.isotope_comps['IDPs'][8:]
    
    isotope_mass =[None]*5
    isotope_mass[0] = [isotopic_data.isotope_mass[0]]
    isotope_mass[1] = [isotopic_data.isotope_mass[1]]
    isotope_mass[2] = [isotopic_data.isotope_mass[2]]
    isotope_mass[3] = isotopic_data.isotope_mass[3:8]
    isotope_mass[4] = isotopic_data.isotope_mass[8:]
    element_mass = isotopic_data.element_mass[1:]
    # Recompute diffusion factors at each time step using Table 2 EUV.
    """
        ----------------------------------------------------------------------------------
    """




    """
        auxiliary calcs ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    Amars=4.*np.pi*Rmars**2 # surface area of mars
    Matm = Patm*1e5*Amars / grav_mars # mass of atmosphere
    # call function to set isotopes and adjust atmosphere to same pressure
    (mole1,mole_elements)=isotopic_data.set_amount_of_element(Matm/MolW_atm[0])
    Natm = np.concatenate([mole1,[0]]) # number of moles of co2, n2, h2o
    # Keep the prescribed initial pressure exact (the helper normalises by
    # number of molecules, rather than mass, so does not quite preserve it).
    initial_scale = Matm / np.dot(Natm, MolW_atm)
    Natm *= initial_scale
    mole_elements *= initial_scale
    # The helper returns approximate reference-isotope counts for C and N;
    # the dynamical elemental inventories must use total elemental counts.
    mole_elements[0] = Natm[0]
    mole_elements[1] = 2.0 * Natm[1]
    Matm = Natm *MolW_atm  # mass of co2, n2, h2o
#     Natm = Matm / MolW_atm     # number of moles of co2, n2, h2o
    # Kurokawa Section 2.1: warm (273.16 K) before collapse, 210 K after.
    Ninit=Natm[1]
    Hscale=Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars) # scale height - m - may change???
    t = np.mgrid[tinit:tfinal+tstep:tstep]

    (euv_max, f_co2_flux_max) = sputtering_co2.sputtering_co2_rate((t[0]-tinit)/1e9)
    nsteps = len(t)
    n_ode=11
    rng=np.random.RandomState(runNo)
    if X_gas_comet is None:
        X_gas_comet = isotopic_data.comet_co2_mass_fraction()
    """
        ----------------------------------------------------------------------------------
    """

    """
        set-up arrays ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    ystore=np.zeros((nsteps,n_ode))
    ystore[0,0]=Patm										# Patm
    ystore[0,1]=(isotopes_sim[0][0]/(3676.5e-6)-1.)*1000.	# d15N
    ystore[0,2]=isotopes_sim[1][0]							# 22Ne/20Ne
    ystore[0,3]=isotopes_sim[2][0]							# 38Ar/36Ar
    ystore[0,4]=isotopes_sim[3][4]							# 86Kr/84Kr
    ystore[0,5]=isotopes_sim[4][7]							# 136Xe/130Xe
    ystore[0,6]=mars_temp_forget.meanTvP1(Patm)				# tglobal
    ystore[0,7]=Hscale										# Hscale
    ystore[0,8]=Matm[0]										# co2
    ystore[0,9]=Matm[1]										# n2
    ystore[0,10]=Matm[2]										# h2o
    
    # Discrete bombardment history. Step weights follow the crater
    # chronology (D >= 20 km crater rate, as before).
    step_weights = np.array([crater_chronology.N20_between(
        Amars, -t[k+1]/1e9, -t[k]/1e9, model=crater_model)
        for k in range(nsteps-1)])
    # Default size cap: the diameter of one body carrying the whole mass
    # budget (~1140 km for 2e21 kg). This reproduces Kurokawa Fig. 3a
    # (~2.3e5 impacts >= 3 km, largest bodies ~300-800 km); the former
    # 2000 km cap gave only ~9e4 impacts.
    if D_max is None:
        D_max = impactor_model.budget_diameter(total_impactor_mass, rho_pr)
    imp_step, imp_diam, imp_mass = impactor_model.sample_impact_history(
        rng, step_weights, total_impactor_mass, rho_pr, D_max)
    imp_vel = impactor_velocity.inv_cdf1(rng.rand(len(imp_diam)))
    # Comet or asteroid: f_comet is the number ratio comets:asteroids after
    # comet_time, so P(comet) = f/(1+f).
    imp_time = t[imp_step]
    imp_comet = (imp_time > -comet_time*1e9) & \
        (rng.rand(len(imp_diam)) < f_comet/(1.+f_comet))
    step_bounds = np.searchsorted(imp_step, np.arange(nsteps))
    
    # escape velocity
    uesc=impactor_atmos_loss_gain_eq4_5.escape_velocity(Ggrav,Mmars,Rmars)
    
    """
        ----------------------------------------------------------------------------------
    """
    if printOutput:
    	print('impacts: %d, total mass %.3g kg' % (len(imp_mass), imp_mass.sum()))
    	print('Note, the pressure for atmospheric collapse is ' + str(Pcollapse))
    #nsteps=0

    """
        main loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    for i in range(nsteps-1):
        # 1-2. collapse and H2O are handled after degassing (below).

        # 3-4. discrete impacts falling in this timestep
        sel = slice(step_bounds[i], step_bounds[i+1])
        diams, mass, vels = imp_diam[sel], imp_mass[sel], imp_vel[sel]
        is_comet = imp_comet[sel]
        num_impacts = len(diams)
        if num_impacts:
            # sum(Matm) / sum(Natm) is the molecular weight
            tmars = _t_atm(collapsed)
            Hscale = Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars)
            rho_0 = Patm*1.e5/tmars/(Rgas/(np.sum(Matm)/np.sum(Natm)))
            xi=np.maximum(impactor_atmos_loss_gain_eq4_5. \
                dimensionless_xi(diams,rho_t,rho_pr, vels, uesc, Hscale, rho_0), 1.e-10)
            deltaM=impactor_atmos_loss_gain_eq4_5.change_atmos(xi,vels,uesc,mass)
            Xpr=np.maximum(impactor_atmos_loss_gain_eq4_5. \
                normalised_projectile_mass(rho_t,rho_pr,vels,uesc,xi),0)
            # Degassed CO2-equivalent mass. Asteroids: X_gas. Comets: their
            # own Fig. 2 volatile content (X_gas_comet, ~2.9 wt%); previously
            # comets also used X_gas, undersupplying all cometary volatiles
            # (C, N and noble gases) by ~3x for X_gas = 1%. The noble-gas and
            # N amounts follow from this CO2 mass via the comet ratios to C.
            retained = (1.-Xpr)*mass
            mass_added_comet = np.sum(retained[is_comet])*X_gas_comet
            mass_added_asteroid = np.sum(retained[~is_comet])*X_gas
            # Impacts within the step act in sequence on the remaining
            # atmosphere: remaining fraction = prod(1 - f_j), each f_j <= 0.99.
            M_now = np.sum(Natm*MolW_atm)
            f_j = np.minimum(deltaM / M_now, 0.99)
            const1 = 1. - np.prod(1. - f_j)
            if impact_log is not None:
                impact_log.extend(zip(np.full(num_impacts, t[i]), diams, vels,
                                      is_comet, f_j))
        else:
            mass_added_comet = mass_added_asteroid = 0.
            const1 = 0.
        mass_added = mass_added_comet + mass_added_asteroid

        # erosion is not isotopically fractionating (Kurokawa §2.7)
        Natm = Natm *(1.- const1)
        mole_elements = mole_elements*(1.-const1)

        # 4a add the carbon dioxide outgassed from asteroids and comets
        Natm[0] = Natm[0]+mass_added / MolW_atm[0]
        # molecular nitrogen
        Natm[1] = Natm[1]+0.5 * mass_added_asteroid / MolW_atm[0]* \
            isotopic_data.abundances1['Asteroids'][1][1] / \
            isotopic_data.abundances1['Asteroids'][1][0] * \
            isotopic_data.solar_abundances[1] / isotopic_data.solar_abundances[0] + \
            0.5*mass_added_comet / MolW_atm[0]* \
            isotopic_data.abundances1['Comets'][1][1] / \
            isotopic_data.abundances1['Comets'][1][0] * \
            isotopic_data.solar_abundances[1] / isotopic_data.solar_abundances[0]
            
        # also for the elemental composition:
        N=mole_elements[1:].copy()
        dN=[mass_added_asteroid / MolW_atm[0] * \
            isotopic_data.abundances1['Asteroids'][1][1] / \
            isotopic_data.abundances1['Asteroids'][1][0] * \
            isotopic_data.solar_abundances[1] / isotopic_data.solar_abundances[0], \
            mass_added_asteroid / MolW_atm[0] * \
            isotopic_data.abundances1['Asteroids'][1][2] / \
            isotopic_data.abundances1['Asteroids'][1][0] * \
            isotopic_data.solar_abundances[2] / isotopic_data.solar_abundances[0], \
            mass_added_asteroid / MolW_atm[0] * \
            isotopic_data.abundances1['Asteroids'][1][3] / \
            isotopic_data.abundances1['Asteroids'][1][0] * \
            isotopic_data.solar_abundances[3] / isotopic_data.solar_abundances[0], \
            mass_added_asteroid / MolW_atm[0] * \
            isotopic_data.abundances1['Asteroids'][1][4] / \
            isotopic_data.abundances1['Asteroids'][1][0] * \
            isotopic_data.solar_abundances[4] / isotopic_data.solar_abundances[0], \
            mass_added_asteroid / MolW_atm[0] * \
            isotopic_data.abundances1['Asteroids'][1][5] / \
            isotopic_data.abundances1['Asteroids'][1][0] * \
            isotopic_data.solar_abundances[5] / isotopic_data.solar_abundances[0] ]
            
        dN2=[mass_added_comet / MolW_atm[0] * \
            isotopic_data.abundances1['Comets'][1][1] / \
            isotopic_data.abundances1['Comets'][1][0] * \
            isotopic_data.solar_abundances[1] / isotopic_data.solar_abundances[0], \
            mass_added_comet / MolW_atm[0] * \
            isotopic_data.abundances1['Comets'][1][2] / \
            isotopic_data.abundances1['Comets'][1][0] * \
            isotopic_data.solar_abundances[2] / isotopic_data.solar_abundances[0], \
            mass_added_comet / MolW_atm[0] * \
            isotopic_data.abundances1['Comets'][1][3] / \
            isotopic_data.abundances1['Comets'][1][0] * \
            isotopic_data.solar_abundances[3] / isotopic_data.solar_abundances[0], \
            mass_added_comet / MolW_atm[0] * \
            isotopic_data.abundances1['Comets'][1][4] / \
            isotopic_data.abundances1['Comets'][1][0] * \
            isotopic_data.solar_abundances[4] / isotopic_data.solar_abundances[0], \
            mass_added_comet / MolW_atm[0] * \
            isotopic_data.abundances1['Comets'][1][5] / \
            isotopic_data.abundances1['Comets'][1][0] * \
            isotopic_data.solar_abundances[5] / isotopic_data.solar_abundances[0] ]
            
        mole_elements[0] = mole_elements[0]+mass_added / MolW_atm[0]
        mole_elements[1] = mole_elements[1]+dN[0] + dN2[0]
        mole_elements[2] = mole_elements[2]+dN[1] + dN2[1]
        mole_elements[3] = mole_elements[3]+dN[2] + dN2[2]
        mole_elements[4] = mole_elements[4]+dN[3] + dN2[3]
        mole_elements[5] = mole_elements[5]+dN[4] + dN2[4]
            
            
        Matm = Natm*MolW_atm
        
        # isotope composition from impact replenishment due to asteroids see equation 11
        isotopes_sim = \
            isotopic_data.impact_replenishment(isotopes_sim, \
                isotope_sources['Asteroids'], \
                isotope_sources['Comets'], N, dN, dN2)
        

        # 5. Sputtering and Photochemical escape
        # sputter each element / isotope: c (Co2), N (N2), Ne, Ar, Kr, Xe       
        (euv_flux, f_co2_flux) = sputtering_co2.sputtering_co2_rate((t[i]-tinit)/1e9)
        Rdiffs_isotopes = MarsSputtering_Calc6_7.calculate_rdiffs(
            element_mass, isotope_mass, euv_flux)
        F_i_sp1 = MarsSputtering_Calc6_7.F_i_sp(
            f_co2_flux, np.concatenate((Natm[:2], mole_elements[2:])), Natm[0],
            ['C','N','Ne','Ar','Kr','Xe'], euv=euv_flux)
        # Escape + continuous supply, integrated exactly over the step.
        # Per-particle escape rates k = F/N are taken at the start of the
        # step. In the collapsed era N2 and noble-gas lifetimes can be far
        # shorter than the 1 Myr step (~0.2 Myr for N2 at 4.1 Ga), so escape
        # and supply are solved TOGETHER for each isotope population
        # (isotopic_data.evolve_open_system). This replaces the former
        # sequence escape -> IDP supply -> volcanic supply with explicit Euler
        # steps and 50% caps, which stripped and over-fractionated small
        # inventories and then reset them to the source composition.
        dt = tstep * 86400 * 365
        NA = 6.02e23
        sp_inv = np.concatenate((Natm[:2], mole_elements[2:]))   # moles
        sputtering_on = (-t[i]/1e9 < dynamo_time) and sputtering_flag
        if sputtering_on:
            k_sp = f_sp * np.where(sp_inv > 0.,
                            np.maximum(F_i_sp1, 0.) / np.maximum(sp_inv*NA, 1e-300), 0.)
        else:
            k_sp = np.zeros(len(sp_inv))

        # 5b photochemical escape rates (per particle)
        if pce_flag:
            # Eq. 8: F_Lyc ~ t^-0.86 with t the SOLAR age (Hu et al.
            # 2015), not time since the model start.
            age0 = sputtering_co2.SOLAR_AGE_AT_START
            Fcph = photo_chemical.carbon_escape((t[i]-tinit)/1e9 + age0,
                (t[i+1]-tinit)/1e9 + age0,
                sputtering_co2.PRESENT_SOLAR_AGE)  # C atoms/s
            # No EUV cap: Kurokawa §2.4 scales these rates as F_EUV and
            # F_EUV^1.5, and their Fig. 14a lifetimes at 4.5 Ga (~300x shorter
            # than today) need the full ~100x saturated EUV. The previous
            # np.minimum(euv_flux, 6) suppressed pre-3.8 Ga N escape.
            (fn2_1,fn2_2,fn2_3)=photo_chemical.nitrogen_escape(
                np.maximum(Natm[0],1e-3),Natm[1],euv_flux,1.,Amars, Hscale)
            fn2 = fn2_1 + fn2_2 + fn2_3   # escaping N atoms/s
            k_ph_C = max(float(Fcph), 0.) / max(float(mole_elements[0])*NA, 1e-300)
            k_ph_N = f_ph_N * max(float(fn2), 0.) / max(float(mole_elements[1])*NA, 1e-300)
        else:
            k_ph_C = k_ph_N = 0.

        # Carbon (CO2): sputtering + photochemical escape. Volcanic CO2 is
        # added below, before the collapse and saturation adjustments.
        lost_C = Natm[0] * -np.expm1(-(k_sp[0] + k_ph_C) * dt)
        Natm[0] -= lost_C
        mole_elements[0] -= lost_C

        # Continuous sources over this step (ledger units: moles).
        h2o=volcano_outgassing.get_outgas_rate(t[i]/1e9,'h2o') * C_vol
        co2=volcano_outgassing.get_outgas_rate(t[i]/1e9,'co2') * C_vol
        n2 =volcano_outgassing.get_outgas_rate(t[i]/1e9,'n2')  * C_vol
        volc_N = 2.*n2*dt/NA                                  # N atoms
        volc_noble = co2*dt/NA * \
            isotopic_data.abundances1['Volcanic degassing'][1][2:] / \
            isotopic_data.abundances1['Volcanic degassing'][1][0] * \
            isotopic_data.solar_abundances[2:] / isotopic_data.solar_abundances[0]
        idp_noble = IDP_rate_moles*tstep / 1e9                # reference isotope

        # Nitrogen: sputtering of N2 (diffusive R) + photochemical escape
        # (R_ph,15/14 = R_diff,15/14 / 1.4), volcanic supply. Total-N basis.
        isotopes_sim[0], N_new, (sp_N, ph_N) = isotopic_data.evolve_open_system(
            isotopes_sim[0], mole_elements[1], [k_sp[1], k_ph_N],
            [Rdiffs_isotopes[0], np.asarray(Rdiffs_isotopes[0])/1.4],
            [(volc_N, isotope_sources['Volcanic degassing'][0])], dt, False)
        mole_elements[1] = N_new
        Natm[1] = N_new / 2.
        # Noble gases: sputtering, volcanic + IDP supply. Reference basis.
        for g in range(1, 5):
            isotopes_sim[g], inv_new, _ = isotopic_data.evolve_open_system(
                isotopes_sim[g], mole_elements[g+1], [k_sp[g+1]],
                [Rdiffs_isotopes[g]],
                [(volc_noble[g-1], isotope_sources['Volcanic degassing'][g]),
                 (idp_noble[g-1], isotope_sources['IDPs'][g])], dt, True)
            mole_elements[g+1] = inv_new
        Matm = Natm * MolW_atm
        if rate_log is not None:
            rate_log.append((t[i], n2, sp_N/2.*NA/dt, ph_N/2.*NA/dt,
                             Matm[1]*grav_mars/Amars/1e5,
                             Matm[0]*grav_mars/Amars/1e5))
        
        
        
        
        # 6. Volcanic degassing of CO2 and H2O (N2 and noble gases were
        # supplied in the open-system solve above).
        Matm = Natm * MolW_atm
        Matm = Matm + np.array([co2,0.,h2o])*MolW_atm*tstep*86400*365/6.02e23
        # Kurokawa §2.1: 273.16 K before collapse, 210 K after, in ALL modes.
        # (Previously the obliquity modes always used 273.16 K, leaving
        # ~6 mbar of water vapour in a collapsed atmosphere.)
        if obliquity_flag == -1:
            collapsed = Patm < 0.5
        water_temperature = 210. if collapsed else 273.16
        ph2o = svp_ice.svp_ice01(water_temperature)
        Matm[2] = np.minimum(Matm[2]*grav_mars/Amars, ph2o)*Amars/grav_mars
        Natm = Matm / MolW_atm
        
        
        
        
        # readjust Matm: firstly it cant go above the svp from Forget et al:
        # Optional upper CO2 limit (Forget et al. 2013 ~3 bar CO2-ice caps).
        # Kurokawa et al. (2018) only discuss this limit: their Fig. 5 shows
        # pressures up to ~6 bar (X_gas = 1%) and ~90 bar (X_gas = 10%), so
        # no cap is applied by default (P_upper_cap=None). The previous code
        # always capped CO2 at 2.75 bar, deleting delivered CO2.
        if P_upper_cap is not None and Matm[0]*grav_mars/Amars/1.e5 > P_upper_cap:
        	Matm[0]=P_upper_cap*Amars*1e5/grav_mars
        	Patm = np.sum(Matm)*grav_mars/(Amars)/1e5
        	mole_elements[0] = Matm[0] / MolW_atm[0]

        # Total pressure consistent with the state (no artificial 6 mbar floor
        # on Patm alone: that left Patm and Matm disagreeing).
        Patm = np.sum(Matm)*grav_mars/Amars/1.e5
        # Obliquity modes (Bradley MPhil §4.3.1): collapse is triggered when the
        # CO2 partial pressure exceeds the CO2 saturation vapour pressure at the
        # WRF 5th-percentile surface temperature for the current total pressure.
        # After collapse the CO2 partial pressure is held at the present-day
        # 6 mbar buffered by surface reservoirs (Kurokawa §2.1; thesis control
        # run). Previously CO2 was instead clipped to svp(T5) itself, which
        # drove pCO2 to ~2 mbar at 90 deg, compared total (not CO2) pressure,
        # and never flagged a collapsed state.
        if obliquity_flag in (1, 2, 3) and not collapsed:
            T5 = {1: mars_temp_forget.minTvP_bradley_0,
                  2: mars_temp_forget.minTvP_bradley_45,
                  3: mars_temp_forget.minTvP_bradley_90}[obliquity_flag](Patm)
            Pcollapse = mars_temperature.co2_vap_press(T5)/1e5
            pCO2 = Matm[0]*grav_mars/Amars/1.e5
            if pCO2 > Pcollapse:
                collapsed = True
        if collapsed and obliquity_flag != -1:
            Matm[0] = 6.e-3*1e5*Amars / grav_mars
            Patm = np.sum(Matm)*grav_mars/(Amars)/1e5
            mole_elements[0] = Matm[0] / MolW_atm[0]
        elif obliquity_flag == -1:
            Pcollapse = 0.5
        if obliquity_flag == -1 and (Patm < Pcollapse):
            Patm = 6.e-3
            Matm[0] = Patm*1e5*Amars / grav_mars # mass of atmosphere
            Patm = np.sum(Matm)*grav_mars/(Amars)/1e5
            #Matm[1:] = 0. # if it's collapsed, just assume it is CO2
            mole_elements[0] = Matm[0] / MolW_atm[0]

        
        
        Natm = Matm / MolW_atm
        # Keep the carbon ledger synchronised with atmospheric CO2 after the
        # collapse and saturation adjustments.
        mole_elements[0] = Natm[0]

        # update scale height
        tmars = _t_atm(collapsed)
        Hscale=Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars)
        Matm = Natm *MolW_atm
        # calculate the temperature global mean from Forget et al.
        tglobal = mars_temp_forget.meanTvP1(Patm)

        ystore[i+1,0]=Patm										# Patm
        ystore[i+1,1]=(isotopes_sim[0][0]/(3676.5e-6)-1.)*1000.	# d15N
        ystore[i+1,2]=isotopes_sim[1][0]						# 22Ne/20Ne
        ystore[i+1,3]=isotopes_sim[2][0]						# 38Ar/36Ar
        ystore[i+1,4]=isotopes_sim[3][4]						# 86Kr/84Kr
        ystore[i+1,5]=isotopes_sim[4][7]						# 136Xe/130Xe
        ystore[i+1,6]=tglobal									# tglobal
        ystore[i+1,7]=Hscale									# Hscale
        ystore[i+1,8]=Matm[0]	# co2
        ystore[i+1,9]=Matm[1]	# n2
        ystore[i+1,10]=Matm[2]  # h2o
        
        if ((t[i]-last_output)>=output_interval):
        	if printOutput:
        		print('time: ' + str(t[i]/1e9) + '; impacts this step: ' + str(num_impacts) + \
					 '; Patm: ' + str(Patm) + \
					 '; d15N: ' + str(ystore[i+1,2]) + '; euv: ' + str(euv_flux) + \
					 '; tg: ' + str(tglobal) )
        	last_output = t[i]
        if printOutput:
        	print('Final time: ' + str(t[-1]/1e9))
    """
        ----------------------------------------------------------------------------------
    """
    
    if output_file is not None:
    	np.savez(output_file,t=t,Patm=ystore[:,0], \
    		isoRat=ystore[:,1:6],\
    		element_mass=isotopic_data.element_mass, \
    		abundances_wrt_sol=mole_elements/isotopic_data.solar_abundances, \
    		isotopes=[isotopic_data.isotopes[0],isotopic_data.isotopes[2],\
    			isotopic_data.isotopes[7], \
    				isotopic_data.isotopes[15]], \
    				temp=ystore[:,6] )
    	print(output_file)
    	
    return t,ystore,mole_elements,isotopes_sim, Amars
    	

