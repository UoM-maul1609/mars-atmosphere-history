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
        X_gas=0.01,f_comet=0.001, C_vol=5.0, crater_model=1,C_Ne_IDP = 10.,\
	dynamo_time = 4.1, output_file=None):
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
    tmars = 210. # temperature of martian surface - may change???
    Rgas = 8.314 # ideal gas constant
    #MolW_atm = 44.01e-3 # molecular weight of martian atmosphere - may change???
    MolW_atm = np.array([44.01e-3, 28.0134e-3, 18.01528e-3])
    
    Obliquity=20. # the obliquity of mars
    
    rho_pr=2600. # density of impactor
    Pcollapse = mars_temp_forget.meanPvObliquity(Obliquity)
    sampleFlag = 0 # 0=once, 1=every time-step; 
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
    Rdiffs_isotopes = MarsSputtering_Calc6_7.calculate_rdiffs(element_mass,isotope_mass)
    Rdiffs_isotopes[0][0]=1.43
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
    ph2o = svp_ice.svp_ice01(tmars)
    Ninit=Natm[1]
    Hscale=Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars) # scale height - m - may change???
    t = np.mgrid[tinit:tfinal+tstep:tstep]

    (euv_max, f_co2_flux_max) = sputtering_co2.sputtering_co2_rate((t[0]-tinit)/1e9)
    nsteps = len(t)
    n_ode=11
    rng=np.random.RandomState(runNo)
    rand_numbers=rng.rand((100000))
    rand_numbers2=rng.rand((100000))
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
    
    # sample the sizes of the impactors
    diams, mass = impactor_model.sample_sizes(rand_numbers, rho_pr)
    # sample the impact velocities - not sure if this should be a different random number
    vels = impactor_velocity.inv_cdf1(rand_numbers2)
    # total number of impacts over the whole history
    num_impactor1_tot=crater_chronology.N1_between(Amars,0.,4.5)
    num_impactor20_tot=crater_chronology.N20_between(Amars,0.,4.5,model=crater_model)
    
    # escape velocity
    uesc=impactor_atmos_loss_gain_eq4_5.escape_velocity(Ggrav,Mmars,Rmars)
    
    """
        ----------------------------------------------------------------------------------
    """
    if printOutput:
    	print(total_impactor_mass/np.sum(mass))
    	print('Note, the pressure for atmospheric collapse is ' + str(Pcollapse))
    #nsteps=0

    """
        main loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    for i in range(nsteps-1):
        if sampleFlag == 1:
            rand_numbers=rng.rand((100000))
            rand_numbers2=rng.rand((100000))
            # sample the sizes of the impactors
            diams, mass = impactor_model.sample_sizes(rand_numbers, rho_pr)
            # sample the impact velocities - not sure if this should be a different random number
            vels = impactor_velocity.inv_cdf1(rand_numbers2)
        
        
        
        # 1. collapse? 0.5 bar and assumed to be 6 mbar
        
        # 2. calculate H2O - two temperatures
        
        # 3. impacts of asteroids and comets. during this time-step
        num_impactor1=crater_chronology.N1_between(Amars,-t[i+1]/1e9,-t[i]/1e9)
        num_impactor20=crater_chronology.N20_between(Amars,-t[i+1]/1e9,-t[i]/1e9,model=crater_model)
        

        # note that results at each time get scaled by:
#         scaling = total_impactor_mass / np.sum(mass) * num_impactor1 / num_impactor1_tot
        scaling = total_impactor_mass / np.sum(mass) * num_impactor20 / num_impactor20_tot
        
        
        
        # 4. loss of the atmosphere following empirical formulae
        # sum(Matm) / sum(Natm) is the molecular weight
        rho_0 = Patm*1.e5/tmars/(Rgas/(np.sum(Matm)/np.sum(Natm)))
        xi=np.maximum(impactor_atmos_loss_gain_eq4_5. \
            dimensionless_xi(diams,rho_t,rho_pr, vels, uesc, Hscale, rho_0)  , 1.e-10) 
        deltaM=impactor_atmos_loss_gain_eq4_5.change_atmos(xi,vels,uesc,mass)
#         deltaM=impactor_atmos_loss_gain_eq4_5.fractional_change_atmos(xi,vels,uesc)
#         deltaM = deltaM*np.sum(Natm*MolW_atm)
        
        # left in the projectiles is (1-Xpr)*mass, but then only a certain % of this 
        Xpr=np.maximum(impactor_atmos_loss_gain_eq4_5. \
            normalised_projectile_mass(rho_t,rho_pr,vels,uesc,xi),0)
        # deltaM is the atmospheric loss for impactor = mass
        # Xpr*mass*0.01 is the mass added from the projectiles
        mass_added = np.sum((1.-Xpr)*mass)*X_gas*scaling
        if(t[i]  > -(dynamo_time)*1e9):
            f_comet1 = f_comet
        else:
            f_comet1 = 0.
            
        mass_added_comet = mass_added * f_comet1
        mass_added_asteroid = mass_added * (1.-f_comet1)
        
        
        
        #deltaM=deltaM-Xpr*mass*0.01
        dM = np.sum(deltaM)*scaling
#         dM = np.sum(Matm)*(1.-dM1)
        const1 = np.minimum(dM / np.sum(Natm*MolW_atm),0.5)
        Natm = Natm *(1.- const1)
        # remove elements in the impact too
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
        F_i_sp1=MarsSputtering_Calc6_7.F_i_sp(f_co2_flux, \
            mole_elements,Natm[0],['C','N','Ne','Ar','Kr','Xe'])
        # CO2 and N2 - only after 4.1 Gyr
        if((-t[i]/1e9< dynamo_time ) and sputtering_flag):
            Natm[0] = Natm[0] - F_i_sp1[0] * tstep *86400*365/6.02e23
            Natm[1] = Natm[1] - 0.5*F_i_sp1[1] * tstep *86400*365/6.02e23
            #Natm = np.maximum(Natm,0.)
            # now elements
            N = mole_elements[1:].copy().tolist()
            dN = [F_i_sp1[1] * tstep *86400*365/6.02e23, \
                F_i_sp1[2] * tstep *86400*365/6.02e23, \
                F_i_sp1[3] * tstep *86400*365/6.02e23, \
                F_i_sp1[4] * tstep *86400*365/6.02e23, \
                F_i_sp1[5] * tstep *86400*365/6.02e23]
            mole_elements[0] = mole_elements[0] - F_i_sp1[0] * tstep *86400*365/6.02e23
            mole_elements[1] = mole_elements[1] - dN[0]
            mole_elements[2] = mole_elements[2] - dN[1]
            mole_elements[3] = mole_elements[3] - dN[2]
            mole_elements[4] = mole_elements[4] - dN[3]
            mole_elements[5] = mole_elements[5] - dN[4]
        else: 
            N = mole_elements[1:].copy().tolist()
            dN =[0]*len(N)  
        
        
        # 5b photochemical escape
        # carbon - rate of escape
        if ((-t[i]/1e9< dynamo_time ) and pce_flag):
                Fcph = photo_chemical.carbon_escape((t[i]-tinit)/1e9,\
                    (t[i+1]-tinit)/1e9,(-tinit)/1e9)  # C atoms/s
                # nitrogen - rate of escape
                (fn2_1,fn2_2,fn2_3)=photo_chemical.nitrogen_escape(\
                    np.maximum(Natm[0],1e-3),Natm[1],np.minimum(euv_flux,6),1.,Amars, Hscale)
                # total loss rate
                fn2 = (fn2_1+fn2_2+fn2_3)*Natm[1]/Ninit
                # it is N2
                dN2 =[fn2*tstep*86400*365/6.02e23 ,0., 0.,0.,0.]
        else:
                Fcph=0.
                # Never reuse the cometary dN2 array as an escape flux.
                dN2 = [0.] * 5

        Natm[0] = Natm[0] - Fcph*tstep*86400*365/6.02e23
        Natm[1] = Natm[1] - dN2[0] # N2
        mole_elements[0] = mole_elements[0] - Fcph*tstep*86400*365/6.02e23
        mole_elements[1] = mole_elements[1] - dN2[0] * 2 # N
        Matm = Natm * MolW_atm

        isotopes_sim = isotopic_data.escape_fr(isotopes_sim, N, dN, dN2, Rdiffs_isotopes)
        
        
        
        # 5c IDPs - just table 3, Kurokawa et al. and fractionation
        N = mole_elements[2:].copy()
        dN = IDP_rate_moles*tstep / 1e9
        mole_elements[2:] = mole_elements[2:] + np.array(dN)
        
        isotopes_sim[1:] = isotopic_data.continuous_sources(isotopes_sim[1:], \
            isotope_sources['IDPs'][1:], N, dN)
        
        
        
        
        # 6. Volcanic degassing - digitise rates and incorporate total
        Matm = Natm * MolW_atm
        h2o=volcano_outgassing.get_outgas_rate(t[i]/1e9,'h2o') * C_vol
        co2=volcano_outgassing.get_outgas_rate(t[i]/1e9,'co2') * C_vol
        n2 =volcano_outgassing.get_outgas_rate(t[i]/1e9,'n2')  * C_vol
        Matm = Matm + np.array([co2,n2,h2o])*MolW_atm*tstep*86400*365/6.02e23
        Matm[2] = np.minimum(Matm[2]*grav_mars/Amars,ph2o)*Amars/grav_mars
        Natm = Matm / MolW_atm
        
        
        
        
        # readjust Matm: firstly it cant go above the svp from Forget et al:
        if(Matm[0]*grav_mars/Amars/1.e5 > mars_temp_forget.data['collapse2']):
        	Matm[0]=mars_temp_forget.data['collapse2']*Amars*1e5/grav_mars
        	Patm = np.sum(Matm)*grav_mars/(Amars)/1e5
        	mole_elements[0] = Matm[0] / MolW_atm[0]

        # readjust Matm: secondly it cant go below 6e-3 (observed):
        Patm = np.maximum(np.sum(Matm)*grav_mars/Amars/1.e5, 6e-3)        
        # if Patm < 0.5 bar - the CO2 condenses at the poles - only the CO2 so will have 
        # to wait to do properly. e.g. need to know how lack of absorption by CO2 will affect 
        # temperature, and then link this to the clausius clapeyron equation
        # same for water vapour really, perhaps having a reservoir
		# readjust Matm: thirdly it cant go below the threshold where it collapses due to
		# radiative effect forget et al:
        if obliquity_flag==1:
            Pcollapse=mars_temperature.co2_vap_press( \
               mars_temp_forget.minTvP_bradley_0(Patm))/1e5
        elif obliquity_flag==2:
            Pcollapse=mars_temperature.co2_vap_press( \
               mars_temp_forget.minTvP_bradley_45(Patm))/1e5
        elif obliquity_flag==3:
            Pcollapse=mars_temperature.co2_vap_press( \
		       mars_temp_forget.minTvP_bradley_90(Patm))/1e5
        elif obliquity_flag==-1:
            Pcollapse=0.5
            
        if((Patm > Pcollapse) and (obliquity_flag != -1)):
            Matm[0] = Pcollapse*1e5*Amars / grav_mars # mass of CO2
            Patm = np.sum(Matm)*grav_mars/(Amars)/1e5
            #Matm[1:] = 0. # if it's collapsed, just assume it is CO2
            mole_elements[0] = Matm[0] / MolW_atm[0]
        if((Patm < Pcollapse) and (obliquity_flag == -1)):
            Patm = 6.e-3
            Matm[0] = Patm*1e5*Amars / grav_mars # mass of atmosphere
            Patm = np.sum(Matm)*grav_mars/(Amars)/1e5
            #Matm[1:] = 0. # if it's collapsed, just assume it is CO2
            mole_elements[0] = Matm[0] / MolW_atm[0]

        
        
        Natm = Matm / MolW_atm
        N= mole_elements[1:].copy()
        dN = [2.*n2*tstep*86400*365/6.02e23] # this is moles of N
        dN.extend(co2*tstep*86400*365/6.02e23* \
            isotopic_data.abundances1['Volcanic degassing'][1][2:] / \
            isotopic_data.abundances1['Volcanic degassing'][1][0] * \
            isotopic_data.solar_abundances[2:] / isotopic_data.solar_abundances[0])
            
        # Atmospheric CO2 was already updated with volcanic CO2 before the
        # collapse adjustment; adding it again here double-counts carbon.
        # Keep the element inventory synchronized with atmospheric CO2.
        mole_elements[0] = Natm[0]
        mole_elements[1] = mole_elements[1] + dN[0] # moles of N - i.e the element
        mole_elements[2:] = mole_elements[2:] + np.array(dN[1:])
        # 7. Isotope fractionation. need to add elemnents also
        isotopes_sim = isotopic_data.continuous_sources(isotopes_sim, \
            isotope_sources['Volcanic degassing'], N, dN)
        
        # update scale height
        Hscale=Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars) # scale height - m - may change???
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
        		print('time: ' + str(t[i]/1e9) + '; num impacts 1: ' + str(num_impactor1) + \
					'; num impacts 20: ' + str(num_impactor20) + '; ratio: ' + \
					 str(num_impactor1/num_impactor20) + '; Patm: ' + str(Patm) + \
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
    	

