import numpy as np
import scipy.interpolate as scint 

data= dict()

# figure 3, fox and dalgarno
# proportional to N2 at exobase
data['phi1'] = np.array([[0.,0.026390532,0.053538147,0.098437286,0.158672989,0.251542364,\
    0.367219551,0.447339199,0.631762776,0.756948626,0.984524335,1.36732528,2.20106072,\
    3.845786827,5.989664298,9.961049426], \
    [0.,1.002557545,1.79028133,2.613810742,3.294117647,3.867007673,4.278772379,4.475703325,\
    4.744245524,4.869565217,5.012787724,5.138107417,5.245524297,5.31713555,5.335038363,\
    5.335038363]])
# N2+ + O -> NO+ + N
data['phi2'] = np.array([[0.,0.026390532,0.035480321,0.04321776,0.055303306,0.071942254,\
    0.098309313,0.123722377,0.174619187,0.234596509,0.290325529,0.377358103,0.53223789,\
    0.714688304,0.959722742,1.421980995,2.005943979,3.021925968,4.478588651,7.088242098,\
    10.00419729], \
    [0.,1.002557545,1.271099744,1.432225064,1.62915601,1.84398977,2.058823529,2.148337596,\
    2.148337596,2.112531969,2.040920716,1.89769821,1.611253197,1.360613811,1.127877238,\
    0.859335038,0.644501279,0.465473146,0.304347826,0.179028133,0.179028133]])
# N2+O++ -> N+ + N + O+
data['phi3'] = np.array([[0.,0.026827119,0.036659341,0.046143565,0.063050033,0.086132761,\
    0.117641429,0.166022728,0.205461865,0.271495116,0.347154507,0.429585921,0.549278502,\
    0.990824239,1.352771744,2.071732406,2.923141181,3.738223413,5.104865195,6.746646745,\
    9.064711872,10.00251816], \
    [0.,1.002557545,1.217391304,1.378516624,1.557544757,1.647058824,1.647058824,1.611253197,\
    1.539641944,1.432225064,1.306905371,1.199488491,1.056265985,0.734015345,0.572890026,\
    0.411764706,0.286445013,0.21483376,0.143222506,0.10741688,0.10741688,0.10741688]])
    

phi1=scint.interp1d(data['phi1'][0],data['phi1'][1],fill_value='extrapolate',kind='linear')
phi2=scint.interp1d(data['phi2'][0],data['phi2'][1],fill_value='extrapolate',kind='linear')
phi3=scint.interp1d(data['phi3'][0],data['phi3'][1],fill_value='extrapolate',kind='linear')

amu=1.66e-27
g_mars=3.721
deltaz=40.*1000. # 40000 m
T=100.
amu_co2=44
amu_n2=28
kBoltz = 1.381e-23


def carbon_escape(age_young,age_old,maxage):
    # see equation 8 from Kurokawa et al
    
    # but integrate and divide by timestep
    Fcph = 7.9e23 * \
        (age_old**0.14 - age_young**0.14)/(age_old-age_young) / 0.14 /\
        maxage**-0.86
    return Fcph
    
# The phi curves are normalised to 1 at x = 0.0264, the present-day BULK
# (homopause) N2/CO2 ratio: Fox & Dalgarno's present-day fluxes already come
# from a model atmosphere with realistic diffusive separation above the
# homopause. Multiplying the bulk ratio by R_diff (~17) before the lookup
# (previous behaviour) counted that enrichment twice: at present-day
# composition it gave phi1 = 4.4, i.e. 4x Fox & Dalgarno's own present flux.
# Hu & Thomas (2023, Eq. 15) likewise scale photochemical N escape with the
# bulk N2 mixing ratio. Set EXOBASE_RATIO_LOOKUP = True to restore the old
# behaviour for comparison.
EXOBASE_RATIO_LOOKUP = False

# Fox & Dalgarno (1983) Table 3 escape fluxes (cm^-2 s^-1) by process:
# phi1 group = photodissociation, photodissociative ionisation, electron
# impact dissociative ionisation, dissociative recombination of N2+;
# phi2 = N2+ + O -> NO+ + N; phi3 = N2 + O++ -> N+ + N + O+.
# Totals: low 2.3e5, high 8.9e5, average 5.6e5. Jakosky et al. (1994) use
# Fox (1993a)'s average-activity total of 4.8e5, close to 'average'.
# Default 'low': with it the photochemical N lifetimes reproduce Kurokawa et
# al. (2018) Fig. 14a (~2 Gyr at 6 mbar today, ~0.55 Gyr at 2 Ga, ~0.15 Gyr at
# 3.5 Ga) within reading accuracy, whereas 'average' gives lifetimes 2-3x
# shorter. Fox & Dalgarno (1980) also predicted a present rate of 2.0e5.
FD83_FLUXES = {
    'low':     ([5.9e4, 1.7e4, 1.4e4, 8.4e4], 4.7e4, 4.8e3),
    'average': ([1.0e5, 3.0e4, 3.1e4, 2.9e5], 8.9e4, 9.4e3),
    'high':    ([1.4e5, 6.0e4, 4.8e4, 5.0e5], 1.3e5, 1.4e4),
}
FD83_SOLAR_ACTIVITY = 'low'


def nitrogen_escape(N_CO2,N_N2,euv,euv_max,Amars,Hatm):
    xval = N_N2 / N_CO2
    if EXOBASE_RATIO_LOOKUP:
        from MarsSputtering_Calc6_7 import upper_atmosphere
        z, upper_T = upper_atmosphere(euv)
        xval = xval * np.exp((amu_co2-amu_n2)*amu*g_mars*z/upper_T/kBoltz)
    phi1_val = phi1(xval)
    phi2_val = phi2(xval)
    phi3_val = phi3(xval)
    
    # from table 3, fox and dalgarno we can see the last two rates are for phi2 and phi3
    # whereas the first 4 are proportional to N2 I think?, so add up the first four 
    # and multiply by phi, this is proportional to euv
    f1, f2, f3 = FD83_FLUXES[FD83_SOLAR_ACTIVITY]
    flux1=sum(f1)*phi1_val*euv/euv_max * Amars * 1e4 
    
    # phi2
    flux2=f2*phi2_val*(euv/euv_max)**1.5 * Amars * 1e4 
    # phi3
    flux3=f3*phi3_val*(euv/euv_max)**1.5 * Amars * 1e4

    return (flux1,flux2,flux3)
    
    