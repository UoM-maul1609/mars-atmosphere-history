"""
digitised data from figure 3, Kurokawa et al (2018)
"""

import numpy as np
import scipy.interpolate as scint 

data = np.array([1.01690620969228, -0.012470382840753302, \
2.0245992260917642, 0.026499563536603432, \
3.0322922424912484, 0.4941389200648487, \
4.039985258890733, 1.5463274722533988, \
4.975700202690254, 4.196283825913454, \
5.983393219089738, 9.301346801346803, \
7.063064308089183, 15.185808704327222, \
7.926801179288743, 17.95267489711934, \
9.006472268288187, 18.186494575383463, \
10.014165284687673, 18.65413393191171, \
11.021858301087159, 18.498254146402296, \
12.029551317486641, 15.497568275346055, \
13.037244333886123, 11.44469385210126, \
13.972959277685645, 9.37928669410151, \
15.052630366685092, 8.989587230327972, \
15.988345310484613, 8.794737498441204, \
17.06801639948406, 8.950617283950619, \
18.003731343283583, 9.37928669410151, \
19.011424359683062, 9.613106372365632, \
20.019117376082548, 9.145467015837387, \
20.954832319882073, 8.13224841002619, \
21.96252533628155, 6.690360394064099, \
22.970218352681037, 5.755081681007608, \
23.977911369080516, 5.0925925925925934, \
24.985604385480002, 4.585983289686997, \
25.993297401879488, 3.9624641476493316, \
26.929012345679013, 3.533794737498443, \
28.008683434678453, 3.105125327347551, \
29.016376451077946, 2.793365756328722, \
30.024069467477425, 2.4426362389325362, \
30.959784411276942, 2.1698466142910604, \
31.967477427676428, 1.9360269360269378, \
32.975170444075914, 1.7022072577628151, \
33.98286346047539, 1.5073575258760457, \
34.99055647687488, 1.3125077939892762, \
35.998249493274365, 1.195597954857213, \
36.93396443707389, 0.9617782765930905, \
38.01363552607334, 0.8838383838383841, \
39.021328542472816, 0.7669285447063245, \
39.95704348627233, 0.6500187055742614, \
40.96473650267182, 0.6110487591969118, \
42.044407591671266, 0.4941389200648487, \
42.90814446287083, 0.1044394562913098, \
43.91583747927031, -0.012470382840753302, \
44.995508568269756, -0.012470382840753302, \
45.93122351206928, -0.012470382840753302, \
46.93891652846876, -0.012470382840753302, \
47.94660954486825, -0.012470382840753302, \
48.88232448866777, -0.012470382840753302])

datx = data[0:-1:2]*1e3
daty = data[1::2]

daty=np.maximum(daty,0.)

# The digitised ordinates represent dN/dv, not probability per data point.
# Integrate over actual velocity spacing to obtain the cumulative distribution.
from scipy.integrate import cumulative_trapezoid
_datyy = cumulative_trapezoid(daty, datx, initial=0.0)
_datyy /= _datyy[-1]
# The zero-density tail gives repeated CDF values, invalid for interp1d.
datyy, _unique_indices = np.unique(_datyy, return_index=True)
_cdf_velocities = datx[_unique_indices]
# Explicit support at zero probability (first positive-density bin starts later).
if datyy[0] != 0.0:
    datyy = np.insert(datyy, 0, 0.0)
    _cdf_velocities = np.insert(_cdf_velocities, 0, datx[0])

"""
    return the velocities
"""
def inv_cdf1(rand_numbers=np.random.rand((1000000))):
    inv_cdf=scint.interp1d(datyy,_cdf_velocities,bounds_error=False,fill_value=(_cdf_velocities[0], _cdf_velocities[-1]),kind='linear')

    return inv_cdf(rand_numbers)
    
    
    
    
if __name__ == "__main__":
    """
        sample and plot
    """
    import matplotlib.pyplot as plt

    inv_cdf=scint.interp1d(datyy,_cdf_velocities,fill_value='extrapolate',kind='linear')
    
    
    rand_numbers=np.random.rand((1000000))
    
    vels = inv_cdf(rand_numbers)
        
     
    bin_edges=np.linspace(0,50000,25)
    (n,bin_edges)=np.histogram(vels,density=True,range=(0,50000),bins=bin_edges)
    plt.ion()
    plt.plot((bin_edges[0:-1]+bin_edges[1:])*0.5,n,'.-')
    plt.xlabel('vels (m s$^{-1}$)')
    plt.ylabel('Histogram - normalized')
    plt.legend(['PDF','Sampled using Monte-Carlo'])


    