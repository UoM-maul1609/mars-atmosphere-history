import numpy as np
def svp_ice01(t):
    """
    calculate the saturation vapour pressure, PH2O, according to Murphy and Koop (2005)
    see also Kurokawa et al. 2018, equation 1
    
    run -i svp_ice.py
    
    """
    pH2O = np.exp(9.550426-5723.265/t + 3.53068*np.log(t)-0.00728332*t)
    
    return pH2O
    
    
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    # create an array of t-values
    T=np.linspace(273.15,223.15,1000)
    pH2O = svp_ice01(T)
    
    plt.ion()
    plt.plot(T-273.15,pH2O)
    plt.xlabel(r'T ($^\circ$C)')
    plt.ylabel('Saturation vapour pressure over ice surface (Pa)')
    plt.show()
    
    
    