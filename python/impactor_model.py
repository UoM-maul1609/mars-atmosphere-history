import numpy as np
import matplotlib.pyplot as plt


dmin=3e3 # the smallest size of an impactor

# start of function
def impactor_sd(D):
    """ 
        Kurokawa et al. (2018)
         Small (D<300 km) impactors have q=3.5
         Large (D>300 km) impactors have q=2
    """
    
    dNdD = np.zeros((len(D)))
    q = np.zeros((len(D)))
    a = np.ones((len(D)))
    
    ind,=np.where(D<300e3)
    q[ind] = 3.5
    ind,=np.where(D>=300e3)
    q[ind] = 2.
    
    
    # Curve matching for the impactor size distribution
    # 300e3**-3.5 = a * 300e3**-2
    # a=(300e3**-3.5)/(300e3**-2)
    # a = 300e3**(-3.5+2=-1.5)
    a[ind] = 300.e3**(-1.5)
        
    
    dNdD = a* D**-q
    
    return (dNdD,a)
    # end of function

def impactor_csd(D,normalisingFactor,a):
    """
        integrate SD from dmin to D 
    """
    if(D<300e3):
        integral1=((D**-2.5)/-2.5 - (dmin**-2.5)/-2.5) 
        integral1 *= normalisingFactor
    else:
        integral1=((300e3**-2.5)/-2.5 - (dmin**-2.5)/-2.5)
        integral2=a*(1./300e3-1./D) 
        integral1 += integral2
        integral1 *= normalisingFactor
        
        
    return integral1
        

def inv_impactor_csd(prob, normalisingFactor,a):
    """
        compute the inverse of impactor_csd
    """
    # first calculate the threshold CP where the form of the equation changes
    K=-(dmin**-2.5)/2.5
    integral1=(300e3**-2.5)/-2.5 - K
    integral1a = integral1 * normalisingFactor
    if(prob <= integral1a):
        # just find the inverse of the D^-2.5 relation
        D=((prob/normalisingFactor+K )*-2.5)**(-1./2.5)
    else:
        # find the inverse of the second part
        D=1./300e3-(prob/normalisingFactor-integral1)/a
        D=1/D
        
    #D=np.minimum(D,1e7)
    return D


def sample_sizes(rand_numbers=np.random.rand((1000000)),rho_pr=2600):

    sampleSize = len(rand_numbers)
    # evaluate the integral of impactor_sd between dmin and infinity
    integral1=(300e3**-2.5)/-2.5 - (dmin**-2.5)/-2.5

    # evaluate the integral of impactor_sd between 1e-6 and inf
    integral2=300.e3**(-1.5)/300e3

    # total area under curve
    integralt=integral1+integral2
    # scale original distribution by area under curve
    nF=1./integralt



    """
        now, test monte carlo generation
        by sampling cumulative distribution to generate sizes of impactors
    """
    diams = np.zeros((sampleSize))
    for i in range(len(rand_numbers)):
        diams[i] = inv_impactor_csd(rand_numbers[i],nF,300.e3**(-1.5)) 
    diams = np.minimum(diams,2000.e3)
    mass = np.pi/6*rho_pr*diams**3  # diameter [m], density [kg/m^3], mass [kg]
    return (diams,mass)

if __name__ == "__main__":

    # starts here
    D=np.logspace(np.log10(3e3),6,10000)
    (imp,a)=impactor_sd(D)

    # evaluate the integral of impactor_sd between dmin and infinity
    integral1=(300e3**-2.5)/-2.5 - (dmin**-2.5)/-2.5

    # evaluate the integral of impactor_sd between 1e-6 and inf
    integral2=a[-1]/300e3

    # total area under curve
    integralt=integral1+integral2
    # scale original distribution by area under curve
    nF=1./integralt
    imp = imp *nF # PDF

    # the cumulative distribution function
    Prob=np.zeros((len(D)))
    for i in range(len(D)):
        Prob[i]=impactor_csd(D[i],nF,a[-1])

    plt.ion()
    plt.show()
    
    fig=plt.figure()
    # plot the PSD - Probability Size Distribution / PDF - area under curve is 1.
    ax1=plt.subplot(211)
    plt.plot(D,imp)
    plt.yscale('log')
    plt.xscale('log')
    x=plt.xlim()
    y=plt.ylim()
    plt.ylabel('PDF from dmin - 10$^6$')
    ax1.text(0.1,0.2,'(a) histogram',transform=ax1.transAxes)

    # plot the CSD - cumulative Size Distribution / CDF - varies between 0 and 1
    ax2=plt.subplot(212)
    plt.plot(D,Prob)
    plt.xscale('log')
    plt.ylabel('Cumulative Frequency')
    plt.xlabel('D (m)')
    ax2.text(0.5,0.2,'(b) cumulative histogram',transform=ax2.transAxes)


    """
        now, test monte carlo generation
        by sampling cumulative distribution to generate sizes of impactors
    """
    rand_numbers=np.random.rand((1000000))
    diams = np.zeros((1000000))
    for i in range(len(rand_numbers)):
        diams[i] = inv_impactor_csd(rand_numbers[i],nF,a[-1]) 

    plt.axes(ax1)
    bin_edges=np.logspace(np.log10(dmin),6,25)
    (n,bin_edges)=np.histogram(diams,density=True,range=x,bins=bin_edges)
    plt.plot((bin_edges[0:-1]+bin_edges[1:])*0.5,n,'.-')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(x)
    plt.ylim(y)
    plt.xlabel('D (m)')
    plt.ylabel('Histogram - normalized')
    plt.legend(['PDF','Sampled using Monte-Carlo'])
    fig.tight_layout()


# ---------------------------------------------------------------------------
# Discrete-impact sampling (October 2026)
# ---------------------------------------------------------------------------
D_BREAK = 300.e3   # m, break between q = 3.5 and q = 2 (Kurokawa §2.2)


def sample_sizes_truncated(rand_numbers, rho_pr=2600., dmax=2000.e3):
    """Vectorised inverse-CDF sampling of the broken power law, truncated at
    dmax and renormalised.

    Unlike sample_sizes, there is no pile-up: sample_sizes integrates the
    q = 2 branch to infinity and then clips every draw above 2000 km to
    exactly 2000 km (a 1.1e22 kg body, 5x Kurokawa's whole 2e21 kg budget).
    """
    u = np.asarray(rand_numbers, dtype=float)
    a = D_BREAK ** -1.5                      # continuity at the break
    d1 = min(dmax, D_BREAK)                  # top of the q = 3.5 branch
    i1 = (dmin ** -2.5 - d1 ** -2.5) / 2.5
    i2 = a * (1. / D_BREAK - 1. / dmax) if dmax > D_BREAK else 0.
    x = u * (i1 + i2)
    small = x <= i1
    d = np.empty_like(x)
    d[small] = (dmin ** -2.5 - 2.5 * x[small]) ** (-1. / 2.5)
    d[~small] = 1. / (1. / D_BREAK - (x[~small] - i1) / a)
    d = np.clip(d, dmin, dmax)
    return d, np.pi / 6. * rho_pr * d ** 3


def budget_diameter(total_mass=2.e21, rho_pr=2600.):
    """Diameter of a single body whose mass equals the whole budget."""
    return (6. * total_mass / (np.pi * rho_pr)) ** (1. / 3.)


def mean_impactor_mass(rho_pr=2600., dmax=2000.e3):
    """Analytic mean mass of the truncated broken power law."""
    a = D_BREAK ** -1.5
    d1 = min(dmax, D_BREAK)
    i1 = (dmin ** -2.5 - d1 ** -2.5) / 2.5
    # integral of D^3 * dN/dD over each branch
    m1 = (d1 ** 0.5 - dmin ** 0.5) / 0.5               # D^3 * D^-3.5
    if dmax > D_BREAK:
        i2 = a * (1. / D_BREAK - 1. / dmax)
        m2 = a * (dmax ** 2 - D_BREAK ** 2) / 2.        # D^3 * a D^-2
    else:
        i2 = m2 = 0.
    return np.pi / 6. * rho_pr * (m1 + m2) / (i1 + i2)


def number_of_impacts(total_mass=2.e21, rho_pr=2600., dmax=2000.e3):
    """Kurokawa §2.2: total number of impacts from the assumed total mass
    and the size distribution, N = M_total / <m>."""
    return int(round(total_mass / mean_impactor_mass(rho_pr, dmax)))


def sample_impact_history(rng, step_weights, total_mass=2.e21, rho_pr=2600.,
                          dmax=2000.e3):
    """Draw one bombardment history.

    The NUMBER of impacts is fixed at N = total_mass / <m> (Kurokawa §2.2),
    so the total mass is 2e21 kg on average but varies between histories,
    being dominated by the few largest bodies. Each impact gets a size from
    the truncated distribution and a timestep drawn with probability
    proportional to step_weights (the crater chronology).
    Returns (step_index, diameters, masses), sorted by step.
    """
    n = number_of_impacts(total_mass, rho_pr, dmax)
    diams, masses = sample_sizes_truncated(rng.rand(n), rho_pr, dmax)
    w = np.asarray(step_weights, dtype=float)
    cdf = np.cumsum(w) / np.sum(w)
    steps = np.minimum(np.searchsorted(cdf, rng.rand(n), side='right'),
                       len(w) - 1)
    order = np.argsort(steps, kind='stable')
    return steps[order], diams[order], masses[order]
