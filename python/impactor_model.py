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
