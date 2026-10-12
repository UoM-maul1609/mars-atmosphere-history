import numpy as np

# Solar ages (Gyr). Kurokawa §2.3: the EUV saturation phase covers the first
# 0.1 Gyr of solar evolution, "the first 0.04 billion years in our
# calculation", so the model start (4.5 Ga) is at a solar age of 0.06 Gyr.
SOLAR_AGE_AT_START = 0.06
PRESENT_SOLAR_AGE = 4.56
SATURATION_AGE = 0.1


# CO2 sputtering rates (s^-1) at 6x, 3x and 1x present EUV.
#  'luhmann_published': Luhmann et al. (1992) Table 1.
#  'jakosky_corrected': Jakosky et al. (1994) Table I, footnote c: Luhmann's
#      rates "corrected" after R. E. Johnson and J. G. Luhmann found a
#      numerical error (net CO2 loss 2.3e42 rather than 6e42 molecules).
# Kurokawa et al. (2018) Fig. 14 lifetimes match the corrected values: N
# and Ar sputtering lifetimes are respectively 3x and 1.5x the published-rate
# values at present, and ~6x and ~3x at 3.5 Ga (with the N2 yield halving of
# Jakosky Eq. 12 accounting for the extra 2x on N).
CO2_SPUTTERING_RATES = {'luhmann_published': [3.e26, 6.e25, 3.e23],
                        'jakosky_corrected': [1.e26, 3.e25, 2.e23]}
CO2_SPUTTERING_SOURCE = 'jakosky_corrected'


def euv_relative(timeGyr):
    """EUV flux relative to present (Ribas et al. 2005, F ~ tau^-1.23, with
    saturation before 0.1 Gyr), as stated in Kurokawa §2.3.

    timeGyr is time since the model start (4.5 Ga). The previous version
    extrapolated Luhmann's 1/3/6x points log-linearly in time, giving 12x at
    4.5 Ga and 9x at 4.1 Ga, against the paper's ~100x and ~17x.
    """
    tau = np.maximum(np.asarray(timeGyr, dtype=float) + SOLAR_AGE_AT_START,
                     SATURATION_AGE)
    return np.maximum((tau / PRESENT_SOLAR_AGE) ** -1.23, 1.0)


def time_of_euv(level):
    """Time since the model start (Gyr) at which euv_relative == level."""
    return PRESENT_SOLAR_AGE * level ** (-1. / 1.23) - SOLAR_AGE_AT_START


def sputtering_co2_rate(timeGyr):
    """
        CO2 sputtering rate F_CO2,sp^tot (s^-1), Luhmann et al. (1992) Table 1
        for 6, 3 and 1x EUV, placed where the Ribas EUV law reaches those
        levels (~0.99, 1.81 and 4.5 Gyr after the model start) and
        interpolated/extrapolated in time (Kurokawa §2.3).

        Interpolation is linear in log10(F): the three values span three
        orders of magnitude. Linear interpolation in F (previous version) gave
        F_CO2 = 2.4e25 s^-1 at 1 Ga, i.e. 80x present for an EUV of only ~1.4x;
        log interpolation gives 2.5e24, consistent with Luhmann's steep
        F-EUV dependence (F ~ EUV^4.8 between 1x and 3x).
        Returns (relative EUV, F_CO2).
    """
    # Luhmann's rates are for 6, 3 and 1x present EUV. Place them at the
    # times the Ribas EUV law reaches those levels (3.49, 2.69 and 0 Ga),
    # as the kinks in Kurokawa Fig. 14b (at ~3.5 and ~2.7 Ga) show they did.
    time_grid = np.array([time_of_euv(6.), time_of_euv(3.), 4.5])
    log_f = np.log10(np.array(CO2_SPUTTERING_RATES[CO2_SPUTTERING_SOURCE]))
    t = np.asarray(timeGyr, dtype=float)
    slope_early = (log_f[1] - log_f[0]) / (time_grid[1] - time_grid[0])
    log_val = np.where(t < time_grid[0],
                       log_f[0] + slope_early * (t - time_grid[0]),
                       np.interp(t, time_grid, log_f))
    val2 = np.maximum(10 ** log_f[-1], 10 ** log_val)
    return (euv_relative(t), val2)


if __name__ == "__main__":
    # test the function by evaluating and plotting
    import matplotlib.pyplot as plt

    t_test_gyr = np.linspace(0., 4.5, 200)
    (euv_flux, f_co2_flux) = sputtering_co2_rate(t_test_gyr)

    plt.ion()
    fig = plt.figure()
    ax1 = plt.subplot(211)
    plt.plot(t_test_gyr, f_co2_flux)
    plt.yscale('log')
    plt.xlabel('time since 4.5 Ga (Gyr)')
    plt.ylabel('F$_{CO2}$')
    plt.text(0.2, 0.8, '(a) ', transform=ax1.transAxes)
    ax2 = plt.subplot(212)
    plt.plot(t_test_gyr, euv_flux)
    plt.yscale('log')
    plt.text(0.2, 0.8, '(b) ', transform=ax2.transAxes)
    plt.xlabel('time since 4.5 Ga (Gyr)')
    plt.ylabel('Relative EUV flux')
    plt.show()
    fig.tight_layout()
