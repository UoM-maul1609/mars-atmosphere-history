import numpy as np
import isotopic_data

# Kurokawa et al. (2018) Table 3 (Flynn 1997): ELEMENT mass deposition rates,
# g Gyr^-1.
IDP_table3 = dict()
IDP_table3['Ne'] = 4.89e12
IDP_table3['Ar'] = 1.11e12
IDP_table3['Kr'] = 1.64e10
IDP_table3['Xe'] = 7.33e9

# Isotope group index in isotopic_data (0=N, 1=Ne, 2=Ar, 3=Kr, 4=Xe) and the
# reference isotope mass of each group.
_GROUP = {'Ne': (1, 20), 'Ar': (2, 36), 'Kr': (3, 84), 'Xe': (4, 130)}
_GROUP_SLICES = {1: slice(1, 2), 2: slice(2, 3), 3: slice(3, 8), 4: slice(8, 16)}


def _reference_fraction_and_mean_mass(element):
    """Reference-isotope atom fraction and mean atomic mass of IDP gas.

    Uses the IDP isotopic composition in Kurokawa Table 4 (21Ne and 40Ar are
    not modelled, as in the paper).
    """
    group, ref_mass = _GROUP[element]
    sl = _GROUP_SLICES[group]
    ratios = np.array(isotopic_data.isotope_comps['IDPs'][sl], dtype=float)
    masses = np.array(isotopic_data.isotope_mass[sl], dtype=float)
    weights = np.concatenate(([1.0], ratios))
    weights = weights / weights.sum()
    mean_mass = np.dot(weights, np.concatenate(([ref_mass], masses)))
    return weights[0], mean_mass


def get_mole(elements, Ne=1):
    """IDP deposition in REFERENCE-ISOTOPE moles per Gyr (20Ne, 36Ar, 84Kr,
    130Xe), matching the basis of the noble-gas element ledger.

    Previously the element mass was divided by the reference-isotope mass and
    added in full, i.e. total-element moles were added to a reference-isotope
    inventory. That overstated the IDP supply by 1/f_ref: about 23x for Xe,
    1.8x for Kr, 1.2x for Ar and 1.07x for Ne.
    """
    rates = np.zeros(len(elements))
    for i, element in enumerate(elements):
        f_ref, mean_mass = _reference_fraction_and_mean_mass(element)
        rates[i] = IDP_table3[element] / mean_mass * f_ref
        if element == 'Ne':
            rates[i] *= Ne
    return rates
