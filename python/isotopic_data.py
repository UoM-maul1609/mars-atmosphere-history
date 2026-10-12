import numpy as np

# solar abundances of these isotopes, table 6, Lodder, 20xx
solar_abundances=np.array([7.001e6, 1.943e6, 1.996e6, 86710, 31.38, 0.236])
carbon_12_pc = 98.8922
nitrogen_14_pc = 99.6337
solar_percents = np.array([98.8922, 99.6337, 92.9431, 85.5946, 56.90, 4.38])


# figure 2, Kurokawa et al
abundances1 = dict()
volx = np.array([12, 14, 20, 36, 84, 130])
voly = np.array([2.257863844,2.757493323,0.000446188,0.949493424,55.30511164,16.66728846])
abundances1['Comets'] = np.array([volx,voly])

volx = np.array([12,14, 20, 36, 84, 130])
voly = np.array([np.nan, np.nan, 0.005999776,0.002360426,0.057817309,1.729585844])
abundances1['IDPs'] = np.array([volx,voly])

volx = np.array([12, 14, 20, 36, 84, 130])
voly = np.array([0.727338094,0.120333307,1.36E-08,1.03771E-06,4.63015E-05,0.000665505])
abundances1['Asteroids'] = np.array([volx,voly])

voly = np.array([3.0137E-06,4.66E-07,2.97E-11,3.37E-09,3.13E-07,2.56E-07])
abundances1['Mars'] = np.array([volx,voly])

voly = np.array([0.000175539, 8.1881E-06, 4.14E-11,3.88E-11,7.49E-09,1.56E-08])
abundances1['Volcanic degassing'] = np.array([volx,voly])

# digitised data from obs (figure 6 kurokawa)
voly = np.array([0.000002820281104255998, \
			4.579506855064379e-7, \
			2.820281104255998e-11,\
			3.359818286283788e-9,\
			3.0575961084169417e-7,\
			2.6723945955333497e-7])
			
abundances1['Obs'] = np.array([volx,voly])

# Present-day observations, Kurokawa Table 4, as [low, central, high] in the
# plotted quantities: d15N, 20Ne/22Ne, 38Ar/36Ar (= 1/(36Ar/38Ar = 4.2 +/- 0.1)),
# 86Kr/84Kr, 136Xe/130Xe. (Previously d15N was 620 [420-800] instead of
# 572 +/- 82, and 38Ar/36Ar was 0.20 +/- 0.05 instead of 0.238 +/- 0.006.)
errorbars1=[[572-82,572,572+82],[9.4,10.1,10.8],[1/4.3,1/4.2,1/4.1],\
     [0.3006-0.0027,0.3006,0.3006+0.0027],[2.294-0.024,2.294,2.294+0.024]]

atom_wgts=volx/1000.


element_mass=volx
del volx,voly

# table 4, Kurokawa et al.
isotopes = ['d15N',
            '20Ne/22Ne',
            '36Ar/38Ar',
            '78Kr/84Kr',
            '80Kr/84Kr',
            '82Kr/84Kr',
            '83Kr/84Kr',
            '86Kr/84Kr',
            '124Xe/130Xe',
            '126Xe/130Xe',
            '128Xe/130Xe',
            '129Xe/130Xe',
            '131Xe/130Xe',
            '132Xe/130Xe',
            '134Xe/130Xe',
            '136Xe/130Xe']
            
isotopes[1] = '22Ne/20Ne'
isotopes[2] = '38Ar/36Ar'

isotope_mass=[15,22,38,78,80,82,83,86,124,126,128,129,131,132,134,136]

isotope_comps = dict()
isotope_comps['Mars'] = [572,10.1,4.2,6.37e-3,4.09e-2,0.2054,\
                         0.2034,0.3006,2.45e-2,2.12e-2,0.4767, \
                         16.400,5.147,6.460,2.587,2.294]
                         
isotope_comps['Volcanic degassing'] = [-30,13.7,5.3,5.962e-3,3.919e-2,\
                    0.20149,0.20141,0.30950,2.851e-2,2.512e-2,0.5073,6.358,\
                    5.043,6.150,2.359,1.988]
                    
isotope_comps['Asteroids'] = [-30, 8.9,5.3,5.962e-3,3.919e-2,\
                    0.20149,0.20141,0.30950,2.851e-2,2.512e-2,0.5073,6.358,\
                    5.043,6.150,2.359,1.988]
                    
isotope_comps['Comets'] = [1000,13.7,5.4,6.470e-3,4.124e-2,0.20629,0.20340,\
                        0.29915,2.947e-2,2.541e-2,0.50873,6.287,4.9958, \
                        6.0479,2.1288,1.6634]
                    
isotope_comps['IDPs'] = [np.nan, 13.7,5.8,6.470e-3,4.124e-2,0.20629,0.20340,\
                        0.29915,2.947e-2,2.541e-2,0.50873,6.287,4.9958, \
                        6.0479,2.1288,1.6634]


isotope_comps['Mars'][1]=1./isotope_comps['Mars'][1]
isotope_comps['Mars'][2]=1./isotope_comps['Mars'][2]

isotope_comps['Volcanic degassing'][1]=1./isotope_comps['Volcanic degassing'][1]
isotope_comps['Volcanic degassing'][2]=1./isotope_comps['Volcanic degassing'][2]

isotope_comps['Asteroids'][1]=1./isotope_comps['Asteroids'][1]
isotope_comps['Asteroids'][2]=1./isotope_comps['Asteroids'][2]

isotope_comps['Comets'][1]=1./isotope_comps['Comets'][1]
isotope_comps['Comets'][2]=1./isotope_comps['Comets'][2]

isotope_comps['IDPs'][1]=1./isotope_comps['IDPs'][1]
isotope_comps['IDPs'][2]=1./isotope_comps['IDPs'][2]



# Kurokawa Fig. 2 normalisation: abundance = (mol of reference isotope per kg
# of body) / (solar mol per kg). Calibrated so that the present-day 6 mbar of
# CO2 reproduces the digitised observed C point (as for the abundance plot):
# mol_i per kg = y_i * solar_abundances[i] * FIG2_MOL_PER_KG.
# Check: Fig. 2 asteroids then contain 0.94 wt% CO2-equivalent, i.e. the
# X_gas = 1% case, and CI-like 0.26 wt% C.
def _fig2_mol_per_kg(p_co2_present=6.e-3, area=4.*np.pi*3389.5e3**2,
                     g=3.721, m_mars=6.39e23):
    co2_moles = p_co2_present*1e5*area/g/44.01e-3
    return co2_moles/(solar_abundances[0]*abundances1['Obs'][1][0])/m_mars


FIG2_MOL_PER_KG = _fig2_mol_per_kg()


def comet_co2_mass_fraction():
    """Degassed CO2-equivalent mass fraction of comets from Kurokawa
    Fig. 2 (~2.9 wt%). Kurokawa §2.2: comet volatiles come from Fig. 2
    (Dauphas 2003 noble gases; C and N from solar C/Ar), independent of the
    asteroid parameter X_gas."""
    return abundances1['Comets'][1][0]*solar_abundances[0]*FIG2_MOL_PER_KG*44.01e-3


def set_amount_of_element(carbon_moles):
    """
        -from moles of total carbon, work out moles of carbon-12
        -from moles of carbon-12 you can calculate moles of others
        -and convert to mass using atomic weights
    """
    
    n_moles = np.zeros(len(atom_wgts))
    # carbon-12
    n_moles[0] = carbon_moles * carbon_12_pc/100.
    for i in range(1,len(atom_wgts)-1):
        # N, Ne, Ar, Kr
        n_moles[i] = carbon_moles * carbon_12_pc/100. * \
            abundances1['Volcanic degassing'][1][i] / abundances1['Volcanic degassing'][1][0] * \
            solar_abundances[i] / solar_abundances[0]
    
    # Xe
    n_moles[len(atom_wgts)-1] = carbon_moles * carbon_12_pc/100. * \
        abundances1['Mars'][1][len(atom_wgts)-1] / abundances1['Mars'][1][0] * \
        solar_abundances[len(atom_wgts)-1] / solar_abundances[0]


    # molecular nitrogen moles
    n2_moles = (n_moles[1]/2.) / (nitrogen_14_pc/100.)
    
    # scale all moles so that they equal the initial number
    initial_number = carbon_moles
    
    denominator = (carbon_moles+n2_moles+np.sum(n_moles[2:]))
    # adjusted carbon
    carbon_moles = carbon_moles * initial_number / denominator
    # adjusted nitrogen
    n2_moles = n2_moles * initial_number / denominator
    # adjusted everything else
    n_moles = n_moles * initial_number / denominator
    
    
    
    
    # now return the total mass of carbon dioxide and molecular nitrogen
    # and also return the total number of moles of the chosen isotopes 
    
    return ([carbon_moles,n2_moles],n_moles)
    
def _weights(ratios):
    """Isotope fractions including the reference isotope (weight 1)."""
    v = np.concatenate(([1.0], np.asarray(ratios, dtype=float)))
    return v / np.sum(v)


# Inventory basis for each isotope group [N, Ne, Ar, Kr, Xe].
# The N ledger (mole_elements[1]) holds TOTAL N atoms. The noble-gas ledgers
# hold REFERENCE-ISOTOPE moles (20Ne, 36Ar, 84Kr, 130Xe): they are built from
# Lodders reference-isotope solar abundances and Kurokawa Fig. 2, which is
# plotted per reference isotope. Mixing must use the matching basis, otherwise
# reservoirs are weighted by 1/sum(ratios), which differs between atmosphere
# and source (for Xe, Mars-like ~34 vs chondritic ~23, a 47% weighting error).
REFERENCE_BASIS = [False, True, True, True, True]


def _counts(ratios, inventory, reference_basis):
    """Isotope populations [reference, others...] for a ledger amount."""
    v = np.concatenate(([1.0], np.asarray(ratios, dtype=float)))
    return inventory * (v if reference_basis else v / np.sum(v))


def _default_basis(n_groups):
    """Groups are assumed to start at N, except that four groups means the
    noble gases only (the IDP call); pass basis explicitly to be safe."""
    if n_groups == len(REFERENCE_BASIS) - 1:   # noble gases only (IDPs)
        return list(REFERENCE_BASIS[1:])
    if 1 <= n_groups <= len(REFERENCE_BASIS):
        return list(REFERENCE_BASIS[:n_groups])
    raise ValueError('Cannot infer inventory basis; pass basis explicitly')


def _mix_one(ratios, source, inventory, amount, reference_basis=False):
    """Exact mixing of isotope populations for a positive source amount.

    inventory and amount must be on the same basis: total element atoms
    (reference_basis=False) or reference-isotope atoms (True).
    """
    if amount <= 0.0:
        return list(ratios)
    if inventory < 0.0:
        raise ValueError('Negative isotope-bearing elemental inventory')
    p = (_counts(ratios, inventory, reference_basis) +
         _counts(source, amount, reference_basis))
    return list(p[1:] / p[0])


def _escape_one(ratios, inventory, loss, factors, reference_basis=False):
    """Finite isotope-selective removal, integrated in safe isotope substeps.

    The escaping composition is proportional to isotope abundance times
    the relative diffusion/photochemical factor (Kurokawa Eq. 9).
    With reference_basis=False, `loss` is the total element loss. With
    reference_basis=True, `inventory` and `loss` are reference-isotope
    amounts (as computed by Eq. 6 from the reference-isotope ledger), so the
    reference isotope loses exactly `loss` and the other isotopes lose
    loss/inventory * N_k * R_k.
    """
    if loss <= 0.0:
        return list(ratios)
    if inventory <= 0.0 or loss >= inventory:
        raise ValueError(f'Escape exceeds inventory: N={inventory:.6g}, loss={loss:.6g}')
    counts = _counts(ratios, inventory, reference_basis)
    relative = np.concatenate(([1.0], np.asarray(factors, dtype=float)))
    if reference_basis:
        loss = loss * np.sum(counts * relative) / counts[0]
        if loss >= np.sum(counts):
            raise ValueError('Escape exceeds inventory')
    remaining_loss = float(loss)
    for _ in range(1000):
        if remaining_loss <= 1.e-12 * loss:
            break
        abundance_weighted = counts * relative
        rate_fraction = abundance_weighted / np.sum(abundance_weighted)
        # Do not remove over half any populated isotope in a substep.
        step = min(remaining_loss, 0.5 * float(np.min(counts / rate_fraction)))
        if step <= 0.0 or not np.isfinite(step):
            raise ValueError('Invalid isotope escape substep')
        counts -= step * rate_fraction
        remaining_loss -= step
    else:
        raise ValueError('Isotope escape requires excessive subcycling')
    return list(counts[1:] / counts[0])


def escape_fr(Imn, N, dN, dN2, Rdiffs_isotopes, basis=None):
    """Kurokawa (2018) Eq. 9: sputtering, then photochemical nitrogen.

    N, dN and dN2 are element ledger amounts: total N atoms for nitrogen,
    reference-isotope atoms for the noble gases (see REFERENCE_BASIS).
    Photochemical N has R_14/15 = 1.4 * diffusive R_14/15, so
    R_15/14 (stored here) is diffusive R_15/14 / 1.4.
    """
    basis = _default_basis(len(Imn)) if basis is None else basis
    result = []
    for m in range(len(Imn)):
        sputtered = _escape_one(Imn[m], N[m], dN[m], Rdiffs_isotopes[m],
                                basis[m])
        if m == 0:
            photofactor = np.asarray(Rdiffs_isotopes[m]) / 1.4
            sputtered = _escape_one(sputtered, N[m] - dN[m], dN2[m],
                                    photofactor, basis[m])
        elif dN2[m] != 0:
            raise ValueError('Photochemical escape only implemented for nitrogen')
        result.append(sputtered)
    return result


def continuous_sources(Imn, Imn0, N, dN, basis=None):
    """Kurokawa Eq. 10 integrated as exact inventory mixing."""
    basis = _default_basis(len(Imn)) if basis is None else basis
    return [_mix_one(r, src, ni, amount, b) for r, src, ni, amount, b
            in zip(Imn, Imn0, N, dN, basis)]


def impact_replenishment(Imn, Imn0_asteroids, Imn0_comets, N, dN, dN2,
                         basis=None):
    """Kurokawa Eq. 11, sequential exact asteroid and comet mixing."""
    mixed = continuous_sources(Imn, Imn0_asteroids, N, dN, basis)
    return continuous_sources(mixed, Imn0_comets,
                              np.asarray(N) + np.asarray(dN), dN2, basis)


def evolve_open_system(ratios, inventory, k_channels, factor_channels,
                       supplies, dt, reference_basis=False):
    """Exact solution over one step of the isotope population balance

        dP_j/dt = s_j - k_j P_j,   k_j = sum_c k_c f_{c,j}  (f = 1 for the
                                                             reference isotope)

    with constant per-particle escape rates k_c (s^-1) for each loss channel
    c, relative isotope factors f_{c,j} (e.g. R_diff, or R_diff/1.4 for N
    photochemistry) and constant supply rates s_j from the listed sources.
    This is Kurokawa Eqs. 9-10 combined and integrated analytically:

        P_j(dt) = P_j(0) exp(-k_j dt) + S_j (1 - exp(-k_j dt)) / (k_j dt),

    where S_j is the total supplied over the step. Unlike applying escape
    and then supply as separate sub-steps, it stays accurate when the escape
    lifetime is much shorter than dt: it then relaxes to the supply-escape
    steady state P_j = s_j / k_j, which is enriched in heavy isotopes by
    1/R_j, instead of stripping the inventory and resetting it to the source
    composition.

    supplies: list of (amount, source_ratios) on the same basis as inventory.
    Returns (new_ratios, new_inventory, removed_by_channel), where
    removed_by_channel[c] is the ledger-basis amount lost through channel c.
    """
    P0 = _counts(ratios, inventory, reference_basis)
    S = np.zeros_like(P0)
    for amount, src in supplies:
        if amount > 0.:
            S = S + _counts(src, amount, reference_basis)
    n = len(P0)
    full_f = [np.concatenate(([1.0], np.asarray(f, dtype=float)))
              for f in factor_channels]
    k = np.zeros(n)
    for kc, fc in zip(k_channels, full_f):
        k = k + kc * fc
    kdt = k * dt
    decay = np.exp(-kdt)
    with np.errstate(divide='ignore', invalid='ignore'):
        growth = np.where(kdt > 1e-12, -np.expm1(-kdt) / kdt, 1. - 0.5 * kdt)
    P = P0 * decay + S * growth
    removed = P0 + S - P
    removed_by_channel = []
    for kc, fc in zip(k_channels, full_f):
        with np.errstate(divide='ignore', invalid='ignore'):
            share = np.where(k > 0., kc * fc / k, 0.)
        r = removed * share
        removed_by_channel.append(r[0] if reference_basis else r.sum())
    if P[0] <= 0.:
        return list(ratios), 0., removed_by_channel
    new_inventory = P[0] if reference_basis else P.sum()
    return list(P[1:] / P[0]), new_inventory, removed_by_channel


if __name__ == "__main__":
    """
        plot
    """
    
    
    import matplotlib.pyplot as plt
    plt.ion()
    plt.plot(abundances1['Comets'][0],abundances1['Comets'][1])
    plt.plot(abundances1['IDPs'][0],abundances1['IDPs'][1])
    plt.plot(abundances1['Asteroids'][0],abundances1['Asteroids'][1])
    plt.plot(abundances1['Mars'][0],abundances1['Mars'][1])
    plt.plot(abundances1['Volcanic degassing'][0],abundances1['Volcanic degassing'][1])

    plt.legend(abundances1.keys())
    plt.yscale('log')   
    plt.xticks([12,14,20,36,84,130])    
    plt.gca().set_xticklabels(['C','N','Ne','Ar','Kr','Xe'])  
