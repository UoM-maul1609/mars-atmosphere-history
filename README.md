# Mars atmosphere evolution model

A Python Monte Carlo model of the atmospheric pressure, elemental inventories and isotope-ratio evolution of Mars from **4.5 billion years ago to the present**, based primarily on **Kurokawa et al. (2018)**. This is an independently maintained research/teaching implementation, **not a validated replacement for published results**.

## Scientific model

The model advances an atmospheric state on nominal **1 Myr** steps, starting with approximately **1.5 bar** of atmosphere. It evolves CO₂, N₂ and H₂O atmospheric masses alongside elemental C, N, Ne, Ar, Kr and Xe inventories and selected isotope ratios. A prescribed initial noble-gas and nitrogen isotope composition is used. The steps implement:

1. **Impact erosion and delivery:** impactor sizes follow a broken power-law distribution with a break at 300 km; impact speeds are sampled from a digitised distribution. Impacts are **discrete events**. Each run draws N = 2×10²¹ kg / ⟨m⟩ impactors (about 2.6×10⁵ with D ≥ 3 km), gives each a time drawn from the crater chronology, and applies it in that 1 Myr step using the atmosphere at that moment. Erosion is calculated using the impact loss formulae, and asteroid/comet material delivers CO₂-equivalent volatiles, nitrogen and trace elements with prescribed source compositions. See *Discrete impacts* below.
2. **Sputtering:** the Luhmann et al. (1992) CO₂ sputtering rate, interpolated in log space through time, is used to estimate elemental escape with mass-dependent isotope fractionation. Sputtering is only applied after the user-specified dynamo shutoff epoch. Upper-atmosphere diffusion uses Kurokawa et al. (2018) Table 2, interpolated against log(EUV).
3. **Photochemical escape:** carbon escape depends on the time-dependent EUV history; nitrogen escape is parameterised using digitised Fox–Dalgarno reaction dependences. Photochemical escape operates throughout the model history, independently of the dynamo; its nitrogen rate is treated as escaping N atoms/s and converted to N₂ molecular inventory using a factor of 1/2. The original extra nitrogen-inventory scaling has been removed because the source ratio is already included in the digitised rate functions.
4. **Interplanetary dust particles:** supply a continuous Ne/Ar/Kr/Xe source, with specified isotopic signatures.
5. **Volcanic degassing:** prescribed time-dependent CO₂, N₂ and H₂O sources contribute atmospheric mass, trace elements and isotopes.
6. **Atmospheric collapse / saturation:** the `obliquity_flag=-1` baseline treats collapse below 0.5 bar as a residual **6 mbar CO₂ partial pressure**, following a simplified approach in Kurokawa et al. Obliquity-dependent modes (`1`, `2`, `3`) follow Bradley (MPhil thesis, §4.3.1): collapse is triggered when the CO₂ partial pressure exceeds the CO₂ saturation vapour pressure at the PlanetWRF 5th-percentile surface temperature for the current pressure; CO₂ is then held at 6 mbar, as in the baseline. See *Obliquity collapse modes* below for known limitations. **There is no prognostic condensed-CO₂ reservoir**; do not interpret this as closed-system carbon conservation across collapse or reinflation.

The atmospheric scale height used in the impact calculations assumes a fixed **210 K** reference temperature. The model also computes a pressure-dependent global-mean temperature for diagnostic output, not as feedback into escape dynamics. H₂O is capped using the ice saturation-vapour-pressure relation at 273.16 K in the warm regime and 210 K in the collapsed 0.5-bar baseline regime (Kurokawa section 2.1).

**Isotope representation:** the code advances ratios by mixing or removing explicit normalised isotope populations (including the reference isotope) corresponding to Kurokawa et al. equations 9–11. These finite-step updates preserve isotope populations during mixing and apply differential fractionation during escape; isotope inventories are not stored as independent state variables. The evolution of C/N is dynamically synchronised to CO₂/N₂ in the corrected version; the noble gases retain the paper-inspired abundance conventions. Model outputs should be compared with observations before asserting agreement.

## Getting started

Python 3.10+ is recommended. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cd python
python -c "import main; t, y, elements, isotope_ratios, area = main.run_model(0, obliquity_flag=2); print(y[-1, :6])"
```

Or work interactively:

```python
import main
import matplotlib.pyplot as plt

t, y, elements, ratios, area = main.run_model(0, obliquity_flag=2)
plt.plot(t / 1e9, y[:, 0])
plt.yscale('log')
plt.xlabel('Time relative to present (Gyr)')
plt.ylabel('Atmospheric pressure (bar)')
plt.show()
```

The first positional argument is the **random seed**. Use the same seed and parameters when comparing model versions. Function options include `sputtering_flag`, `pce_flag`, `X_gas`, `f_comet`, `C_vol`, `crater_model`, `C_Ne_IDP`, and `dynamo_time`. `obliquity_flag` is `-1` for a fixed threshold, `1` for 0°, `2` for 45°, and `3` for 90°. These cases are *alternative parameterisations*, not a time series of actual Mars obliquity.

`run_model` returns `(t, ystore, mole_elements, isotopes_sim, Amars)`. `ystore` columns are:

| Column | Variable | Units |
|---|---|---|
| 0 | total pressure | bar (numerically close to atm but not identical) |
| 1 | δ¹⁵N | ‰ |
| 2 | ²²Ne/²⁰Ne | dimensionless (invert for ²⁰Ne/²²Ne) |
| 3 | ³⁸Ar/³⁶Ar | dimensionless |
| 4 | ⁸⁶Kr/⁸⁴Kr | dimensionless |
| 5 | ¹³⁶Xe/¹³⁰Xe | dimensionless |
| 6 | diagnostic global temperature | K |
| 7 | atmospheric scale height | m |
| 8–10 | CO₂, N₂, H₂O atmospheric masses | kg |

`t` contains years relative to the present, from -4.5 × 10⁹ to 0. The elemental inventory in `mole_elements` is in moles, and `isotopes_sim` holds a nested list of isotope ratios in the order defined by `isotopic_data.isotopes`.

## Running batches and reproducing the figures

Edit the configuration variables at the top of `python/batch_runs.py` (`NUM_RUNS`, `NUM_CORES`, `OBLIQUITY`, and process switches). Then, from `python/`:

```bash
python batch_runs.py
```

In IPython, use `%run batch_runs.py` (rather than `%run -i`); the script does not depend on the interactive namespace. **Restart the kernel after updating any module.** `%run` re-executes `batch_runs.py` but keeps previously imported modules such as `isotopic_data` cached. On macOS the parallel workers start fresh and use the new code while the plotting process does not, which produces mismatched plots. Each trajectory takes about 3–4 s.

Results are saved in `/tmp/$USER/` on Linux and macOS (created if missing). `tempfile.gettempdir()` is deliberately not used there, because on macOS it resolves to a per-session `/var/folders/...` path. On Windows the system temporary directory is used, under a user-named subfolder:

- `atmospheric_evolution.png`: six-panel atmospheric pressure and isotope-ratio histories, with present-day observations shown as red points and error bars.
- `elemental_abundances.png`: final model element abundances compared with digitised observations (Kurokawa et al., figure 6).
- `batch_results.npz`: numerical time series, elemental inventories, run seeds and collapse classification; use this to replot without rerunning.

**Blue curves** are runs ending below the *plotting classification threshold* of 0.1 bar; **grey curves** finish above it. That 0.1-bar diagnostic is **not** the physical model's CO₂ collapse condition. The observed present-day δ¹⁵N and other isotope constraints are plotted in red. There is no attempt to fit parameters to the observations automatically.

The supplied original reference plots are in `images/`. They are historical illustrations, not regression-test targets for this corrected version.

## Corrections in this version

Compared with the original pre-review archive:

- Corrected the extra Avogadro conversion in photochemical carbon escape.
- Eliminated stale cometary contributions being reused as photochemical isotope escape when that process is disabled.
- Corrected impactor mass from metre diameters and integrated the velocity distribution over its nonuniform velocity grid.
- Reconciled initial CO₂/N₂ molecular masses with specified initial pressure and dynamic elemental C/N inventories.
- Removed a second addition of volcanic carbon to the elemental ledger after atmospheric CO₂ had already been updated.
- Updated batch plotting to use raw math-text strings (no invalid escape warnings), avoid double counting collapsed trajectories, label units and colours, and write plots/data into `/tmp/$USER/`.
- `main.run_model(..., output_file='path.npz')` is now an explicit option; importing the module no longer interprets unrelated command-line arguments as output file paths.

These corrections preserve the published-inspired parameterisations as far as practical; they can change the numerical results. Validate against the tagged original using matched seeds and inputs.

## Limitations and outstanding scientific checks

- **Nitrogen escape** has been verified against Fox & Dalgarno (1980, 1983), Jakosky et al. (1994) and Luhmann et al. (1992), and reproduces Kurokawa Fig. 14. See *October 2026: nitrogen*.
- **Collapse bookkeeping**: CO₂ can be removed/assigned by saturation rules without tracking an atmospheric–surface transfer; a physically conserved surface reservoir is a future model extension.
- **Time stepping**: escape and continuous supply of N and the noble gases are integrated exactly over each step for constant rates, and erosion compounds per impact. The rates themselves are held fixed within each 1 Myr step, so test against smaller steps before publication.
- **Impactor distribution**: the size distribution is truncated at `D_max` (default: the diameter of a single body carrying the whole budget, about 1140 km) and renormalised.
- **Observational fit**: with discrete impacts and the `run_model` defaults, 4/40 runs collapse with the fixed threshold and 2/40 at 45°. All of them collapse at 4.26–4.44 Ga, before the 4.1 Ga constraint from ALH84001. N abundance and ¹³⁶Xe/¹³⁰Xe are not reproduced. These are scientific discrepancies to investigate, not plotting bugs.
- **Batch configuration**: `batch_runs.py` uses Kurokawa's best-fit settings: `X_gas = 0.01`, `C_vol = 5`, `C_Ne_IDP = 10`, `f_comet = 0.001` and `crater_model = 2`.
- The legacy `python/batch_process.py` is retained for historical reference but is not the supported batch workflow; use `batch_runs.py`.

## Repository layout

- `python/main.py`: model timestep integration and state bookkeeping.
- `python/batch_runs.py`: seed-controlled ensembles, diagnostics and plots.
- `python/isotopic_data.py`: composition tables and isotope equations.
- `python/impactor_model.py`, `impactor_velocity.py`, `impactor_atmos_loss_gain_eq4_5.py`, `crater_chronology.py`: impacts and chronology.
- `python/sputtering_co2.py`, `MarsSputtering_Calc6_7.py`, `photo_chemical.py`: atmospheric escape.
- `python/IDPs.py`, `volcano_outgassing.py`: continuous sources.
- `python/mars_temp_forget.py`, `mars_temperature.py`, `svp_ice.py`: temperature and vapour pressure.
- `references/`: bundled PDFs (Kurokawa 2018, Forget 2013, Gierasch & Toon 1973, Fox & Dalgarno 1980 and 1983, Jakosky et al. 1994, Luhmann et al. 1992, Morbidelli et al. 2012); `images/`: historical example plots.

## References

- Kurokawa et al. (2018), *Icarus*, **299**, 443–459, https://doi.org/10.1016/j.icarus.2017.08.020. Principal atmospheric and isotope model framework.
- Forget et al. (2013), *Icarus*, **222**, 81–99, https://doi.org/10.1016/j.icarus.2012.10.019, supplied as `references/Forget_et_al_2013.pdf` (previously mislabelled 2018): global-mean temperature vs pressure (Fig. 1) and permanent-CO₂-cap regimes vs obliquity (Fig. 14).
- Bradley, MPhil thesis (University of Manchester): PlanetWRF 0°/45°/90° obliquity simulations at 0.006, 0.5 and 1.0 bar (Fig. 4.3), the source of the obliquity collapse lookup tables in `mars_temp_forget.py`.
- Gierasch and Toon (1973), *Journal of the Atmospheric Sciences*, supplied as `references/Gierasch_and_Toon_JAS1973.pdf`: early Martian climate context.

Further rate curves and abundance datasets originate from works identified in comments in the respective Python helpers (including Fox and Dalgarno for nitrogen photochemical escape). Consult those sources before changing physical coefficients.

## Paper-alignment revision (October 2026)

Checked directly against the supplied Kurokawa et al. (2018) PDF:

| Paper location | Current code |
|---|---|
| §2.1, Eq. 1 | Warm/cold H₂O saturation temperatures (273.16/210 K) for the fixed-collapse baseline |
| Eq. 6, Table 1 | CO₂ (44 amu) and N₂ (28 amu) in molecular sputtering yield and diffusion |
| Eq. 7, Table 2 | EUV-dependent homopause–exobase separation and temperature; log-EUV interpolation between tabulated points is an explicit numerical choice |
| §2.3–2.4 | Sputtering begins after dynamo termination; photochemical escape runs throughout |
| §2.4 | Nitrogen escape based on digitised N₂/CO₂ functions, with no additional initial-N multiplier |
| Eqs. 9–11, §2.7 | Explicit reference isotope in mixing/escape normalisation; for photochemical nitrogen, R(15/14)=Rdiff(15/14)/1.4 |

The sputtering Eq. 6 flux uses CO₂ and N₂ **molecular** inventories, converted to two N atoms for elemental loss and isotope fractionation. (The former 1 Myr Euler step with a 50% per-species escape cap has been replaced by an exact open-system solve; see *October 2026: nitrogen*.) The code uses finite-step isotope population mixing/removal rather than the exact forward-Euler transcription of the published differential equations. This is numerically preferable but requires convergence tests for 1 Myr steps. The published EUV upper-atmosphere Table 2 is sampled by interpolation in log(EUV); the paper does not specify precisely which interpolation was used. The `obliquity_flag=1,2,3` collapse modes are not part of Kurokawa et al.; see *Obliquity collapse modes* below.

**Limit of source validation:** Fox & Dalgarno (1980, 1983), Jakosky et al. (1994) and Luhmann et al. (1992) have now been checked (see *October 2026: nitrogen*). Morbidelli et al. (2012) has also been checked (see *October 2026: crater chronology*). Pepin (1994), Terada et al. (2009) and Zhang et al. (1993) have not been verified against primary copies. Avoid claiming exact independent reproduction of those papers.

## Tests for this revision

Run `python -m unittest discover -s tests -v` from the root. Twenty automated checks cover:

- Kurokawa Table 2 nodes, molecular sputtering reference masses and the Table 1 CO₂ yield;
- the CO₂-only sputtering limit, fractionation signs and photochemical nitrogen selectivity;
- isotope-conserving mixing on both inventory bases and the IDP reference-isotope conversion;
- the truncated size sampler, the impact count against Kurokawa Fig. 3a, chronology-weighted impact timing, and seed reproducibility;
- the exact open-system escape and supply solve (steady state, mass balance, mixing limit), the bulk-ratio photochemical N lookup at both solar-activity levels, the Ribas-timed sputtering nodes, and N and Ar lifetimes against Kurokawa Fig. 14.

Ensembles of 40 seeds were also run for `obliquity_flag = -1` and `2` (each run about 3.5 s). These are execution smoke tests, **not** quantitative reproduction of the paper’s published ensembles.

## Reproducibility warning

The physical corrections materially alter predicted isotope trajectories. Do not treat the old fitted parameter combinations as calibrated for this branch. Benchmark against the original version using matching seeds, evaluate isotope/molecular elemental conservation, timestep convergence and Monte Carlo uncertainty, then compare ensemble statistics with the paper's figures.

## October 2026: isotope diffusion sign correction

The preceding paper-aligned version mistakenly calculated isotope diffusion mass difference as `reference_mass - isotope_mass`. With `R_diff = exp(-delta_m*g*delta_z/(k_B*T))`, that preferentially removed **heavier** isotopes (e.g. 15N over 14N), leading to artificial light-isotope enrichment. The revised function now uses `isotope_mass - reference_mass`, so the heavy isotope has `0 < R_diff < 1` for positive height and temperature. A regression test checks all modeled isotope groups against Kurokawa Equation 7. This is a targeted correction, **not** full validation or a fitted reproduction of the published ensemble.

## October 2026: obliquity collapse modes

### Corrections

Checked against Bradley's thesis (§4.3.1, Fig. 4.3), Kurokawa et al. (2018) §2.1 and Forget et al. (2013):

- **Collapse outcome.** Before, CO₂ was clipped every step to svp(T₅(P)) instead of collapsing. The thesis describes a collapse to present-day conditions. With the digitised tables, that clip settled at pCO₂ ≈ 2 mbar at 90° and ≈ 5 mbar at 45°, below the observed 6 mbar. Now, once triggered, the collapse is sticky and holds CO₂ at 6 mbar, as in the baseline.
- **Partial pressure.** The trigger compared *total* pressure (CO₂ + N₂ + H₂O) with the CO₂ saturation pressure. It now uses the CO₂ partial pressure. The total pressure is still used to look up T₅.
- **Water vapour after collapse.** The obliquity modes always capped H₂O at the 273.16 K saturation pressure (≈ 6 mbar), even after collapse. Collapsed runs therefore ended with as much water vapour as CO₂ (e.g. 8–11 mbar total with only 2–5 mbar CO₂). All modes now use 210 K after collapse (Kurokawa §2.1).
- **Patm floor.** `Patm = max(P, 6 mbar)` changed `Patm` without changing `Matm`, so the T₅ lookup could use a pressure that did not match the state. This floor has been removed. The baseline's 6 mbar collapse is unchanged.
- Removed the dead `Obliquity = 20` / `meanPvObliquity` initialisation; it was always overwritten.

### Known methodological limitations (not changed)

These are scientific issues with the 5th-percentile approach. Address them before drawing obliquity conclusions:

1. **The 5th percentile misses the cold trap.** Permanent CO₂ caps form at the coldest, polar, grid points. The area poleward of 80° is about 3% of the surface, so an area-weighted 5th percentile does not sample the polar caps. This is probably why the 0° tables never trigger collapse at any pressure (svp(T₅) > P everywhere, e.g. 46 mbar at P = 6 mbar). That contradicts Forget et al. (2013, §5.4 and Fig. 14), where *low* obliquity favours permanent caps: a permanent cap forms at 0.5 bar with 10° obliquity but not at 25°. The thesis conclusion that low obliquity resists collapse is most likely an artefact of this metric.
2. **The criterion is circular where CO₂ ice exists.** The flat floors of the 45° and 90° 5th-percentile curves (≈ −125 °C at 6 mbar, ≈ −82 to −85 °C at 0.5 bar) are the CO₂ frost point (−124 °C and −88 °C). Surfaces carrying CO₂ ice are pinned there. svp(T₅) ≈ P is then true by construction, and whether collapse triggers depends on a few kelvin of digitisation error.
3. **Seasonal ice versus permanent caps.** Forget et al. diagnose collapse from CO₂ ice that survives the year. A seasonal-temperature percentile cannot separate seasonal frost from perennial caps.
4. **Sparse table, linear extrapolation.** Only three pressures (0.006, 0.5, 1.0 bar) are used, with linear extrapolation outside them (for example to the 1.5 bar initial state). T₅ is linear in P between 6 mbar and 0.5 bar, although temperature depends roughly on log P.
5. **No CO₂ reservoir.** Collapsed CO₂ is deleted, not stored, so the atmosphere cannot reinflate when conditions change.

### Suggested improvements

- From the WRF output, diagnose the **annual-minimum (or annual-mean) temperature of the coldest grid cell**, or the **net annual CO₂ surface-ice accumulation** (thesis Figs 4.9–4.10). Trigger collapse when perennial CO₂ ice accumulates, which is the Forget et al. criterion.
- Alternatively, use Forget et al. (2013) Fig. 14 directly. It is already digitised as `data['collapse1']` / `meanPvObliquity` in `mars_temp_forget.py` but unused. It gives a collapse pressure threshold as a function of obliquity, and replaces Kurokawa's fixed 0.5 bar with an obliquity-dependent value.
- Interpolate the temperature tables in log P, and add WRF runs at 0.1–0.3 bar and 2 bar, where the thresholds lie.
- Treat obliquity as time-varying (chaotic, mean ≈ 38° ± 14°; Laskar et al. 2004) and sample it stochastically per run or per epoch, rather than fixing it for 4.5 Gyr.
- Add a polar-cap CO₂ reservoir so that collapse and reinflation conserve carbon.

## October 2026: equation audit against Kurokawa et al. (2018)

### Verified as matching the paper

- **Eq. 2, 4 and 5** (impact ξ, atmospheric loss polynomial, projectile loss). Eq. 3 is normalised by projectile mass, following Shuvalov (2009); the printed `m_atm` in Kurokawa Eq. 3 appears to be a typo.
- **Eq. 6 and 7** (sputtering with dilution factor α, diffusive fractionation with the corrected sign), Table 1 yields, and the Table 2 upper atmosphere.
- **Eq. 8** photochemical carbon escape, and the R_ph = 1.4 R_diff nitrogen photochemical fractionation.
- **Eq. 9–11** isotope escape, continuous sources and impact replenishment, given the basis fix below.
- **Table 3** IDP rates, given the conversion fix below.
- **Table 4** isotope compositions, every entry. The 20/22 and 36/38 ratios are stored inverted on purpose.
- **Volcanic degassing:** cumulative CO₂ is 0.41 bar for `C_vol = 1`, matching the paper's 0.4 bar. The N₂ rates for `C_vol = 5` match Fig. 12a.
- **Size distribution:** q = 3.5 below 300 km and q = 2 above, continuous at the break.
- **Crater chronology N(1):** Ivanov/Neukum form. δ¹⁵N is referenced to air ¹⁵N/¹⁴N = 3.6765 × 10⁻³.

### Fixed in this revision

1. **IDP supply units.** Element mass rates (g/Gyr) were divided by the reference-isotope mass and added in full to reference-isotope ledgers (²⁰Ne, ³⁶Ar, ⁸⁴Kr, ¹³⁰Xe). That overstated the IDP supply about 23× for Xe, 1.75× for Kr, 1.18× for Ar and 1.08× for Ne. `IDPs.get_mole` now returns reference-isotope moles, using the mean atomic mass of the IDP composition.
2. **Isotope mixing basis.** The noble-gas ledgers hold reference-isotope moles, but mixing treated them as total-element moles. Each reservoir was therefore weighted by 1/Σratios, which differs between the atmosphere and the source: a 47% weighting error for Mars-like versus chondritic Xe. `isotopic_data.REFERENCE_BASIS` now records the basis of each group (N is total, the noble gases are reference isotope), and mixing and escape use the matching basis.
3. **CO₂ sputtering yield.** `Y_CO2` was 0.67 while Table 1 gives 0.7, which raised every non-CO₂ sputtering rate by 4%.
4. **CO₂ sputtering history.** Luhmann's F_CO₂ (3 × 10²⁶, 6 × 10²⁵ and 3 × 10²³ s⁻¹) was interpolated linearly in time. That gave 2.4 × 10²⁵ s⁻¹ at 1 Ga (80× present) for an EUV of only about 1.4×. It is now interpolated in log F, giving 2.5 × 10²⁴ s⁻¹, about 10× less sputtering through the collapsed era.
5. **EUV history.** EUV is now the Ribas et al. (2005) τ^−1.23 law with saturation before a solar age of 0.1 Gyr (§2.3), instead of a log-linear extrapolation of Luhmann's points (12× instead of about 100× at 4.5 Ga, 9× instead of 17× at 4.1 Ga). This affects the Table 2 lookup; photochemical N still caps EUV at 6× (see below).
6. **Lyman-continuum clock (Eq. 8).** F_Lyc ∝ t^−0.86 now uses the solar age, not time since 4.5 Ga, which made the flux singular at the start.
7. **Comet onset.** Comets were switched on at `dynamo_time`, so changing the dynamo epoch also moved the LHB. They now use a separate `comet_time = 4.1` (§2.2).
8. **Plotted observations.** δ¹⁵N was drawn as 620 (420–800) instead of 572 ± 82, and ³⁸Ar/³⁶Ar as 0.20 ± 0.05 instead of 1/(4.2 ± 0.1) = 0.238 ± 0.006.
9. **Abundance-plot normalisation.** The hard-coded 1e16 made the model curves sit about 2.7× high against the digitised observations. The factor is now calibrated so that 6 mbar of CO₂ reproduces the observed C point (about 2.7e16). N and ³⁶Ar give 3.4e16 and 2.1e16 as cross-checks.

Effect on a collapsed best-fit run (seed 0, `C_vol = 5`, `C_Ne,IDP = 10`, `f_comet = 0.1%`):

| Quantity | Before | After | Observed |
|---|---|---|---|
| δ¹⁵N | 441‰ | 599‰ | 572 ± 82‰ |
| ²⁰Ne/²²Ne | 11.0 | 10.2 | 10.1 ± 0.7 |
| ³⁶Ar/³⁸Ar | 4.49 | 4.68 | 4.2 ± 0.1 |

### Outstanding discrepancies (not changed; need decisions or primary sources)

1. ~~Impacts are not discrete events.~~ **Resolved:** see *Discrete impacts* below.
2. ~~Nitrogen under-retained by about 100×.~~ **Resolved:** see *October 2026: nitrogen*.
   - The φ curves are normalised to 1 at x = 0.0264, the present-day *bulk* N₂ fraction. The code multiplies by the exobase enrichment R_diff (about 17) before the lookup, which looks like double counting. Using the bulk ratio improves N/C 14×, but it remains about 9× low.
   - N₂ sputtering then dominates.
   - The original code hid this by multiplying photochemical escape by `N2/N2_initial`, a factor with no physical basis that the previous review removed.
   - Fox & Dalgarno (1983) and Jakosky et al. (1994) are needed to fix the normalisation. Until then, treat N abundances and δ¹⁵N as provisional.
3. **Comet volatiles scale with `X_gas`.** Comet CO₂ delivery is `f_comet × X_gas × retained mass`, but Kurokawa sets comet volatile content from Fig. 2 (C about 3× CI per unit mass), independent of `X_gas`. With `X_gas = 1%` this probably undersupplies cometary noble gases by more than an order of magnitude, so calibrated `f_comet` values are not comparable with the paper.
4. ~~Default chronology.~~ **Resolved:** checked against Morbidelli et al. (2012); see *October 2026: crater chronology*. Collapses still occur at about 4.3–4.4 Ga, before the 4.1 Ga ALH84001 constraint; this is not caused by the chronology.
5. ~~Photochemical N EUV cap.~~ **Resolved:** the cap was removed (October 2026: nitrogen).
6. **Euler erosion.** Erosion removes Σf of the atmosphere in one step instead of 1 − exp(−Σf). The 50% cap never binds in the tested seeds; the largest step loss was 48%, against 38% if compounded.
7. **Minor approximations.** Surface gravity is used at the homopause–exobase level (g is about 8% lower at 150 km). Reference-isotope masses stand in for mean element masses in Eq. 6.

## October 2026: discrete impacts

Kurokawa treats impacts as instantaneous stochastic events (§2.2; Fig. 12 caption, after Kurosawa 2015). The previous code applied one fixed sample of 10⁵ impactors in **every** 1 Myr step, scaled by the crater flux, so each impactor was smeared over 4.5 Gyr. Seeds then differed only in the size sample, and a sample that lacked a giant body was rescaled to 2×10²¹ kg using many more small, highly erosive bodies.

### Method (`impactor_model.sample_impact_history`, `main.run_model`)

1. **Number of impacts.** Following Kurokawa §2.2, N = M_total / ⟨m⟩, with ⟨m⟩ evaluated analytically for the truncated broken power law (q = 3.5 below 300 km, q = 2 above, D ≥ 3 km). The total mass is 2×10²¹ kg on average but varies between runs: the median is about 0.9×10²¹ kg, because a few large bodies carry most of the mass.
2. **Size cap.** The default `D_max` is the diameter of one body carrying the whole budget (about 1140 km). This reproduces Kurokawa Fig. 3a: about 2.3–2.6×10⁵ impacts with D ≥ 3 km, about 1.2×10⁴ with D ≥ 10 km, and largest bodies of a few hundred to about 1000 km.
   - The old 2000 km clip piled every draw above 2000 km onto a 1.1×10²² kg body, five times the whole budget. Truncating at 2000 km instead would give only 9×10⁴ impacts.
   - `D_max` can be set explicitly.
3. **Timing.** Each impact gets a timestep drawn with probability proportional to the crater-chronology rate in that step (`crater_model`), plus a velocity from the digitised distribution.
4. **Comets.** After `comet_time`, each impact is a comet with probability f/(1 + f), so `f_comet` is now the comet-to-asteroid **number** ratio as Kurokawa defines it. Before, it was the fraction of delivered volatile mass.
5. **Application.** All impacts in a step use the atmosphere at the start of that step. Erosion compounds as ∏(1 − f_j), each f_j ≤ 0.99, which replaces the old linear sum with a 50% cap. No single impact removes more than about 3×10⁻⁴ of a 0.5–1.5 bar atmosphere.
6. **Reproducibility.** Every random draw comes from `RandomState(runNo)`, so a seed reproduces its history exactly. Pass `impact_log=[]` to `run_model` to receive (time, D, v, is_comet, eroded fraction) for every impact.

Runs are also about 7× faster (about 3.5 s), because only the impacts in each step are evaluated, not 10⁵ per step.

### Results (40 seeds, `run_model` defaults)

| Mode | Collapsed | Collapse time | Collapsed runs above 0.5 bar at 4.1 Ga |
|---|---|---|---|
| Fixed 0.5 bar threshold | 4/40 | 4.26–4.38 Ga | 0 |
| 45° obliquity | 2/40 | 4.30–4.36 Ga | 0 |

- **What decides collapse.** It is now governed by the delivery side. Collapsed runs drew no large volatile-rich body: their median largest impactor is about 340 km, against about 790 km in surviving runs, and their median total impactor mass is about 0.5×10²¹ kg, against 1.7×10²¹ kg. Erosion is the sum of about 2.6×10⁵ small impacts, so it varies little between runs.
- **Ratios in collapsed runs** (crater model 1): ²⁰Ne/²²Ne is about 10.2, which matches. δ¹⁵N (about 1000‰) and ³⁶Ar/³⁸Ar (4.3–5.0) scatter widely around the observations, and ¹³⁶Xe/¹³⁰Xe is about 1.95.

## October 2026: nitrogen

A 100-run batch showed collapsed runs ending with N 10²–10⁴× below observed, and δ¹⁵N shooting to 4,000–10,000‰ after volcanism stopped at 0.31 Ga. The cause was N₂ being lost far too fast. The N ledger itself is consistent: N atoms equal 2 × N₂ to machine precision throughout.

Every N-related number has now been checked against the primary sources:
- Fox & Dalgarno (1980, *PSS* 28, 41);
- Fox & Dalgarno (1983, *JGR* 88, 9027);
- Jakosky, Pepin, Johnson & Fox (1994, *Icarus* 111, 271);
- Luhmann, Johnson & Zhang (1992, *GRL* 19, 2151);
- and against Kurokawa et al.'s Fig. 12a (N₂ budget terms) and Fig. 14 (N and Ar lifetimes).

### Fixes

1. **Photochemical lookup uses the bulk N₂/CO₂ ratio** (`photo_chemical.nitrogen_escape`).
   - Fox & Dalgarno (1983), Fig. 3 caption: "Relative escape fluxes as a function of the N₂/CO₂ mixing ratio **at the turbopause**. The flux at the current value of the N₂ mixing ratio, 0.026, is arbitrarily set equal to 1.0". Jakosky et al. (1994) also use the homopause ratio.
   - The code multiplied the bulk ratio by R_diff (about 17) before the lookup, counting diffusive enrichment twice. Present-day composition gave φ₁ = 4.4. `EXOBASE_RATIO_LOOKUP = True` restores the old behaviour.
2. **N₂ sputtering yield halved.** Jakosky et al. (1994) Eq. 12 uses the Table II yield ratio "with the N yield divided by 2 in this case since we are tracking the loss of N₂ molecules". The Table II yields are in N *atoms* ejected. The code used 2.4/0.7 per N₂ molecule, which doubled N₂ sputtering.
3. **Corrected Luhmann CO₂ sputtering rates** (`sputtering_co2.CO2_SPUTTERING_SOURCE`).
   - Luhmann et al. (1992) Table 1 gives 3×10²³, 6×10²⁵ and 3×10²⁶ CO₂/s at 1×, 3× and 6× EUV. Jakosky et al. (1994) Table I gives values "corrected" after Johnson and Luhmann found a numerical error: 2×10²³, 3×10²⁵ and 1×10²⁶.
   - Kurokawa's Fig. 14 lifetimes match the corrected values. Together with fix 2 they reproduce the N sputtering lifetimes (about 3× and 6× the old values at 0 and 3.5 Ga) and the Ar ones (1.5× and 3×).
   - The default is `'jakosky_corrected'`; `'luhmann_published'` is available.
4. **Fox & Dalgarno solar-activity level** (`photo_chemical.FD83_SOLAR_ACTIVITY`).
   - The code's six process fluxes are a correct transcription of the Fox & Dalgarno (1983) Table 3 *solar-cycle averages*, with the right φ₁/φ₂/φ₃ grouping. The table gives low, high and average values (totals 2.3, 8.9 and 5.6 ×10⁵ cm⁻² s⁻¹).
   - With the **low-activity** values, the photochemical lifetimes reproduce Kurokawa Fig. 14a (about 2, 0.55 and 0.15 Gyr at 0, 2 and 3.5 Ga, 6 mbar). The average values are 2–3× shorter. Fox & Dalgarno (1980) also give a present rate of 2.0×10⁵.
   - The default is `'low'`. Jakosky et al. (1994) used Fox (1993a)'s average-activity total of 4.8×10⁵, so `'average'` is available.
5. **No EUV cap on photochemical N escape.** `np.minimum(euv_flux, 6)` has been removed. Kurokawa §2.4 scales the rates as F_EUV (photon and electron processes, including dissociative recombination) and F_EUV^1.5 (ion–molecule reactions). Fig. 14a needs the full saturated EUV of about 100× at 4.5 Ga.
6. **Luhmann sputtering nodes at Ribas times** (3.49 and 2.69 Ga for 6× and 3× EUV), matching the kinks in Kurokawa Fig. 14b. Luhmann's Fig. 4 shows the sputtering rate linear in log against age, consistent with the log interpolation.
7. **Escape and supply solved together, exactly** (`isotopic_data.evolve_open_system`).
   - For N and each noble gas, every isotope population follows dP/dt = S − kP over the step, with constant escape and volcanic/IDP supply rates, giving P(Δt) = P₀e^(−kΔt) + (S/k)(1 − e^(−kΔt)).
   - This replaces escape, then IDP supply, then volcanic supply as separate explicit-Euler sub-steps with 50% caps. When escape lifetimes were shorter than the 1 Myr step, the old sequence stripped each inventory to an extremely fractionated remnant and then reset it to the source composition. This caused the δ¹⁵N spikes of 4,000–10,000‰.
   - Jakosky Eq. 5 gives the steady state this solution reaches (ratio = source ratio / R).
   - CO₂ escape is also integrated exponentially.

Also verified against Jakosky Table III: R(28/44) = 17.4, R(15/14) = 0.837, R(38/36) = 0.700 and R(36/44) = 4.2 at Δz/T = 0.4. The code gives the same values.

### Lifetimes against Kurokawa Fig. 14 (6 mbar, N₂/CO₂ = 10⁻², Ar/CO₂ = 10⁻⁶)

| | 3.5 Ga | 2 Ga | 0 Ga |
|---|---|---|---|
| N sputtering, Kurokawa (read-off) | ~0.005 Gyr | ~0.04 | ~2 |
| N sputtering, code | 0.0053 | 0.044 | 1.96 |
| Ar sputtering, Kurokawa | ~0.013 | – | ~6 |
| Ar sputtering, code | 0.014 | 0.14 | 6.0 |
| N photochemical, Kurokawa | ~0.15 | ~0.55 | ~2 |
| N photochemical, code | 0.20 | 0.73 | 1.6 |

### Results (40 seeds, fixed threshold, `run_model` defaults, no calibration factors; 4 collapsed)

| | Before | After | Observed |
|---|---|---|---|
| N abundance, model/obs | 10⁻⁴–10⁻² | 0.88 | 1 |
| N/C atoms | ~10⁻⁴ | 0.039 | 0.043 |
| δ¹⁵N | spikes to 4,000–10,000‰ | 465‰ (361–469) | 572 ± 82‰ |
| ³⁶Ar/³⁸Ar | spikes | 4.76 | 4.2 ± 0.1 |

### Diagnostics and options

- `run_model(..., rate_log=[])` records per step the time, volcanic, sputtered and photochemical N₂ fluxes (molecules/s), pN₂ and pCO₂, for comparison with Fig. 12a.
- `f_sp` (all sputtering) and `f_ph_N` (photochemical N) are optional multipliers, as in Hu & Thomas (2023, arXiv:2302.04241). Both default to 1 and are no longer needed to match N.

### Remaining

- **δ¹⁵N** in collapsed runs is about 100‰ below the observed range.
- **Short transients after collapse and dynamo shutdown,** within a few Myr of 4.1 Ga: ³⁸Ar/³⁶Ar up to about 1, δ¹⁵N up to about 1,200‰. They are Rayleigh fractionation of the inventory left from the thick atmosphere as it is sputtered away, before supply takes over.
- **Ne and Ar are low in collapsed runs.** This is mostly the comet-delivery issue: comet volatiles scale with `X_gas`.

## October 2026: crater chronology

`crater_chronology.py` was checked against Morbidelli et al. (2012, *EPSL* 355–356, 144), Fig. 3 caption. In the formulas below t is age in Gyr and rates are craters with D > 20 km per km² per Gyr.

| `crater_model` | Rate for t < 4.1 Gyr | Rate for t > 4.1 Gyr | Used by |
|---|---|---|---|
| 1 | 2.7×10⁻¹⁶e^(6.93t) + 5.9×10⁻⁷ ✓ | same exponential (Neukum & Ivanov 1994, the blue line in Morbidelli Fig. 3) | Bradley's thesis ("the blue curve") |
| 2 | as above ✓ | sawtooth, solid curve: 2.5×10⁻² exp[−((4.5 − t)/0.003)^0.34] ✓ | Kurokawa et al. (2018) |
| 3 | as above ✓ | sawtooth, dashed curve: 2.0×10⁻² exp[−((4.5 − t)/0.01)^0.5] ✓ | – |

Also checked: N(1) = 5.44×10⁻¹⁴(e^(6.93t) − 1) + 8.38×10⁻⁴t ✓ (Morbidelli Fig. 1).

Changes:
- **Cap removed.** Models 2 and 3 capped the rate at the Neukum & Ivanov value for 4.5 Ga. That cap is not in the paper and changed the post-4.1 Ga integral by 2–4%.
- **Default changed.** The default `crater_model` is now **2**, the Morbidelli sawtooth that Kurokawa used, in `run_model` and `batch_runs.py`. Use `crater_model=1` to reproduce Bradley's thesis setup.

Model 2 puts 19% of impacts after 4.1 Ga. This is consistent with Morbidelli's estimate that the sawtooth event accounts for about a quarter of the lunar bombardment; their Fig. 3 cumulative N₂₀ gives about 21%.

**Effect.** In 40 fixed-threshold seeds, 3 collapse, at 4.28–4.41 Ga, against 4 at 4.26–4.38 Ga with model 1. All three chronologies put 70–78% of impacts before 4.3 Ga. Runs that draw no large volatile-delivering body therefore lose their atmosphere then, whichever chronology is chosen.

Kurokawa's collapses cluster around 4 Ga and the ALH84001 constraint requires more than 0.5 bar at 4.1 Ga, so a difference in collapse timing remains. Candidates are the comet volatile delivery (see the equation audit) and impact-erosion details such as the fixed 210 K scale height. It is not the chronology.
