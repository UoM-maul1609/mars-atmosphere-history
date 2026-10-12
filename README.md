# Mars atmosphere evolution model

A Python Monte Carlo model of the atmospheric pressure, elemental inventories and isotope-ratio evolution of Mars from **4.5 billion years ago to the present**, based primarily on **Kurokawa et al. (2018)**. This is an independently maintained research/teaching implementation, **not a validated replacement for published results**.

## Scientific model

The model advances an atmospheric state on nominal **1 Myr** steps, starting with approximately **1.5 bar** of atmosphere. It evolves CO₂, N₂ and H₂O atmospheric masses alongside elemental C, N, Ne, Ar, Kr and Xe inventories and selected isotope ratios. A prescribed initial noble-gas and nitrogen isotope composition is used. The steps implement:

1. **Impact erosion and delivery:** impactor sizes follow a broken power-law distribution with a break at 300 km; impact speeds are sampled from a digitised distribution. An integrated crater chronology weights the sampled projectiles through time. Erosion is calculated using the impact loss formulae, and asteroid/comet material delivers CO₂-equivalent volatiles, nitrogen and trace elements with prescribed source compositions. The present implementation uses a fixed impactor sample for the full trajectory and clips diameters at 2000 km.
2. **Sputtering:** prescribed CO₂ loss versus EUV history is used to estimate elemental escape with mass-dependent isotope fractionation. Sputtering is only applied after the user-specified dynamo shutoff epoch.
3. **Photochemical escape:** carbon escape depends on the time-dependent EUV history; nitrogen escape is parameterised using digitised Fox–Dalgarno reaction dependences. The implementation retains its historical nitrogen rate and molecular/atomic convention pending verification against the source rates.
4. **Interplanetary dust particles:** supply a continuous Ne/Ar/Kr/Xe source, with specified isotopic signatures.
5. **Volcanic degassing:** prescribed time-dependent CO₂, N₂ and H₂O sources contribute atmospheric mass, trace elements and isotopes.
6. **Atmospheric collapse / saturation:** the `obliquity_flag=-1` baseline treats collapse below 0.5 bar as a residual **6 mbar CO₂ partial pressure**, following a simplified approach in Kurokawa et al. Obliquity-dependent modes (`1`, `2`, `3`) use digitised Forget et al. temperature information to estimate a CO₂ saturation-pressure limit. These branches retain the existing code's experimental approximations. **There is no prognostic condensed-CO₂ reservoir**; do not interpret this as closed-system carbon conservation across collapse or reinflation.

The atmospheric scale height used in the impact calculations assumes a fixed **210 K** reference temperature. The model also computes a pressure-dependent global-mean temperature for diagnostic output, not as feedback into escape dynamics. H₂O is capped using an ice saturation-vapour-pressure relation at the reference temperature.

**Isotope representation:** the code advances ratios using the formulae corresponding to Kurokawa et al. equations 9–11, not individually conserved isotope inventories. The evolution of C/N is dynamically synchronised to CO₂/N₂ in the corrected version; the noble gases retain the paper-inspired abundance conventions. Model outputs should be compared with observations before asserting agreement.

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

In IPython, use `%run batch_runs.py` (rather than `%run -i`); the new script does not depend on the interactive namespace. A batch may take time: each trajectory integrates 4500 steps and evaluates a large impactor sample.

Results are saved in your operating system’s temporary directory, under a user-specific subfolder (typically `/tmp/<username>/` on macOS/Linux; Python resolves the actual location through `tempfile.gettempdir()`):

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
- Updated batch plotting to use raw math-text strings (no invalid escape warnings), avoid double counting collapsed trajectories, label units and colours, and write plots/data into the user-specific temporary directory.
- `main.run_model(..., output_file='path.npz')` is now an explicit option; importing the module no longer interprets unrelated command-line arguments as output file paths.

These corrections preserve the published-inspired parameterisations as far as practical; they can change the numerical results. Validate against the tagged original using matched seeds and inputs.

## Limitations and outstanding scientific checks

- **Photochemical nitrogen escape**: confirm whether each helper return rate counts escaping N atoms or N₂ molecules and whether the extra factor `Natm[1] / Ninit` is justified. Until that is resolved, δ¹⁵N predictions are provisional.
- **Collapse bookkeeping**: CO₂ can be removed/assigned by saturation rules without tracking an atmospheric–surface transfer; a physically conserved surface reservoir is a future model extension.
- **Time stepping**: forward-Euler source and isotope updates are not adaptively error-controlled. Test 1 Myr versus smaller steps and impactor sample size/seed convergence.
- **Impactor distribution**: the 2000 km cap creates a pile-up at the upper diameter limit, rather than a formally truncated distribution.
- **Observational fit**: the illustrative 45° batch showed 3/10 cases reaching low pressures; nitrogen and xenon isotope ratios did not simultaneously reproduce observed values. These are scientific discrepancies to investigate, not plotting bugs.
- The legacy `python/batch_process.py` is retained for historical reference but is not the supported batch workflow; use `batch_runs.py`.

## Repository layout

- `python/main.py`: model timestep integration and state bookkeeping.
- `python/batch_runs.py`: seed-controlled ensembles, diagnostics and plots.
- `python/isotopic_data.py`: composition tables and isotope equations.
- `python/impactor_model.py`, `impactor_velocity.py`, `impactor_atmos_loss_gain_eq4_5.py`, `crater_chronology.py`: impacts and chronology.
- `python/sputtering_co2.py`, `MarsSputtering_Calc6_7.py`, `photo_chemical.py`: atmospheric escape.
- `python/IDPs.py`, `volcano_outgassing.py`: continuous sources.
- `python/mars_temp_forget.py`, `mars_temperature.py`, `svp_ice.py`: temperature and vapour pressure.
- `references/`: bundled PDFs; `images/`: historical example plots.

## References

- Kurokawa et al. (2018), *Icarus*, **299**, 443–459, https://doi.org/10.1016/j.icarus.2017.08.020. Principal atmospheric and isotope model framework.
- Forget et al. (2018), supplied as `references/Forget_et_al_2018.pdf`: Martian climate/CO₂ stability constraints.
- Gierasch and Toon (1973), *Journal of the Atmospheric Sciences*, supplied as `references/Gierasch_and_Toon_JAS1973.pdf`: early Martian climate context.

Further rate curves and abundance datasets originate from works identified in comments in the respective Python helpers (including Fox and Dalgarno for nitrogen photochemical escape). Consult those sources before changing physical coefficients.
