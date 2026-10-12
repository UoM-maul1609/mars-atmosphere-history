"""Run seeded model ensembles and plot observational comparisons.

Run from the python/ directory with `python batch_runs.py` or in IPython
with `%run batch_runs.py`. Edit the CONFIG section for your experiment.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import getpass
import tempfile

import matplotlib.pyplot as plt
import numpy as np

from isotopic_data import abundances1, element_mass, errorbars1, solar_abundances
from main import run_model

# ------------------------------------------------------------------
# CONFIG: edit here. Seeds are always 0, 1, ..., NUM_RUNS - 1.
# ------------------------------------------------------------------
NUM_RUNS = 10
NUM_CORES = 10
OBLIQUITY = 2  # -1: fixed threshold, 1: 0°, 2: 45°, 3: 90°
SPUTTERING = True
PHOTOCHEMICAL_ESCAPE = True
COLLAPSED_PRESSURE_THRESHOLD = 0.1  # bar, classification for plotting ONLY
# Match the original script: /tmp/<username>/ on macOS/Linux.
OUTPUT_DIRECTORY = Path(tempfile.gettempdir()) / getpass.getuser()
MODEL_KWARGS = dict(
    sputtering_flag=SPUTTERING,
    pce_flag=PHOTOCHEMICAL_ESCAPE,
    C_Ne_IDP=1.0,
    f_comet=0.0005,
    X_gas=0.01,
    crater_model=1,
    dynamo_time=4.1,
)


def run_one(seed):
    return run_model(seed, OBLIQUITY, **MODEL_KWARGS)


def plot_trajectories(output, collapsed, folder):
    labels = [
        'Pressure (bar)',
        r'$\delta^{15}\mathrm{N}$ (‰)',
        r'$^{20}\mathrm{Ne}/^{22}\mathrm{Ne}$',
        r'$^{38}\mathrm{Ar}/^{36}\mathrm{Ar}$',
        r'$^{86}\mathrm{Kr}/^{84}\mathrm{Kr}$',
        r'$^{136}\mathrm{Xe}/^{130}\mathrm{Xe}$',
    ]
    fig, axes = plt.subplots(6, 1, figsize=(9, 14), sharex=True,
                             constrained_layout=True)
    obs_indices = [None, 0, 1, 2, 3, 4]
    for i, (t, y, elements, isotopes, area) in enumerate(output):
        age_gyr = t / 1e9
        values = [y[:, 0], y[:, 1], 1 / y[:, 2], y[:, 3], y[:, 4], y[:, 5]]
        style = dict(color='tab:blue', lw=1.3, alpha=.7) if collapsed[i] else dict(
            color='0.5', lw=.65, alpha=.45)
        for ax, value in zip(axes, values):
            ax.plot(age_gyr, value, **style)
    for k, ax in enumerate(axes):
        ax.set_ylabel(labels[k])
        if k == 0:
            ax.set_yscale('log')
            ax.text(.02, .08, f'{sum(collapsed)}/{len(output)} collapsed',
                    transform=ax.transAxes)
        else:
            obs = errorbars1[obs_indices[k]]
            yerr = np.array([[obs[1] - obs[0]], [obs[2] - obs[1]]])
            ax.errorbar([0.0], [obs[1]], yerr=yerr, fmt='ro', capsize=5,
                        label='Observed present day' if k == 1 else None)
        ax.grid(alpha=.18)
        ax.set_xlim(-4.5, .12)  # leave room for present-day observation bars
    axes[-1].set_xlabel('Time relative to present (Gyr)')
    obliquity_label = {-1: 'Fixed collapse threshold', 1:'0°', 2:'45°', 3:'90°'}.get(
        OBLIQUITY, str(OBLIQUITY))
    fig.suptitle(f'Mars atmospheric evolution | obliquity: {obliquity_label} | '
                 f'sputtering: {SPUTTERING} | PCE: {PHOTOCHEMICAL_ESCAPE}')
    filename = folder / 'atmospheric_evolution.png'
    fig.savefig(filename, dpi=200)
    plt.close(fig)
    return filename


def plot_abundances(output, collapsed, folder):
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for i, (t, y, elements, isotopes, area) in enumerate(output):
        abundance = elements / solar_abundances / 1e16
        style = dict(color='tab:blue', lw=1.3, alpha=.7) if collapsed[i] else dict(
            color='0.5', lw=.65, alpha=.45)
        ax.plot(element_mass, abundance, '-s', ms=3, **style)
    ax.plot(element_mass, abundances1['Obs'][1], 'ro-', lw=1.5, ms=6,
            label='Observed (digitised Kurokawa et al. fig. 6)')
    ax.set_yscale('log')
    ax.set_xlabel('Reference isotope mass number')
    ax.set_ylabel('Abundance relative to solar (model normalisation)')
    ax.set_xticks(element_mass, ['C', 'N', 'Ne', 'Ar', 'Kr', 'Xe'])
    ax.set_title('Present-day element abundances')
    ax.grid(alpha=.3)
    ax.legend()
    filename = folder / 'elemental_abundances.png'
    fig.savefig(filename, dpi=200)
    plt.close(fig)
    return filename


def main():
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    workers = max(1, min(NUM_CORES, NUM_RUNS))
    if workers == 1:
        output = [run_one(i) for i in range(NUM_RUNS)]
    else:
        # Executor.map returns results in input seed order, independent of runtime.
        with ProcessPoolExecutor(max_workers=workers) as executor:
            output = list(executor.map(run_one, range(NUM_RUNS)))
    collapsed = [bool(y[-1, 0] < COLLAPSED_PRESSURE_THRESHOLD)
                 for _, y, _, _, _ in output]
    paths = [plot_trajectories(output, collapsed, OUTPUT_DIRECTORY),
             plot_abundances(output, collapsed, OUTPUT_DIRECTORY)]
    np.savez_compressed(
        OUTPUT_DIRECTORY / 'batch_results.npz',
        t=output[0][0],
        ystore=np.stack([result[1] for result in output]),
        mole_elements=np.stack([result[2] for result in output]),
        seeds=np.arange(NUM_RUNS),
        collapsed=np.array(collapsed),
        obliquity=OBLIQUITY,
    )
    print(f'Completed {NUM_RUNS} runs; {sum(collapsed)} classified as collapsed '
          f'(final pressure < {COLLAPSED_PRESSURE_THRESHOLD} bar).')
    for path in paths:
        print(f'Saved: {path}')
    print(f'Saved: {OUTPUT_DIRECTORY / "batch_results.npz"}')


if __name__ == '__main__':
    main()
