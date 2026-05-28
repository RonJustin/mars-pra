"""
Mars Mission Probabilistic Risk Assessment — v9.3
==================================================
Independent exploratory PRA framework for external technical discussion.

Author : Patrick Jauslin — Quantix Analytics Consulting (quantixac.com)
Version: 9.3 — May 2026
License: MIT

IMPORTANT DISCLAIMER
--------------------
This framework is an independent analytical exercise. Prior distributions
are parameterised from publicly available analogue data and expert
elicitation. They are NOT empirically calibrated against Mars-class
operational records (none exist). All outputs are structured plausibility
estimates, not validated engineering specifications.

Dependencies
------------
    pip install numpy scipy matplotlib

Usage
-----
    python mars_pra_v93.py              # run full simulation, save figures + CSV
    python mars_pra_v93.py --seed 123   # reproducible with custom seed
    python mars_pra_v93.py --n 50000    # faster run with fewer samples
"""

import argparse
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import beta, lognorm, norm as snorm, t as t_dist, gaussian_kde
import csv

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# 0.  CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

DEFAULTS = dict(seed=42, n=100_000, out_dir='.')

PALETTE = dict(
    navy='#1C2030', gold='#C8A86C', blue1='#185FA5', blue2='#378ADD',
    blue3='#B5D4F4', amber='#BA7517', red='#A32D2D', green='#3B6D11',
    purp='#534AB7', teal='#0F6E56', white='#FFFFFF', dgray='#555555',
)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': PALETTE['white'], 'axes.facecolor': PALETTE['white'],
    'axes.edgecolor': '#CCCCCC', 'axes.labelcolor': PALETTE['navy'],
    'xtick.color': PALETTE['navy'], 'ytick.color': PALETTE['navy'],
    'text.color': PALETTE['navy'],
    'grid.color': '#E8E8E8', 'grid.linewidth': 0.8,
})

# ─────────────────────────────────────────────────────────────────────────────
# 1.  PRIOR DISTRIBUTIONS
#     All parameterisations are documented with source rationale.
#     See Appendix A of the companion memorandum for full derivation.
# ─────────────────────────────────────────────────────────────────────────────

PRIOR_SPEC = {
    # ── Solar / Radiation ──────────────────────────────────────────────────
    # Distribution: Beta (bounded [lo, hi])
    # Mode: 3.5% — basic shielding (polyethylene + water wall, ~20 g/cm²),
    #              mission timing to solar minimum.
    # Lower (1.0%): enhanced shielding + active magnetic, solar minimum.
    # Upper (15.0%): minimal shielding (~5 g/cm² Al-equiv.), solar maximum.
    # Sources: Zeitlin et al. (2004) [MARIE, ~0.22 mGy/day GCR cruise];
    #          Hassler et al. (2013) [CRATER, 1.84 mSv/day MSL cruise];
    #          Cucinotta (2014) [LD50/30 ~3–4 Gy, SEP lethality model];
    #          Schwenn (2006) [solar cycle event frequency].
    # Note: Full LOC derivation requires GEANT4/HZETRN transport model
    #       for Starship-specific shielding geometry. These values are
    #       scenario approximations.
    'solar': dict(
        type='beta', mode=3.5, lo=1.0, hi=15.0, conc=12.0,
        scenario_b=dict(mode=10.0, lo=4.0, hi=20.0),
    ),

    # ── Engineering / Propulsion ───────────────────────────────────────────
    # Distribution: Beta
    # Mode: 4.5% — novel propulsion + structural system, some ISS analogue
    #              heritage but no Starship deep-space precedent.
    # Lower (1.5%): analogous to mature shuttle-era systems [Fragola 1996].
    # Upper (13.0%): novel system with no deep-space flight heritage.
    # Sources: Fragola, Maggio & Collins (1996) [ISS PRA];
    #          NASA/TM-2011-216468 [ISS PRA overview].
    'engineering': dict(
        type='beta', mode=4.5, lo=1.5, hi=13.0, conc=12.0,
        scenario_b=dict(mode=9.0, lo=4.0, hi=16.0),
    ),

    # ── ECLSS / Life-support ───────────────────────────────────────────────
    # Distribution: Beta (with rank correlation to Engineering, rho=0.30)
    # Mode: 2.0% — redundant architecture assumed; ISS ECLSS uptime >95%
    #              for 90-day missions.
    # Lower (0.5%): fully mature heritage, multiple demonstrated closures.
    # Upper (9.0%): epistemic upper bound reflecting absence of any
    #               270-day closed-loop deep-space operational record.
    # Sources: Hendrickx et al. (2006) [MELiSSA Phase IV];
    #          ISS ECLSS operational data (90-day maximum validated duration);
    #          Biosphere 2 (1991) as failure-mode precedent for O2 closure loss.
    # CRITICAL: No 270-day continuous closed-loop deep-space ECLSS record
    #           exists as of May 2026. This prior cannot be updated via Bayes
    #           until such data is collected.
    'eclss': dict(
        type='beta', mode=2.0, lo=0.5, hi=9.0, conc=12.0,
        corr_with='engineering', rho_a=0.30, rho_b=0.35,
        scenario_b=dict(mode=4.5, lo=2.0, hi=10.0),
    ),

    # ── Human Factors / Psychology ────────────────────────────────────────
    # Distribution: Beta
    # Mode: 0.8% — optimised crew selection, psychological pre-screening,
    #              270-day isolation with structured task schedule.
    # Lower (0.2%): near-ideal conditions, high crew cohesion.
    # Upper (8.0%): full 20-minute communication delay, irreversibility
    #               stress, prolonged confinement without abort option.
    # Sources: Mars-500 analogue (520-day isolation, 2010–2011);
    #          Antarctic winter-over cohort studies (conflict incidence).
    'human_factors': dict(
        type='beta', mode=0.8, lo=0.2, hi=8.0, conc=12.0,
        scenario_b=dict(mode=5.0, lo=1.5, hi=12.0),
    ),

    # ── Debris Collision (4–10 cm) ────────────────────────────────────────
    # Distribution: LogNormal (heavy-tailed, reflecting unmeasured regime)
    # Median: 0.024% — phase-space flux integral, heliocentric gradient
    #                   r^-1.5 corrected (from v9.x framework, Equations 3.1–3.2).
    # p95: 1.8% — dense population upper bound (power-law alpha=3.5).
    # LogNormal chosen over Beta because:
    #   (a) the 4–10 cm regime in interplanetary space has never been
    #       directly measured; the distribution tail is genuinely unknown;
    #   (b) the ratio p95/median ≈ 75 spans orders of magnitude, inconsistent
    #       with the bounded support of Beta distributions;
    #   (c) LogNormal is the standard choice for failure rates with high
    #       epistemic uncertainty in PRA literature [NRC 2008].
    # Sources: Grün et al. (1985) [meteoroid flux model, d<1mm];
    #          Fechtig et al. (1981) [heliocentric gradient validation];
    #          NRC (2008) [LogNormal for epistemic uncertainty in PRA].
    # Validation gap: gradient r^-1.5 empirically validated for d<1mm only.
    #                 Extrapolation to 4–10 cm is unvalidated.
    'debris': dict(
        type='lognormal', median=0.024, p95=1.8,
        scenario_b=dict(median=0.8, p95=2.5),
    ),
}

# ─────────────────────────────────────────────────────────────────────────────
# 2.  SAMPLING FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def _beta_params(mode_pct, lo_pct, hi_pct, conc=12.0):
    """Fit Beta(a,b) to mode and range via moment matching."""
    lo, hi, mode = lo_pct/100, hi_pct/100, mode_pct/100
    span = hi - lo
    mode_s = np.clip((mode - lo)/span, 0.05, 0.95)
    a = mode_s*(conc-2)+1
    b = conc - a
    return max(a, 1.01), max(b, 1.01), lo, span


def sample_beta(mode_pct, lo_pct, hi_pct, conc=12.0, n=100_000):
    """Draw n samples from Beta prior fitted to mode and range."""
    a, b, lo, span = _beta_params(mode_pct, lo_pct, hi_pct, conc)
    return (lo + beta(a, b).rvs(n)*span)*100  # → % units


def sample_lognormal(median_pct, p95_pct, n=100_000):
    """Draw n samples from LogNormal prior with given median and p95."""
    mu = np.log(median_pct)
    sigma = (np.log(p95_pct) - mu)/1.645
    return np.clip(lognorm(s=sigma, scale=np.exp(mu)).rvs(n), 1e-4, 6.0)


def apply_gaussian_copula(s1, s2, rho, n=100_000):
    """Apply Gaussian copula rank correlation between two sample arrays."""
    u1 = np.argsort(np.argsort(s1))/n
    u2 = np.argsort(np.argsort(s2))/n
    z1 = snorm.ppf(np.clip(u1, 1e-6, 1-1e-6))
    z2 = snorm.ppf(np.clip(u2, 1e-6, 1-1e-6))
    z2c = rho*z1 + np.sqrt(1-rho**2)*z2
    u2n = snorm.cdf(z2c)
    return np.sort(s2)[(u2n*n).astype(int).clip(0, n-1)]


def apply_t_copula(s1, s2, rho, nu=4, n=100_000):
    """Apply t-copula (heavier tail dependence than Gaussian)."""
    u1 = np.argsort(np.argsort(s1))/n
    u2 = np.argsort(np.argsort(s2))/n
    z1 = t_dist.ppf(np.clip(u1, 1e-6, 1-1e-6), df=nu)
    z2 = t_dist.ppf(np.clip(u2, 1e-6, 1-1e-6), df=nu)
    z2c = rho*z1 + np.sqrt(1-rho**2)*z2
    u2n = t_dist.cdf(z2c, df=nu)
    return np.sort(s2)[(u2n*n).astype(int).clip(0, n-1)]


def total_loc(*factors_pct):
    """Compute total LOC using P_total = 1 - prod(1 - p_i)."""
    s = np.ones(len(factors_pct[0]))
    for f in factors_pct:
        s *= (1 - f/100)
    return (1 - s)*100


# ─────────────────────────────────────────────────────────────────────────────
# 3.  MAIN SIMULATION
# ─────────────────────────────────────────────────────────────────────────────

def run_simulation(seed=42, n=100_000, out_dir='.'):
    np.random.seed(seed)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f'Running simulation: seed={seed}, N={n:,}')

    # ── Draw prior samples ──────────────────────────────────────────────────
    spec = PRIOR_SPEC
    solar_A  = sample_beta(spec['solar']['mode'],
                           spec['solar']['lo'], spec['solar']['hi'], n=n)
    eng_A    = sample_beta(spec['engineering']['mode'],
                           spec['engineering']['lo'], spec['engineering']['hi'], n=n)
    eclss0   = sample_beta(spec['eclss']['mode'],
                           spec['eclss']['lo'], spec['eclss']['hi'], n=n)
    hf_A     = sample_beta(spec['human_factors']['mode'],
                           spec['human_factors']['lo'], spec['human_factors']['hi'], n=n)
    debris_A = sample_lognormal(spec['debris']['median'], spec['debris']['p95'], n=n)

    # Scenario B (pessimistic shift)
    sb=spec['solar']['scenario_b']; solar_B=sample_beta(sb['mode'],sb['lo'],sb['hi'],n=n)
    eb=spec['engineering']['scenario_b']; eng_B=sample_beta(eb['mode'],eb['lo'],eb['hi'],n=n)
    eb2=spec['eclss']['scenario_b']; eclss0_B=sample_beta(eb2['mode'],eb2['lo'],eb2['hi'],n=n)
    hb=spec['human_factors']['scenario_b']; hf_B=sample_beta(hb['mode'],hb['lo'],hb['hi'],n=n)
    db=spec['debris']['scenario_b']; debris_B=sample_lognormal(db['median'],db['p95'],n=n)

    # ── Apply Gaussian copula correlation: Engineering–ECLSS ──────────────
    rho_A = spec['eclss']['rho_a']
    rho_B = spec['eclss']['rho_b']
    eclss_A = apply_gaussian_copula(eng_A, eclss0,   rho_A, n=n)
    eclss_B = apply_gaussian_copula(eng_B, eclss0_B, rho_B, n=n)

    # ── Compute LOC ────────────────────────────────────────────────────────
    loc_A = total_loc(solar_A, eng_A, eclss_A, hf_A, debris_A)
    loc_B = total_loc(solar_B, eng_B, eclss_B, hf_B, debris_B)

    # ── Percentile summary ─────────────────────────────────────────────────
    results = {}
    for label, arr in [('Scenario_A', loc_A), ('Scenario_B', loc_B)]:
        pcts = {p: float(np.percentile(arr, p)) for p in [5,10,25,50,75,90,95]}
        pcts['mean'] = float(arr.mean())
        pcts['std']  = float(arr.std())
        results[label] = pcts
        print(f'{label}: mean={pcts["mean"]:.2f}%  '
              f'p5={pcts[5]:.2f}%  p50={pcts[50]:.2f}%  p95={pcts[95]:.2f}%')

    # ── Variance decomposition (squared Spearman rank correlation) ─────────
    factor_samples_A = [solar_A, eng_A, eclss_A, hf_A, debris_A]
    factor_labels    = ['Solar/Radiation', 'Engineering', 'ECLSS', 'HumanFactors', 'Debris']
    var_raw = [np.corrcoef(s, loc_A)[0,1]**2 for s in factor_samples_A]
    var_sum = sum(var_raw)
    var_pct = {lbl: float(v/var_sum*100) for lbl, v in zip(factor_labels, var_raw)}
    results['variance_decomposition_A'] = var_pct
    print('Variance shares (Sc A):', {k: f'{v:.1f}%' for k,v in var_pct.items()})

    # ── Save CSV ────────────────────────────────────────────────────────────
    csv_path = out_dir / 'mars_pra_v93_samples.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['solar_A','eng_A','eclss_A','hf_A','debris_A','loc_A',
                         'solar_B','eng_B','eclss_B','hf_B','debris_B','loc_B'])
        for i in range(n):
            writer.writerow([
                f'{solar_A[i]:.4f}', f'{eng_A[i]:.4f}', f'{eclss_A[i]:.4f}',
                f'{hf_A[i]:.4f}', f'{debris_A[i]:.6f}', f'{loc_A[i]:.4f}',
                f'{solar_B[i]:.4f}', f'{eng_B[i]:.4f}', f'{eclss_B[i]:.4f}',
                f'{hf_B[i]:.4f}', f'{debris_B[i]:.6f}', f'{loc_B[i]:.4f}',
            ])
    print(f'CSV saved: {csv_path}')

    # ── Save JSON summary ───────────────────────────────────────────────────
    json_path = out_dir / 'mars_pra_v93_summary.json'
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'JSON saved: {json_path}')

    return {
        'loc_A': loc_A, 'loc_B': loc_B,
        'factors_A': factor_samples_A,
        'factor_labels': factor_labels,
        'var_pct': var_pct,
        'results': results,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 4.  COPULA SENSITIVITY ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def copula_sensitivity(seed=42, n=100_000, out_dir='.'):
    """Compare Gaussian vs t-copula across rho range 0.0–0.5."""
    np.random.seed(seed)
    out_dir = Path(out_dir)

    solar  = sample_beta(3.5, 1.0, 15.0, n=n)
    eng    = sample_beta(4.5, 1.5, 13.0, n=n)
    eclss0 = sample_beta(2.0, 0.5,  9.0, n=n)
    hf     = sample_beta(0.8, 0.2,  8.0, n=n)
    debris = sample_lognormal(0.024, 1.8, n=n)

    rhos = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    sensitivity = {'gaussian': {}, 't_nu4': {}}

    for rho in rhos:
        for key, fn in [('gaussian', apply_gaussian_copula),
                        ('t_nu4',    lambda s1, s2, r, n=n: apply_t_copula(s1, s2, r, nu=4, n=n))]:
            ec = fn(eng, eclss0.copy(), rho)
            loc = total_loc(solar, eng, ec, hf, debris)
            sensitivity[key][rho] = {
                'p5':  float(np.percentile(loc, 5)),
                'p50': float(np.percentile(loc, 50)),
                'p95': float(np.percentile(loc, 95)),
                'mean': float(loc.mean()),
            }

    json_path = out_dir / 'copula_sensitivity.json'
    with open(json_path, 'w') as f:
        json.dump(sensitivity, f, indent=2)
    print(f'Copula sensitivity saved: {json_path}')
    return sensitivity


# ─────────────────────────────────────────────────────────────────────────────
# 5.  ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Mars Mission PRA v9.3 — Quantix Analytics Consulting')
    parser.add_argument('--seed', type=int, default=DEFAULTS['seed'])
    parser.add_argument('--n',    type=int, default=DEFAULTS['n'],
                        help='Number of Monte Carlo samples')
    parser.add_argument('--out',  type=str, default=DEFAULTS['out_dir'],
                        help='Output directory for CSV, JSON, figures')
    args = parser.parse_args()

    sim = run_simulation(seed=args.seed, n=args.n, out_dir=args.out)
    sens = copula_sensitivity(seed=args.seed, n=min(args.n, 50_000), out_dir=args.out)

    print('\n── Copula sensitivity summary (p50 LOC) ──')
    print(f'{"rho":>5}  {"Gaussian p50":>14}  {"t-copula p50":>14}')
    for rho in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        g = sens['gaussian'][rho]['p50']
        t = sens['t_nu4'][rho]['p50']
        print(f'{rho:>5.1f}  {g:>13.2f}%  {t:>13.2f}%')

    print('\nDone. All outputs written to:', args.out)
