# Mars Mission PRA v9.3

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Independent%20Technical%20Memo-orange)]()
[![Monte Carlo](https://img.shields.io/badge/Monte%20Carlo-N%3D100%2C000-purple)]()
[![Version](https://img.shields.io/badge/Version-v9.3-red)]()

**Independent Exploratory Probabilistic Risk Assessment for Crewed Mars Missions**

> Prior distributions are expert-elicitation estimates, not empirically calibrated against Mars-class operational records. Results are structured plausibility estimates intended for external technical discussion — not engineering certification.

**Author:** Patrick Jauslin · [Quantix Analytics Consulting](https://quantixac.com) · May 2026  
**Companion memo:** `Mars_Risk_v9_3_Quantix_Final.docx`

---

## The Core Finding at a Glance

Published Mars collision risk figures (0.01–0.1%) are derived from **robotic probes** with cross-sections of 5–100 m². A crewed Starship has an effective collision area of approximately **2,600 m²**. Without explicit scale correction, published figures may understate crewed vessel exposure by up to **130×**.

| Mission Class | Cross-Section | Cumulative Exposure | Empirical Precedent |
|---|---|---|---|
| Robotic probes (historical) | 5–100 m² | ~250 m²·years | Yes — 50+ missions |
| Crewed Starship (100 pax) | ~2,626 m² | ~1,970 m²·years | **None** |
| Colony vessel (500 pax) | ~12,921 m² | ~9,691 m²·years | **None** |

A 500-person Starship mission carries approximately **38.8× the cumulative collision exposure** of all prior Mars robotic missions combined.

---

## Monte Carlo Results (N = 100,000)

| Scenario | p5 | p25 | **p50 (median)** | p75 | p95 | Mean |
|---|---|---|---|---|---|---|
| **Scenario A** (conservative-optimistic) | 8.8% | 10.9% | **12.6%** | 14.4% | 17.3% | 12.7% |
| **Scenario B** (conservative-pessimistic) | 22.7% | 25.4% | **27.2%** | 29.1% | 31.9% | 27.3% |
| Crew Dragon design standard | — | — | **~0.37%** | — | — | ~0.37% |

> Scenario A median is approximately **34× above the Crew Dragon 1-in-270 standard**.

---

## Variance Decomposition

Squared Spearman rank correlations (approximating first-order Sobol sensitivity indices):

| Risk Factor | Mode (Sc. A) | Distribution | Variance Share |
|---|---|---|---|
| Engineering / propulsion | 4.5% | Beta | **~33%** |
| Solar / radiation | 3.5% | Beta | **~29%** |
| ECLSS / life-support | 2.0% | Beta (correlated, ρ=0.30) | **~21%** |
| Human factors | 0.8% | Beta | ~6% |
| Debris (4–10 cm regime) | 0.024% | LogNormal (heavy tail) | **~10%** |

The LogNormal debris prior reflects **epistemic uncertainty in the unmeasured 4–10 cm interplanetary regime** — no in-situ measurements exist at 1.0–1.5 AU for this particle size class.

---

## Method

### Architecture

```
N = 100,000 Monte Carlo samples
├── Risk factors: Beta distributions (bounded probabilities)
├── Debris: LogNormal (heavy-tailed epistemic uncertainty)
├── Correlation: Gaussian copula, ρ = 0.30 (Engineering ↔ ECLSS)
└── Total LOC: P_total = 1 − ∏(1 − pᵢ)
```

### Debris Flux Model

Phase-space flux integration:

```
λ = ∫_trajectory ∫∫ A_eff(v_rel) · v_rel · f(r,v,d) d³v dt
```

Mean encounter velocity ~20 km/s (Öpik-Wetherill). Heliocentric radial gradient:

```
ρ(r) = ρ₁AU · (1 AU / r)^1.5
f_gradient ≈ 0.72 for 1.0–1.52 AU
```

Empirically validated for d < 1 mm (Grün 1985; Fechtig 1981). Extrapolation to 4–10 cm is **unvalidated** — the LogNormal prior captures this.

### Copula Sensitivity

| Copula | ρ = 0.0 | ρ = 0.1 | ρ = 0.3 | ρ = 0.5 |
|---|---|---|---|---|
| Gaussian (p50) | 12.6% | 12.4% | 12.2% | 12.1% |
| Gaussian (p95) | 17.3% | 17.3% | 17.3% | 17.2% |
| t-copula ν=4 (p50) | 12.6% | 12.6% | 12.5% | 12.5% |
| t-copula ν=4 (p95) | 17.3% | 17.3% | 17.3% | 17.4% |

**The p50 LOC varies by less than 0.5 percentage points across the full ρ range and copula family.** The dominant uncertainty source is the prior distributions themselves, not the correlation structure.

---

## Indicative Investment Priority Order

These reflect analytical inferences from the risk-weight structure — not design prescriptions.

| Priority | Component | Variance Addressed | Indicative Cost |
|---|---|---|---|
| **1A** | ECLSS closed-loop maturation (2-yr ISS analogue, MELiSSA extension) | ~21% of LOC variance | USD 2–5B |
| **1B** | Engineering heritage (3–5 uncrewed deep-space Starship flights, 270-day FMEA) | ~33% of LOC variance | USD 1–3B |
| **1C** | Radiation shielding + solar minimum timing (zero-cost timing lever: ~40% SEP reduction) | ~29% of LOC variance | USD 0.5–2B |
| **2** | Pathfinder debris missions (piezoelectric arrays at 1.0, 1.2, 1.4 AU) | ~10% of LOC variance + gradient validation | USD 0.75–3B |

This ordering is robust across copula family, correlation magnitude, and scenario.

---

## Epistemic Status of Key Inputs

| Input | Empirical Basis | Primary Source | Key Gap |
|---|---|---|---|
| GCR dose rate | Good (MARIE/CRATER data) | Zeitlin (2004); Hassler (2013) | Starship shielding transport model |
| SEP lethality | Partial (event freq. only) | Schwenn (2006); Cucinotta (2014) | Full LET spectrum + dose-rate model |
| Engineering failure rate | Partial (ISS/Shuttle analogue) | Fragola (1996); NASA PRA (2011) | No Starship deep-space heritage |
| **ECLSS 270-day deep-space** | **None direct** | MELiSSA (Hendrickx 2006) | **No 270-day closed-loop record exists** |
| Human factors (9-month) | Analogue only | Mars-500; Antarctic winter-over | Comm delay + irreversibility unmodelled |
| **Debris (4–10 cm)** | **None (interplanetary)** | Grün (1985) for d < 1 mm | **No in-situ measurements at 1–1.5 AU** |
| Gradient r⁻¹·⁵ | Partial (d < 1 mm only) | Grün (1985); Fechtig (1981) | Not validated for 4–10 cm |

The two `None direct` entries represent the **primary structural gaps**. Because no Mars-class likelihood function exists for these factors, the priors are effectively the posteriors until dedicated data programs are completed.

---

## Crewed Spaceflight LOC Context

| Programme | LOC (est.) | Empirical Record | Primary Risk Driver |
|---|---|---|---|
| Mercury (early) | ~40% | 0 LOC / 6 crewed | Unknown systems; no heritage |
| Apollo (lunar) | ~1–15% | 1 LOC / 17 crewed | Launch vehicle; mission complexity |
| Space Shuttle | ~0.5–1.5% | 2 LOC / 135 flights | O-ring failure; foam impact |
| ISS missions | ~0.5–1% | 0 LOC (to date) | Mature systems; abort-to-Earth available |
| Crew Dragon | < 0.37% | 0 LOC / 7+ missions | CtCap 1-in-270 standard (Fragola 1996) |
| **Mars Sc. A (v9.3)** | **median 12.6% [8.8–17.3%]** | No precedent | Engineering + solar/radiation |
| **Mars Sc. B (v9.3)** | **median 27.2% [22.7–31.9%]** | No precedent | Exceeds all documented tolerances |

---

## Usage

```bash
pip install numpy scipy matplotlib
python mars_pra_v93.py --seed 42 --n 100000 --out ./output/
```

**Dependencies:**

```
numpy>=1.24
scipy>=1.11
matplotlib>=3.7
```

**Output files** (written to `--out` directory):

| File | Contents |
|---|---|
| `loc_distributions.png` | Full Monte Carlo LOC distributions (Sc. A + B) |
| `prior_derivation.png` | Prior distribution panels with analogue anchors |
| `variance_decomposition.png` | Squared Spearman rank correlations (tornado chart) |
| `copula_sensitivity.png` | Gaussian vs. t-copula comparison across ρ range |
| `loc_percentiles.csv` | Full percentile table (p5–p95, both scenarios) |
| `copula_sensitivity.csv` | Sensitivity table (Table 5.1 from memo) |

---

## Version History

| Version | Primary Correction | Remaining Limitation |
|---|---|---|
| v1.0–v6.0 | Circular inverse calculation; forward power-law estimation | Static volume; uniform density |
| v7.0 | Equation system; citations resolved | V = A·d geometry; homogeneous density |
| v8.0–v9.0 | Phase-space integration; heliocentric gradient r⁻¹·⁵; risk-weight recalibration | Solar model underdeveloped; point estimates |
| v9.1 | Solar model expanded; uncertainty ranges; tone moderated | No Monte Carlo; no Bayesian framing |
| v9.2 | Monte Carlo (N=100k); Gaussian copula; 18 references | No copula sensitivity; priors not documented |
| **v9.3** | Copula sensitivity (Gaussian vs. t, ρ range); full prior derivation appendix; reproducible code | No FMEA/fault tree; priors not empirically calibrated |

---

## References

- Cucinotta, F.A. et al. (2014). Space radiation cancer risk projections and uncertainties. *NASA/TP-2011-216155*.
- Fechtig, H. & Grün, E. (1981). Interplanetary dust in the solar system. *Space Research*, 21.
- Fragola, J.R. et al. (1996). Crew survival on board Space Station. *Probabilistic Safety Assessment and Management*.
- Grøn, E. et al. (1985). Collisional balance of the meteoritic complex. *Icarus*, 62.
- Hassler, D.M. et al. (2013). Mars' surface radiation environment measured with the Mars Science Laboratory's Curiosity Rover. *Science*, 343.
- Hendrickx, L. et al. (2006). Microbial ecology of the closed artificial ecosystem MELiSSA. *Res. Microbiol.*, 157.
- NASA (2011). Probabilistic Risk Assessment Procedures Guide for NASA Managers and Practitioners. *NASA/SP-2011-3421*.
- Nelsen, R.B. (2006). *An Introduction to Copulas* (2nd ed.). Springer.
- Opik, E.J. (1951). Collision probabilities with the planets. *Proc. R. Irish Acad.*, 54A.
- Schwenn, R. (2006). Space weather: The solar perspective. *Living Reviews in Solar Physics*, 3.
- Wetherill, G.W. (1967). Collisions in the asteroid belt. *J. Geophys. Res.*, 72.
- Zeitlin, C. et al. (2004). Energetic particle radiation in deep space: Results from the Marie instrument. *Adv. Space Res.*, 33.

---

## Repository Contents

```
mars-pra/
├── mars_pra_v93.py                   # Full Monte Carlo simulation (N=100,000)
├── Mars_Risk_v9_3_Quantix_Final.docx # Companion technical memorandum
├── requirements.txt                  # Python dependencies
├── CHANGELOG.md                      # Version history with corrections
├── .gitignore                        # Python / output ignores
└── output/                           # Generated figures and CSV outputs (gitignored)
```

---

## License

MIT — see [LICENSE](LICENSE) file.  
**Author:** Patrick Jauslin · Quantix Analytics Consulting · [quantixac.com](https://quantixac.com)

---

## Disclaimer

This is an **independent exploratory** PRA framework for external technical discussion. It is not:
- Engineering certification or flight safety documentation
- A peer-reviewed publication
- Affiliated with SpaceX, NASA, or any space agency

Prior distributions are expert-elicitation estimates anchored on analogue data. No Mars-class likelihood function exists for the dominant risk factors. All quantitative outputs carry substantial epistemic uncertainty and should be read as structured plausibility estimates.
