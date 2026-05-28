# Changelog — Mars Mission PRA

All notable changes to this probabilistic risk assessment framework are documented here.

---

## v9.3 — May 2026 *(current)*

**Added**
- Copula sensitivity analysis: Gaussian vs. t-copula (ν=4) across ρ ∈ [0.0, 0.5]
- Full prior parameter derivation appendix (Appendix A) with analogue data anchors
- Reproducible Python code (`mars_pra_v93.py`) with CLI flags `--seed`, `--n`, `--out`
- CSV output files for all key tables (percentiles, copula sensitivity)
- Epistemic status matrix (Figure 9)

**Confirmed**
- p50 LOC is robust to copula family and ρ magnitude (varies < 0.5 pp)
- Dominant uncertainty source is prior distributions, not correlation structure

**Remaining limitations**
- No FMEA / fault tree
- Priors not empirically calibrated against Mars class operational records

---

## v9.2 — April 2026

**Added**
- Full Monte Carlo simulation (N = 100,000)
- Gaussian copula implementation (ρ = 0.30 for Engineering ↔ ECLSS)
- 18 literature references

**Limitations at release**
- No copula sensitivity analysis
- Prior distributions not documented in main text

---

## v9.1 — March 2026

**Added**
- Expanded solar / radiation model with scenario-based shielding estimates
- Explicit uncertainty ranges for all point estimates
- Moderated language to reflect exploratory status

**Limitations at release**
- No Monte Carlo; point estimates only
- No Bayesian framing

---

## v9.0 — February 2026

**Added**
- Phase-space flux integration for debris model
- Heliocentric radial gradient r⁻¹·⁵ (Grün 1985; Fechtig 1981)
- Risk-weight recalibration across all factors

**Limitations at release**
- Solar model underdeveloped
- Point estimates only; no uncertainty propagation

---

## v8.0 — January 2026

**Added**
- Phase-space integration replaces prior static volume approach
- Initial heliocentric gradient implementation

---

## v7.0 — December 2025

**Fixed**
- Corrected volume geometry (V = A·d; homogeneous density assumption documented)
- All citations resolved and verified

**Limitations at release**
- V = A·d geometry remains a simplification
- Homogeneous density assumption not validated

---

## v1.0–v6.0 — October–November 2025

**Initial versions**
- Circular inverse calculation identified and corrected
- Forward power-law estimation introduced

**Limitations**
- Static volume assumption
- Uniform density without heliocentric gradient
- Circular reference in original collision probability derivation

---

*Quantix Analytics Consulting · Patrick Jauslin · quantixac.com*
