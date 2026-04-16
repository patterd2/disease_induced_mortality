# Disease-Induced Mortality

Phase-space analysis of an SIR epidemic model in which the per-infection mortality rate $f(I)$ decreases with infection prevalence — a signature of pathogen-mediated selection for reduced virulence.

---

## Model

$$S' = \Lambda - \beta S I - \mu S + \delta R$$

$$I' = \beta S I - (\gamma + \mu) I - f(I)\,I$$

$$R' = \gamma I - (\mu + \delta) R$$

where $\beta = \mathcal{R}_0 \,(\mu/\Lambda)\,(\gamma + \mu + b)$ and the disease-induced mortality function is

$$f(I) = b + \frac{a I^2}{1 + d I^2}, \qquad d = \frac{-a\,m}{(m-1)\,b}$$

For $a < 0$, $f$ decreases monotonically from $f(0) = b$ to $f(\infty) = b/m$.

---

## Scripts

| File | Purpose |
|------|---------|
| `single_trajectory.m` | Integrates the ODE system and plots 2D/3D phase portraits and time series |
| `plot_mortality_function.m` | Plots $f(I)$ (rational form) for varying values of $a$ |
| `plot_mortality_function_alt.m` | Plots $f(I)$ (logistic form) for varying inflection points $\theta$ |
| `set_ggplot_style.m` | Shared ggplot2-style MATLAB figure defaults and colour palette |

---

## Baseline parameters

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| $\Lambda$ | $1/2600$ day$^{-1}$ | Recruitment rate |
| $\mu$ | $1/2600$ day$^{-1}$ | Natural mortality rate |
| $\mathcal{R}_0$ | $1.0023$ | Basic reproduction number |
| $\gamma$ | $1$ day$^{-1}$ | Recovery rate |
| $\delta$ | $2$ day$^{-1}$ | Rate of waning immunity |
| $b$ | $0.02$ day$^{-1}$ | Baseline disease-induced mortality |
| $a$ | $-10^5$ day$^{-1}$ | Mortality curvature parameter |
| $m$ | $100$ | Fold reduction in mortality at high prevalence |

---

## Dependencies

MATLAB (no additional toolboxes required). `set_ggplot_style.m` is bundled in the repository and must be on the MATLAB path before running any plotting script.
