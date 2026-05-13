# infoprocap

![IPC Visualisations](./images/visualisations.jpg)
*Figure 1: IPC Visualizations

A MATLAB toolkit for computing **Information Processing Capacity (IPC)** of stationary physical systems with limited experimental data using a multivariate Legendre polynomial basis.

All classes are accessed via the `infoprocap` namespace (e.g., `infoprocap.IPC`).

## Installation

Copy the `+infoprocap` folder into your MATLAB project's root directory.

## Requirements

- MATLAB R2019b or later — no external toolboxes needed

## Package Structure

```
+infoprocap/
├── IPC.m       % Core class — basis construction & capacity computation
├── Plotter.m   % Static plotting utilities
└── Utils.m     % Helper functions
```

---

## `IPC` — Core Class

Constructs a multivariate orthonormal Legendre polynomial basis and computes IPC from readouts.

### Constructor


- `u` —  input matrix (values in `[-1, 1]`)
- `max_deg` — maximum total degree of the Legendre polynomial product basis

### Key Properties

| Property | Description |
|---|---|
| `sample_size` / `basis_size` / `dimn` | No. of input samples, No. of basis terms, Input dimension |
| `u` | Input signal `[sample_size × dimn]` |
| `y` | Evaluated product basis terms `[sample_size × basis_size]` |
| `degrees` | Multi-index degree table `[basis_size × dimn]` (Used to generate product basis from individual polynomials) |


### Methods

**`C = calcCap(obj, X, sample_idxs, basis_idxs, use_bias)`**
Raw capacity from readouts `X`

**`C = fitCap(obj, X, sample_idxs, alg)`**
Fitted capacity after asymptotic fitting and false-positive thresholding.
- `alg = 1` — theoretical threshold
- `alg = 2` — threshold by minimum-negative value 

**`[C_hat, dC_hat] = estCap(obj, X, alg)`**
Fitted Capacity with split-half uncertainty estimate (`dC_hat`).

**`[samps_arr,Cm_arr] = scanCap(obj, X)`**
Scans raw capacity vs. number of samples. Useful to visualise the asymptotic form of capacities.

**`initThresholds(obj, X, sample_idxs)`**
Pre-computes theoretical thresholds per basis term (required for `alg = 1`).

---

## `Plotter` — Visualization

Static class. Both methods accept a `filename` argument:

| Value | Behaviour |
|---|---|
| `"no_plot"` | Skip plotting |
| `"no_save"` | Display without saving |
| Any string | Save to file at 300 DPI |

**`Cm = Plotter.capMat(obj, C, filename)`**
2D heatmap of capacity over polynomial degrees. Requires `dimn = 2`. [Figure 1(a)]

**`Cd = Plotter.capBar(obj, C, filename,y_lim)`**
Stacked bar chart of capacity by total degree, split into single-feature, 2-feature cross, and higher-order interaction terms.  [Figure 1(b)]

---

## `Utils` — Helpers

**`Utils.dispPerc(i, len)`** — Prints progress percentage every 1% increment.

**`Utils.stars_and_bars(n, k)`** — Enumerates all degree multi-indices via stars-and-bars combinatorics. Used internally by `IPC` to generate product basis.

**`Utils.W3j(j1,j3)`** — Wigner 3j function when j1=j2 and m1=m2=m3=0.

---

## Example
```matlab
ipc1 = infoprocap.IPC(u,8);              % initialize IPC object with input u and maximum degree of basis=8
[C,dC] = infoprocap.estCap(X,1);        % estimate capacities C with readouts X and Algorithm 1. dC is an approximate uncertainty measure of capacities
infoprocap.Plotter.capBar(ipc1,C);      % plot capacity bar plot
```
See example_photonic_elm.m for a basic example with a photonic system.

---

## References

[1] Joni Dambre, David Verstraeten, Benjamin Schrauwen, and Serge Massar. ''[Information processing capacity of dynamical systems.](https://www.nature.com/articles/srep00514)'' Scientific reports 2.1 (2012): 1-7.
