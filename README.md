# infoprocap

A MATLAB toolkit for computing and visualizing **information-processing capacity (IPC)** of dynamical systems using a multivariate Legendre polynomial basis.

All classes are accessed via the `infoprocap` namespace (e.g., `infoprocap.IPC`).

## Requirements

- MATLAB R2019b or later — no external toolboxes needed

## Package Structure

```
+infoprocap/
├── IPC.m       % Core class — basis construction & capacity computation
├── Plotter.m   % Static plotting utilities
└── Utils.m     % Math and progress helper functions
```

---

## `IPC` — Core Class

Constructs a multivariate orthonormal Legendre polynomial basis and computes IPC from reservoir readout signals.

### Constructor

```matlab
obj = infoprocap.IPC(u, max_deg)
```

- `u` — `[N × d]` input matrix (values in `[-1, 1]`)
- `max_deg` — maximum polynomial degree

### Key Properties

| Property | Description |
|---|---|
| `u` | Input signal `[N × d]` |
| `y` | Evaluated product basis terms `[N × basis_size]` |
| `degrees` | Multi-index degree table `[basis_size × d]` |
| `max_deg` / `basis_size` / `dimn` | Degree, number of basis terms, input dimension |

### Methods

**`C = calcCap(obj, X, sample_idxs, basis_idxs, use_bias)`**
Raw capacity from readouts `X`. Pass `0` to `sample_idxs` or `basis_idxs` to use all.

**`C = fitCap(obj, X, sample_idxs, alg)`**
Bias-corrected capacity with negative values zeroed.
- `alg = 1` — threshold by minimum-negative value
- `alg = 2` — theoretical threshold (call `initThresholds` first)

**`[Chat, dChat] = estCap(obj, X, alg)`**
Capacity with split-half uncertainty estimate (`dChat`).

**`[Cm_arr, Cs_arr] = scanCap(obj, X)`**
Scans capacity vs. number of training samples; returns mean and std across independent subsets.

**`initThresholds(obj, X, sample_idxs)`**
Pre-computes theoretical thresholds per basis term (required for `alg = 2`).

---

## `Plotter` — Visualization

Static class. Both methods accept a `filename` argument:

| Value | Behaviour |
|---|---|
| `"no_plot"` | Skip plotting |
| `"no_save"` | Display without saving |
| Any string | Save to file at 300 DPI |

**`Cm = Plotter.Cap_mat(obj, C, K, filename)`**
2D heatmap of capacity over polynomial degrees. Requires `d = 2`.

**`Cd = Plotter.Cap_deg(obj, C, K, cap_lim, filename)`**
Stacked bar chart of capacity by total degree, split into single-feature, 2-feature cross, and higher-order interaction terms.

---

## `Utils` — Helpers

**`Utils.dispPerc(i, len)`** — Prints progress percentage every 1% increment.

**`Z = Utils.stars_and_bars(n, k)`** — Enumerates all degree multi-indices via stars-and-bars combinatorics. Used internally by `IPC`.

---

## Example

```matlab
% 1. Build IPC object
u = 2*rand(500, 2) - 1;           % [N x d] inputs in [-1, 1]
ipc = infoprocap.IPC(u, 5);       % max degree 5

% 2. Compute bias-corrected capacity from reservoir readouts X
C = ipc.fitCap(X, 0, 1);
disp("Total capacity: " + sum(C));

% 3. Visualize
infoprocap.Plotter.Cap_mat(ipc, C, ipc.basis_size, "cap_matrix.png");
infoprocap.Plotter.Cap_deg(ipc, C, ipc.basis_size, 25, "cap_degree.png");
```
