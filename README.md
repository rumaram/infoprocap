# infoprocap

![IPC Visualisations](./images/visualisations.jpg)
*Figure 1: Example IPC visualizations.*

A MATLAB package for computing **Information Processing Capacity (IPC)** of stationary memory-less systems from finite data using a multivariate orthonormal Legendre polynomial basis.

All classes are accessed through the `infoprocap` namespace, for example `infoprocap.IPC`.

## Installation

Copy the `+infoprocap` folder into the root of your MATLAB project.

## Requirements

- MATLAB R2019b or later
- No external toolboxes required

## Package structure

```text
+infoprocap/
├── IPC.m       % Core IPC class
├── Plotter.m   % Visualization helpers
└── Utils.m     % Utility functions
```

## IPC class

`IPC` constructs a multivariate orthonormal Legendre product basis from the input signal and provides methods for raw capacity computation, thresholded capacity estimation, sample-size scans, basis lookup, and CSV export.

### Constructor

```matlab
ipc = infoprocap.IPC(u, max_deg)
```

#### Inputs

- `u` — input matrix of size `[sample_size x dimn]`, with values in `[-1, 1]`
- `max_deg` — maximum total degree used to generate the polynomial product basis

### Properties

| Property | Description |
|---|---|
| `u` | Input samples |
| `u_basis` | Individual Legendre basis values for every sample, dimension, and degree |
| `y` | Evaluated multivariate product basis terms |
| `max_deg` | Maximum total polynomial degree |
| `degrees` | Multi-index degree table for the product basis |
| `basis_size` | Number of basis terms |
| `sample_size` | Number of input samples |
| `dimn` | Number of input dimensions |
| `threshs` | False-positive capacity thresholds |
| `K` | Number of readouts used in the most recent estimation or scan |
| `basis_names` | Human-readable names of basis terms |
| `disp_prog` | Flag for command-window progress display |

### Methods

#### `C = calcCap(readouts, sample_idxs, basis_idxs, use_bias)`

Computes raw capacities from a readout matrix.

**Arguments**
- `readouts` — readout matrix
- `sample_idxs` — optional subset of sample indices, default: all samples
- `basis_idxs` — optional subset of basis indices, default: all basis terms
- `use_bias` — optional flag to append a bias column of ones, default: `1`

#### `C = estCap(X, alg, sample_idxs)`

Estimates capacities using two-sample extrapolation followed by thresholding.

**Arguments**
- `X` — readout matrix
- `alg` — thresholding method, default: `1`
  - `1` — theoretical threshold
  - `2` — minimum-negative threshold
- `sample_idxs` — optional subset of sample indices used for estimation

#### `initThresholds(X)`

Computes theoretical false-positive thresholds for each basis term and stores them in `threshs`.

#### `[samps_arr, C_arr] = scanCap(X)`

Computes raw capacities while increasing the number of samples from `K + 1` up to `sample_size`.

#### `idx = nameToidx(name)`

Returns the basis index corresponding to a basis-term name.

#### `name = idxToname(idx)`

Returns the basis-term name corresponding to a basis index.

#### `exportCSV(C, filename)`

Exports IPC results to a CSV file containing a summary block, total capacity statistics, and one row per basis function with capacity and threshold values.

## Plotter class

`Plotter` is a static utility class for visualizing IPC results.

Both plotting methods accept a `filename` argument with the following behavior:

| Value | Behavior |
|---|---|
| `"no_plot"` | Skip plotting |
| `"no_save"` | Display the figure without saving |
| Any other string | Save the figure at 300 DPI |

### `Cm = Plotter.capMat(C, filename)`

Creates a 2D heatmap of capacity over polynomial degrees. Requires `dimn = 2`.

### `Cd = Plotter.capBar(C, filename, y_lim)`

Creates a stacked bar chart of capacity by total degree, split into single-feature, two-feature cross, and higher-order interaction terms.

## Utils class

`Utils` contains helper routines used by the package.

- `Utils.dispPerc(i, len)` — prints progress updates every 1%
- `Utils.stars_and_bars(n, k)` — enumerates degree multi-indices for basis construction
- `Utils.W3j(j1, j3)` — evaluates the required Wigner 3j term used in threshold calculations

## Example

```matlab
% Input samples in [-1, 1]
ipc = infoprocap.IPC(u, 8);

% Raw capacities
C_raw = ipc.calcCap(X);

% Estimated capacities with theoretical thresholding
C_est = ipc.estCap(X, 1);

% Basis lookup
idx = ipc.nameToidx("P_1(u_1)");
name = ipc.idxToname(idx);

% Export results
ipc.exportCSV(C_est, "ipc_results.csv");

% Plot capacities
infoprocap.Plotter.capBar(C_est, "no_save", []);
```

See `example_photonic_elm.m` for a basic photonic-system example.

## Reference

Joni Dambre, David Verstraeten, Benjamin Schrauwen, and Serge Massar, ["Information processing capacity of dynamical systems"](https://www.nature.com/articles/srep00514), *Scientific Reports* 2, 514 (2012).
