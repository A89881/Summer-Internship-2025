# χ<sup>zz</sup> — Longitudinal Susceptibility Solver for Real-Space Dyson Equation

![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository implements a general-purpose Python framework to compute the **longitudinal spin susceptibility** $\( \chi^{zz} \)$ in real space. The solver uses a Dyson-like equation to compute the response per site based on bare static susceptibilities. It is designed to be general for **any crystal system**, given that appropriate lattice transformations, sublattice shifts, and interaction parameters are defined.

---

## Overview

The pipeline consists of the following steps:

1. Format raw RSPt-style input files (bare response data).
2. Generate candidate and valid k-sites for each j-site.
3. Map j-sites to valid k-sites based on distance and lattice transform.
4. Solve a Dyson-like system for $\( \chi^{zz} \)$ site-by-site.
5. Optionally visualize results and verify symmetry.

The script supports systems with:
- Arbitrary Bravais lattices
- Multiple sublattices and site-dependent shifting
- Multi-site interaction kernels $\( U_{↑↑}, U_{↓↓}, U_{↑↓}, U_{↓↑} \)$

---

## Repository Structure

```text
root/
├── main.py                     # Main execution and configuration script
│
├── format.py                   # Step 1: format input file to CSV
├── determine_k_site.py         # Step 2: find valid k-sites per j
├── solving_equation_mult.py    # Step 3: Dyson equation solver for χᶻᶻ
├── decay_plot.py               # Step 4: decay plots
├── thr_d_plot.py               # Step 5: optional 3D/FT visualization
│
├── input-file.dat              # Raw RSPt output (user-provided)
│
├── formated_data.csv           # Output: parsed χ⁰ data
├── k_pot_coords.csv            # Output: candidate k-sites
├── neighbouring_k_to_j.json    # Output: final k → j mapping
├── xzz_output_iter.csv         # Output: solved χᶻᶻ per j-site
├── formated_output.dat         # Output: combined χ⁰ and χᶻᶻ
├── debug_iter.txt              # Optional debug log
│
├── static_decay_plot.png       # Plot: χ⁰ decay
├── xzz_decay_plot.png          # Plot: χᶻᶻ decay
├── comparison_decay_plot.png   # Overlay plot: χ⁰ vs χᶻᶻ
````

---

## Input Format

The main input file (`input-file.dat`) must be **space-separated** and include this header:

```
i    j    dx    dy    dz    Jij        χ⁰↑           χ⁰↓
```

Where:

* `i, j`: site indices (used for applying sublattice shift rules)
* `dx, dy, dz`: integer displacements from site `j` to `i`
* `χ⁰↑`, `χ⁰↓`: bare static spin susceptibility components
* `Jij`: optional and unused in the solver

Example line:

```
1  1  0  0  1  0.0000  0.01234  0.01234
```

---

## Customizing for Your System

The solver is fully general, allowing you to provide custom physical parameters for your system in `main.py`.

### 1. Choose or define your material dataset

```python
url = r"YourFolder/your-input.dat"
```

### 2. Provide a lattice transformation matrix

For example, for a cubic system:

```python
base_change_matrix = [
    [1.0, 0.0, 0.0],  e_x
    [0.0, 1.0, 0.0],  e_y
    [0.0, 0.0, 1.0],  e_z
]
```

Or define your own matrix that maps from integer basis vectors to Cartesian coordinates.

### 3. Define sublattice shift rules (if applicable)

```python
shift_rules = {
    (1, 2): (-0.5, -0.5, -0.5),
    (2, 1): (0.5, 0.5, 0.5)
}
```

Leave empty if the system has no sublattice logic:

```python
shift_rules = {}
```

### 4. Set your U-kernel parameters

Example (uniform across all sites):

```python
U_params = [(2.0, 2.0, 0.0, 0.0)]
```

Or define one per site type (multi-site kernel support):

```python
U_params = [
    (2.0, 2.0, 0.0, 0.0),  # Site type 1
    (2.5, 2.5, 0.0, 0.0),  # Site type 2
]
```

---

## Running the Solver

1. Configure `main.py` with:

   * Path to input file (`url`)
   * Lattice transformation matrix
   * Sublattice shift rules
   * Kernel parameters

2. Run:

```bash
python main.py
```

This will:

* Format the data
* Compute neighborhoods
* Solve for χᶻᶻ
* Optionally generate plots

---

## Output Files

| File                        | Purpose                              |
| --------------------------- | ------------------------------------ |
| `formated_data.csv`         | Parsed and shifted input             |
| `k_pot_coords.csv`          | All spatial candidates for k-sites   |
| `neighbouring_k_to_j.json`  | Mapping of j-sites to valid k-sites  |
| `xzz_output_iter.csv`       | Resulting χᶻᶻ values for all j-sites |
| `formated_output.dat`       | Combined input and output            |
| `debug_iter.txt`            | Optional log                         |
| `xzz_decay_plot.png`        | Decay of χᶻᶻ                         |
| `static_decay_plot.png`     | Decay of χ⁰                          |
| `comparison_decay_plot.png` | Overlay: χ⁰ vs χᶻᶻ                   |

---

## Dependencies

Use `pip` to install required packages:

```bash
pip install -r requirements.txt
```

**requirements.txt:**

```
numpy
pandas
matplotlib
tqdm
chardet
```

Tested with Python 3.8+

---

## Notes

* Supports general materials and lattice symmetries
* Accepts arbitrary sublattice shifting and kernel structure
* Input must be symmetric and consistent for physical accuracy
* Output visualizations help verify decay and symmetry

---

## License

MIT License. See [LICENSE](LICENSE) for terms.
