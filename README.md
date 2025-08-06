# χ<sup>zz</sup> - Longitudinal Susceptibility Solver for Real-Space Dyson Equation

![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository provides a Python implementation for computing the **longitudinal spin susceptibility** \( \chi^{zz} \) in real space using a Dyson-like self-consistent framework. The solver is designed to support multi-site systems such as BCC Chromium or Iron, with site- and spin-dependent exchange kernels.

## Overview

The code performs the following steps:

1. Parses raw input files (`.dat`) containing bare response values \( \chi^{0}_{\uparrow}, \chi^{0}_{\downarrow} \).
2. Constructs valid \( j \rightarrow k \) neighborhoods based on spatial geometry.
3. Solves the Dyson equation per site to compute the longitudinal response \( \chi^{zz}_j \).
4. Generates optional visualizations, including decay plots and 3D Fourier projections.
5. Validates symmetry by comparing \(\vec{r}\) vs \(-\vec{r}\) response.

---

## Repository Structure

```text
root/
├── input-file.dat                  # Required: raw χ⁰ data (with dx, dy, dz, χ⁰↑, χ⁰↓)
│
├── format.py                      # Step 1: formats input into CSV
├── determine_k_site.py           # Step 2: finds all valid k-sites within cutoff
├── solving_equation_mult.py      # Step 3: solves Dyson-like equations for χᶻᶻ
├── decay_plot.py                 # Step 4: generates static vs dynamic decay plots
├── thr_d_plot.py                 # Step 5: optional 3D or Fourier-space visualization
├── main.py                       # Driver script, configure and run everything
│
├── formated_data.csv             # Output: formatted χ⁰ data with shifts
├── k_pot_coords.csv              # Output: full list of k-vectors in cutoff radius
├── neighbouring_k_to_j.json      # Output: final k-mapping for Dyson solver
├── xzz_output_iter.csv           # Output: solved χᶻᶻ per j-site
├── formated_output.dat           # Output: merged χ⁰ and χᶻᶻ per site
├── debug_iter.txt                # Optional: iteration log
│
├── xzz_decay_plot.png            # Plot: χᶻᶻ spatial decay
├── static_decay_plot.png         # Plot: χ⁰ static decay
├── comparison_decay_plot.png     # Plot: overlay of χ⁰ vs χᶻᶻ decay
````

---

## Input Format

The main input file (`input-file.dat`) must follow the structure:

```
# i    j    dx    dy    dz    Jij        χ⁰↑           χ⁰↓
```

Each line corresponds to a pair of sites with:

* `i, j`: site indices (ignored in solver, used for sublattice logic)
* `dx, dy, dz`: integer lattice displacements between sites
* `χ⁰↑`, `χ⁰↓`: diagonal static response values (spin-up and spin-down)

---

## Recommended Inputs

* Ensure uniform and symmetric sampling over spatial shells
* Lattice positions should follow integer-based Cartesian grid units
* Sub-sublattice shifts (e.g. for AFM systems) must be handled via `shift_rules` in `main.py`

---

## Dependencies

Install required packages with:

```bash
pip install numpy pandas matplotlib tqdm
```

Tested with Python 3.8+

---

## Running the Solver

1. Edit `main.py` and set the correct input path:

```python
url = r"AFM-Cr\AFM-chfile.dat"  # or your own file
```

2. Choose the correct:

   * `base_change_matrix_*` depending on your crystal symmetry
   * `shift_rules_*` for sublattice configurations (e.g., AFM)

3. Run the solver:

```bash
python main.py
```

The script will:

* Format the input
* Determine k-site neighborhoods
* Solve for longitudinal χᶻᶻ
* Produce visualizations

---

## Visualization Options

Enable or disable plotting in `main.py` by modifying the relevant block near the end of the script. The plotting scripts support real-space decay and 3D projections.

---

## Notes

* Only real-valued, diagonal χ⁰ values are currently supported.
* The mapping and solver can handle different site types (multi-site logic).
* For Fourier-space comparisons, ensure proper normalization and consistent grid definitions.

---

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for more information.

