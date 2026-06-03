# Main Text Figures — Hydrogen Bonding in Water under Extreme Confinement

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20534179.svg)](https://doi.org/10.5281/zenodo.20534179)

**Contact:** Jordan A. Hachtel (hachtelja@ornl.gov)

**Manuscript:** *Hydrogen bonding in water under extreme confinement*
- arXiv preprint: https://arxiv.org/abs/2402.17989
- Published manuscript link: *(to be updated upon publication)*

**Data:** All datasets are archived on Zenodo: https://doi.org/10.5281/zenodo.20534179

---

## Overview

This repository contains the Python notebooks, data, and analysis code needed to reproduce Figures 1, 2, and 3 of the main text. Each notebook is self-contained and produces a publication-quality figure as output.

Large spectrum image datasets (>39 MB) are hosted on Zenodo and are downloaded automatically the first time each notebook is run.

---

## Setup

### Requirements

Python 3.8+ with the following packages:

```
numpy
scipy
matplotlib
Pillow
jupyterlab
```

Install with:

```bash
pip install -r requirements.txt
```

or with conda:

```bash
conda install numpy scipy matplotlib pillow jupyterlab
```

### Running the notebooks

1. Clone or download this repository.
2. Launch Jupyter from anywhere:
   ```bash
   jupyter notebook
   ```
3. Open a notebook from the `Notebooks/` folder and run all cells.

The notebooks automatically detect the repository root — no path configuration required.

---

## Repository Structure

```
HydrogenBonding_MainTextFigures/
├── Notebooks/                      # Jupyter notebooks (one per figure)
│   ├── Main_Text_Figure_1.ipynb
│   ├── Main_Text_Figure_2.ipynb
│   └── Main_Text_Figure_3.ipynb
├── Code/
│   └── CNT_Analysis_Functions.py   # Shared analysis function library
├── Data/
│   ├── vEEL_Point_Spectra/         # Preprocessed vibrational EELS point spectra
│   ├── vEEL_Spectrum_Images/       # Spectrum image datasets (see note below)
│   ├── Images/                     # TEM and STEM reference images
│   └── Schematics/                 # Schematic diagrams
├── requirements.txt
└── README.md
```

> **Data hosting:** Small files (< 1 MB) are stored directly in this repository. Large spectrum image datasets (39–150 MB each) are hosted on Zenodo and are downloaded automatically the first time the relevant notebook is run. All datasets are also available for direct download at https://doi.org/10.5281/zenodo.20534179.

---

## Data Files

The table below uses **GitHub** and **Zenodo** to indicate where each file is stored. GitHub files are present immediately after cloning. Zenodo files are downloaded automatically on the first notebook run.

### vEEL Point Spectra (`Data/vEEL_Point_Spectra/`)

Preprocessed vibrational EELS point spectra with calibrated energy axes. Each dataset has a paired `_E.npy` (energy axis in eV) and `_S.npy` (signal intensity) file. All files are in **GitHub**.

| File prefix | Description |
|---|---|
| `EmptyCNT_RT` | Empty CNT, room temperature |
| `FilledCNT_1pt4nm_RT` | 1.4 nm water-filled CNT, room temperature |
| `FilledCNT_2pt3nm_RT` | 2.3 nm water-filled CNT, room temperature |
| `FilledCNT_2pt3nm_Cryo` | 2.3 nm water-filled CNT, cryogenic |
| `FilledCNT_2pt3nm_RT_PostCryo` | 2.3 nm water-filled CNT, room temperature after cryo cycle |
| `LiquidCell_RT` | Liquid cell reference, room temperature |

### vEEL Spectrum Images (`Data/vEEL_Spectrum_Images/`)

Hyperspectral EELS datasets. Nion Swift `.npy` + `.json` pairs store 3D spectrum images with spatial and spectral calibration metadata. Preprocessed 1D energy-axis files (`_E.npy`) are provided separately for datasets that are not in Nion Swift format.

| File | Location | Description |
|---|---|---|
| `EmptyCNT_Cryo_SI.json` | GitHub | Nion Swift metadata for empty CNT cryo SI |
| `EmptyCNT_Cryo_SI.npy` | **Zenodo** (39 MB) | Empty CNT 3D EELS spectrum image, cryogenic |
| `EmptyCNT_Cryo_SI_Z.npy/.json` | GitHub | Simultaneously acquired HAADF image |
| `FilledCNT_2pt3nm_Cryo_SI.json` | GitHub | Nion Swift metadata for 2.3 nm CNT cryo SI |
| `FilledCNT_2pt3nm_Cryo_SI.npy` | **Zenodo** (49 MB) | 2.3 nm filled CNT 3D EELS spectrum image, cryogenic |
| `FilledCNT_2pt3nm_Cryo_SI_E.npy` | GitHub | 2.3 nm CNT cryo SI energy axis (eV) |
| `FilledCNT_2pt3nm_Cryo_SI_Z.npy/.json` | GitHub | Simultaneously acquired HAADF image |
| `FilledCNT_1pt4nm_Cryo_2DSI_E.npy` | GitHub | 1.4 nm CNT cryo 2D-SI energy axis (eV) |
| `FilledCNT_1pt4nm_Cryo_2DSI_SI.npy` | **Zenodo** (150 MB) | 1.4 nm filled CNT cryo 2D spectrum image |
| `FilledCNT_pt8nm_Cryo_SISeq_E.npy` | GitHub | 0.8 nm CNT cryo SI sequence energy axis (eV) |
| `FilledCNT_pt8nm_Cryo_SISeq_MAADF.npy` | GitHub | 0.8 nm CNT cryo simultaneously acquired MAADF image |
| `FilledCNT_pt8nm_Cryo_SISeq_SI.npy` | **Zenodo** (57 MB) | 0.8 nm filled CNT cryo spectrum image sequence |

### Images (`Data/Images/`)

TEM and STEM reference images in Nion Swift format (`.npy` + `.json`) or `.png`. All files are in **GitHub**.

| File | Description |
|---|---|
| `SmallFOV_STEM.npy/.json` | Small field-of-view HAADF STEM image |
| `LargeFOV_TEM.png` | Large field-of-view TEM image |

### Schematics (`Data/Schematics/`)

All files are in **GitHub**.

| File | Description |
|---|---|
| `EELS_Schematic.png` | Schematic of the vibrational EELS experiment (Figure 1) |
| `H2OCNTschematic.tiff` | Water-in-CNT hydrogen bonding schematic (Figure 2) |

---

## Analysis Code (`Code/CNT_Analysis_Functions.py`)

Shared Python library imported by all notebooks. Key functions:

| Function | Description |
|---|---|
| `LoadEELS_SI(fname)` | Load a Nion Swift EELS spectrum image (.npy + .json) |
| `LoadImage(fname)` | Load a Nion Swift STEM image (.npy + .json) |
| `GetSpectraEnergyAxes(dat, disp)` | ZLP centering / energy axis calibration |
| `GetCalibratedSpectra(dat, ens)` | Align spectra to a common energy axis |
| `GetExponentialFit_2R(en, spec, n, ...)` | Power-law background subtraction |
| `GetFilteredSpectrum(spec, w, order)` | Savitzky-Golay smoothing |
| `Get_i(arr, val)` | Find array index closest to a value |
| `NormArray(arr)` | Normalize array to 0–1 range |

---

## Notes

- Extended data and supplementary figures are available from the authors upon reasonable request.
- Methodological details and experimental conditions are described in the Methods section of the manuscript.
- The spectrum image `.npy` + `.json` format follows the Nion Swift software export convention. The JSON contains spatial and spectral calibrations under the key `spatial_calibrations` (list of dicts with `scale`, `offset`, and `units` fields; the last entry is always the spectral/energy axis).
