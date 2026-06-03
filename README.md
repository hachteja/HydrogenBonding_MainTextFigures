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
pandas
jupyter
```

Install with:

```bash
pip install -r requirements.txt
```

or with conda:

```bash
conda install numpy scipy matplotlib pillow pandas jupyter
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
│   ├── vEEL_Spectrum_Images/       # Vibrational EELS spectrum image datasets
│   ├── Images/                     # TEM and STEM reference images
│   ├── Schematics/                 # Schematic diagrams (PNG)
│   ├── DFT_vDOS/                   # DFT vibrational density of states (CSV)
│   ├── DFT_Snapshots/              # MD simulation snapshots (TIFF)
│   └── DFT_Heatmaps/               # MD statistical heatmaps (CSV and PNG)
│       └── Preprocessed/
├── requirements.txt
└── README.md
```

---

## Data Files

### vEEL Point Spectra (`Data/vEEL_Point_Spectra/`)

Preprocessed vibrational EELS point spectra with calibrated energy axes. Each sample has a paired `_E.npy` (energy axis in eV) and `_S.npy` (signal) file.

| File prefix | Description |
|---|---|
| `EmptyCNT_RT` | Empty CNT, room temperature |
| `FilledCNT_1pt4nm_RT` | 1.4 nm water-filled CNT, room temperature |
| `FilledCNT_2pt3nm_RT` | 2.3 nm water-filled CNT, room temperature |
| `FilledCNT_2pt3nm_Cryo` | 2.3 nm water-filled CNT, cryogenic |
| `FilledCNT_2pt3nm_RT_PostCryo` | 2.3 nm water-filled CNT, room temperature after cryo cycle |
| `LiquidCell_RT` | Liquid cell reference, room temperature |

### vEEL Spectrum Images (`Data/vEEL_Spectrum_Images/`)

Nion Swift format `.npy` + `.json` pairs (3D EELS spectrum images and HAADF images) plus preprocessed 2D arrays for the 1.4 nm CNT datasets.

| File | Format | Description |
|---|---|---|
| `EmptyCNT_Cryo_SI` | `.npy` + `.json` | Empty CNT 3D EELS spectrum image, cryogenic |
| `EmptyCNT_Cryo_SI_Z` | `.npy` + `.json` | Simultaneously acquired HAADF image |
| `FilledCNT_2pt3nm_Cryo_SI` | `.npy` + `.json` | 2.3 nm filled CNT 3D EELS spectrum image, cryogenic |
| `FilledCNT_2pt3nm_Cryo_SI_Z` | `.npy` + `.json` | Simultaneously acquired HAADF image |
| `FilledCNT_1pt4nm_Cryo_2DSI_E` | `.npy` | 1.4 nm filled CNT cryo 2D-SI energy axis (eV) |
| `FilledCNT_1pt4nm_Cryo_2DSI_SI` | `.npy` | 1.4 nm filled CNT cryo 2D-SI spectrum image |
| `FilledCNT_pt8nm_Cryo_SISeq_E` | `.npy` | 0.8 nm filled CNT cryo SI sequence energy axis (eV) |
| `FilledCNT_pt8nm_Cryo_SISeq_SI` | `.npy` | 0.8 nm filled CNT cryo SI sequence |
| `FilledCNT_pt8nm_Cryo_SISeq_MAADF` | `.npy` | Simultaneously acquired MAADF image |

### Images (`Data/Images/`)

TEM and STEM reference images. `.npy`/`.json` pairs use Nion Swift format; one file is provided as `.png`.

| File | Description |
|---|---|
| `CNT1_TEM.npy/.json` | TEM image of CNT 1 |
| `CNT2_TEM.npy/.json` | TEM image of CNT 2 |
| `SmallFOV_STEM.npy/.json` | Small field-of-view STEM image |
| `LargeFOV_TEM.png` | Large field-of-view TEM image |

### Schematics (`Data/Schematics/`)

| File | Description |
|---|---|
| `EELS_Schematic.png` | Schematic of the EELS experiment |
| `Molecular_Schematic.png` | Schematic of molecular hydrogen bonding |
| `H2OCNTschematic.tiff` | Water-in-CNT hydrogen bonding schematic (Figure 2) |

### DFT vDOS (`Data/DFT_vDOS/`)

Vibrational density of states from molecular dynamics simulations. CSV files with three columns: Frequency (eV), Frequency (meV), vDOS.

| File | Description |
|---|---|
| `bulk_300K_1pt0_vDOS.csv` | Bulk water, 300 K, 1.0 g/cc |
| `MixedPhaseIce_vDOS.csv` | Mixed-phase ice |
| `rCNT_300K_0pt5_vDOS.csv` | Rigid CNT, 300 K, 0.5 g/cc |
| `vCNT_AllT_0pt5_vDOS.csv` | Vibrating CNT, all temperatures, 0.5 g/cc |
| `vCNT_AllT_pt75_vDOS.csv` | Vibrating CNT, all temperatures, 0.75 g/cc |
| `vCNT_AllT_1pt0_vDOS.csv` | Vibrating CNT, all temperatures, 1.0 g/cc |

### DFT Snapshots (`Data/DFT_Snapshots/`)

Single-frame MD simulation snapshots exported from VESTA as TIFF files.

| File | Description |
|---|---|
| `bulk_300K_1pt0_snapshot.tif` | Bulk water, 300 K, 1.0 g/cc |
| `MixedPhaseIce_snapshot.tif` | Mixed-phase ice |
| `rCNT_300K_0pt5_snapshot.tif` | Rigid CNT, 300 K, 0.5 g/cc |
| `vCNT_300K_0pt5_snapshot.tif` | Vibrating CNT, 300 K, 0.5 g/cc |

### DFT Heatmaps (`Data/DFT_Heatmaps/`)

Statistical heatmaps of molecular configurations throughout MD runs. Raw CSV files have four columns: O-O distance (dOO), bond angle (theta), O-H distance (dOH), intermolecular H-O distance (dHO). Preprocessed PNG files were generated in Origin.

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
