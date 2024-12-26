This readme file was generated on 2024-12-26 by Jordan Hachtel

If these data or notebooks have been helpful to you please cite our published manuscript. Thank you.

Title: Main Text Figures and Data for 'Hydrogen bonding in water under extreme confinemnt.'

Contact Information: Jordan A. Hachtel (hachtelja@ornl.gov)

Link to Published Manuscript (To be updated upon final publication)

Link to arXiv Paper: https://arxiv.org/abs/2402.17989

File List: Listed in relevant groupings
Jupyter Notebook used to process and visualize the data for Figure 1, 2 and 3 of the main text in the manuscript. 
	-Main_Text_Figure_1.ipynb
	-Main_Text_Figure_2.ipynb
	-Main_Text_Figure_3.ipynb 

Library of python functions is used in jupyter notebooks.
	-CNT_Analysis_Functions.py

vEEL Point Spectrum Datasets. These spectra are preprocessed and uploaded as numpy files with calibrated energy axes, hence the identical name with _E and _S suffixes.
	-EmptyCNT_RT_E.npy, EmptyCNT_RT_S.npy
	-FilledCNT_1pt4nm_Cryo_E.npy, FilledCNT_1pt4nm_Cryo_S.npy
	-FilledCNT_1pt4nm_RT_E.npy, FilledCNT_1pt4nm_RT_S.npy
	-FilledCNT_2pt3nm_Cryo_E.npy, FilledCNT_2pt3nm_Cryo_S.npy
	-FilledCNT_2pt3nm_RT_E.npy, FilledCNT_2pt3nm_RT_S.npy
	-FilledCNT_2pt3nm_RT_PostCryo_E.npy, FilledCNT_2pt3nm_RT_PostCryo_S.npy
	-LiquidCell_RT_E.npy, LiquidCell_RT_S.npy

vEEL Spectrum Image Datasets. These datasets are presented as acquired from the Nion Swift software exported as .npy datafile with metadata in a .json file of the same name. Processing performed in the Python Notebooks. Each dataset has a fully hyperspectral 3D EELS dataset and a corresponding simultaneously acquired HAADF image with at _Z suffix exported in the same manner.
	-EmptyCNT_Cryo_SI.npy, EmptyCNT_Cryo_SI.json
	-EmptyCNT_Cryo_SI_Z.npy, EmptyCNT_Cryo_SI_Z.json
	-FilledCNT_2pt3nm_Cryo_SI.npy, FilledCNT_2pt3nm_Cryo_SI.json
	-FilledCNT_2pt3nm_Cryo_SI_Z.npy, FilledCNT_2pt3nm_Cryo_SI_Z.json

Images. TEM and STEM Reference Images. One TEM image is only presented as a .png, all other images are presented as .npy datafiles with .json metadata.
	-CNT1_TEM.npy, CNT1_TEM.json
	-CNT2_TEM.npy, CNT2_TEM.json
	-SmallFOV_STEM.npy, SmallFOV_STEM.json
	-LargeFOV_TEM.png

Schematics. Schematics representing critical aspects of experiments or calculations. All exported as .png files.
	-EELS_Schematic.png
	-Molecular_Schematic.png

DFT vDOS Datasets. Vibrational density of states obtained from molecular dynamics run. Exported as .csv files with three columns Frequency (eV), Frequency (meV), vDOS. Processing conducted within Python Notebooks. vCNT data presented with room temperature (columns 1-3) cryo temperature (columns 4-6).
	-bulk_300K_1pt0_vDOS.csv
	-MixedPhaseIce_vDOS.csv
	-rCNT_300K_0pt5_vDOS.csv
	-vCNT_AllT_0pt5_vDOS.csv
	-vCNT_AllT_pt75_vDOS.csv
	-vCNT_AllT_1pt0_vDOS.csv
DFT Snapshot Datasets. Single snapshots from the MD simulations exemplifying structure of water under different conditions. Files are loaded into VESTA as .cif files and exported as .tif files for inclusion in the main text figures. Here we provide the exported .tif files. 
	-bulk_300K_1pt0_snapshot.tif
	-MixedPhaseIce_snapshot.tif
	-rCNT_300K_0pt5_snapshot.tif
	-vCNT_300K_0pt5_snapshot.tif

DFT Heatmap Datasets (Preprocessed). Statistical heatmaps of molecular configurations throughout MD runs. Main text figures use data processed in Origin and exported as .png files. We present the processed .png files and the raw data. Colorbar from Origin also included as .png file.
	-bulk_300K_1pt0_heatmap.png
	-rCNT_300K_0pt5_heatmap.png
	-vCNT_300K_0pt5_heatmap.png
	-Origin_ColorBar.png

DFT Heatmap Datasets (Unprocessed). Output of molecular configurations. Used for additional processing in Python Notebooks. Exported as .csv files with four columns containing oxygen-oxygen distances (dOO), oxygen-oxygen/oxygen-hydrogen bond angles (theta), oxygen-hydrogen distances (dOH), intermolecular hydrogen-oxygen distances (dHO).
	-bulk_300K_1pt0_heatmap.csv
	-rCNT_300K_0pt5_heatmap.csv
	-rCNT_300K_pt75_heatmap.csv
	-rCNT_300K_1pt0_heatmap.csv
	-vCNT_300K_0pt5_heatmap.csv
	-vCNT_300K_pt75_heatmap.csv
	-vCNT_300K_1pt0_heatmap.csv
	-vCNT_100K_0pt5_heatmap.csv
	-vCNT_100K_pt75_heatmap.csv
	-vCNT_100K_1pt0_heatmap.csv

Relationship between files: CNT_Analysis_Functions.py must be imported into all the Python Notebooks for them to run. All datasets and images are loaded in the Notebooks and processed or visualized.

Additional data and Python Notebooks are used for Extended Data and Supplementary Information in the manuscript. The authors will provide these on reasonable request. Methodological information and experimental conditions are described in detail in the Methods section of the manuscript. Processing methodologies are highlighted in the Python Notebook.
