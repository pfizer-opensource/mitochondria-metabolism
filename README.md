# In Silico Mitochondria Muscle Metabolism QSP Model

## Description
-----------
 This is a mechanistic model of myocardial mitochondrial metabolism including oxidative phosphorylation, $\beta$-oxidation, TCA cycle, metaboilte transport and electrophysiology. The model is an extension of previously published model code from  Wu et al., J Physiol 586 (17), 4193 (2008) and van Eunen et al., PLoS computational biology vol. 9,8 (2013). The model is used to evaluate how substrate selection influences myocardial energetics at rest and during exercise under healthy and heart failure conditions.

   ## Primary Results  
   -----

This code reproduces figures from *Assessing Potential Mitochondrial Targets in the Failing Heart Using Quantitative Systems Pharmacology*
 

### Prerequisites
---
MATLAB

This code was written in Matlab (v. 2019b)

## Setup
---
Add all files/directories in this repository to the MATLAB working directory/path.

## Contents
---
### *in vivo* Folder
1. `Run_Final_Figurs.m` -> Script to generate figures from the manuscript
2. `setup.m` -> Defines constants and initial conditions
3. `define_globals.m` -> Defines global variables
4. `Cell_dxdt.m` -> Diff Eqs. for cell version of the model
5. `Cell_Flux.m` -> Calculates Flux 
6. `PDHphosdephos.m` -> Calculates phosphorylation status of Pyruvate dehydrogenase
7.  `param14.mat` -> Parameter file
8. `Varstruc_default.mat` -> Variables that change upon pertubations for heart failuree and interventions
9. `Run_substrate_select.m` -> Function to sweep malonyl-coA and Pyruvate
10. `ATP_Hyrdo_Sweep.m` -> Function to increase workrate in mitochodria by increasing ATP hydrolysis
11. `Genreate_Fig4b.m` -> Function to reproduce figure 4 of the  manuscript
12. `Fig6B.m` -> Function to reproduce figur 6B of the manuscript

  ### *in vitro* Folder
  1. `setup1.m` -> Defines constants and initial conditions for optimization experiment
  2. `define_global1.m` -> Defines global variables
  3. `define_global_opt.m`   - > Defines other global variable options
  4. `IsoMito_dXdT`-> Diff Eqs. for Isolated Mitocondria
  5. `Mito_dXdT.m` -> Diff Eqs. for Mitochondria (similar to IsoMito) for fitting specifity factors in vitro (sf)
  6. `Mito_Flux.m` -> Calculates flux values
  7. `param14.mat` -> Parameter values
  8. `Reproduce_InVitro` -> Function to reproduce *in vitro* experiments
  9. `InVitroExpCond.m` -> Function to define *in vitro* experimental conditions
  10. `Run_sf_optimization.m` -> Function to optimize specificity factors in $\beta$-oxidation enzymes
  11. `Objfun.m` -> Function to define the objective function for the minimization problem
  12. `Plot_variability` -> Plotting functioon
  13. 'dxdt_electrode` -> Function to mimic O2 sensor delay in MVO2 experiments
  14. `Varstruc_default.mat` -> Variables that change upon pertubations for heart failuree and interventions
  15.  `s_F.mat` ->  Specificity factors
  16. `sf_out.mat` -> Fitted specificity factors
  17. `fluxes.m` -> Calculates flux for VO2 for IsoMito
  18.  `Table1.mat` -> Digitized data

## Usage
---
To generate figures from the manuscript:

`Run_Final_Figures.m `

Figures are saved in pdf format

## Authors
---
Lyndsey F. Meyer*, Neda Nourabadi, CJ Musante, Daniel A. Beard, Anna Sher

*Correspondence to Lyndsey.Meyer@pfizer.com or Dabeard@med.umich.edu
[![DOI](https://zenodo.org/badge/558893205.svg)](https://zenodo.org/badge/latestdoi/558893205)


![alt text](https://github.com/openPfizer/DigitalHealthData/blob/master/img/osbypfizer.png)
