# Solenoid Optimizer for MR Microscopy

## Overview

This MATLAB package optimizes micro-solenoid coil parameters for high-resolution magnetic resonance (MR) microscopy. It computes field homogeneity, coil resistance, inductance, capacitive and dielectric losses, and evaluates the Signal-to-Noise Ratio (SNR) for varying coil geometries. The tool is based on analytical models and empirical enhancement factor tables from Minard & Wind (2001).

## Features

- **Field Homogeneity Modeling**: Calculates `ΔB_xy/B_xy` within the sample volume using Biot-Savart law.
- **SNR Optimization**: Identifies optimal wire diameter and number of turns for maximum SNR.
- **Loss Analysis**: Models resistive losses from coil, leads, sample magnetic/dielectric interactions, and capacitive losses.
- **Temperature Compensation**: Computes SNR with and without temperature effects.
- **Graphical Visualization**: Field deviation maps and SNR-vs-geometry plots.

## File Structure
.
├── coilOptimizer.m                  # Master script: runs entire optimization, plots results
├── solenoidOptimizer.m           # Class definition for solenoid coil modeling and SNR computation
├── NturnOptimization.m           # Function to optimize number of turns and wire diameter for SNR
├── RFfield_map.m                 # Function to calculate and visualize Bxy field uniformity (sagittal slice)
├── deltaBxyCalc (inline)         # Local helper function inside main_script for ΔBxy/Bxy calculation


## Dependencies

- MATLAB R2018 or newer (recommended)
- No external toolboxes required

## How to Use

1. **Edit Sample Parameters**  
   Modify `lsample`, `dsample`, `capthickness`, etc. in `main_script.m` to reflect your sample setup.

2. **Run Optimization**  
   Execute `main_script.m`. It will:
   - Find the minimum coil length-to-diameter ratio for desired field homogeneity.
   - Optimize the number of turns and wire diameter for maximum SNR.
   - Plot SNR vs wire diameter and field deviation maps.

3. **Interpret Outputs**
   - **Console Output**: Displays best coil geometry (`n`, `d_wire`, `l_coil`) and SNR.
   - **Figure 1**: Plot of field inhomogeneity vs. `l_coil/d_coil`.
   - **Figure 2**: SNR vs wire diameter for multiple `n`.
   - **Figure 3**: B-field deviation map in the sagittal plane (YZ plane).


## Citation

If you use this tool in your research, please cite:

> Minard, K. R., & Wind, R. A. (2001). *Solenoidal microcoil design. Part II: Optimizing winding parameters for maximum signal-to-noise performance.* Concepts in Magnetic Resonance, 13(3), 190–210.

## Author

Bibek Dhakal  
Ph.D. Candidate, Department of Physics  
Vanderbilt University Institute of Imaging Science (VUIIS)  
Email: bibek.dhakal@vanderbilt.edu

Special thanks to Ben M. Hardy (Remcom Inc. and VUIIS Alumni) for the initial version of the solenoidOptimizer code.

---

**Last updated**: March 26, 2025
