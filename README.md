# microSolenoidGeometryOptimizer

## Overview
This MATLAB package provides tools to optimize micro-solenoid coil geometry for high-resolution magnetic resonance (MR) microscopy.  
It evaluates **field homogeneity**, **coil resistance and inductance**, **capacitive and dielectric losses**, and **signal-to-noise ratio (SNR)** across different coil geometries.  
The approach is based on analytical models and empirical enhancement factor tables from Minard & Wind (2001).

## Features
- **Field Homogeneity Modeling**: Calculates `ΔB_xy/B_xy` within the sample volume using Biot–Savart law.  
- **SNR Optimization**: Identifies optimal wire diameter and number of turns for maximum SNR.  
- **Loss Budgeting**: Includes resistive losses from the coil, lead wires, capacitors, and sample (magnetic + dielectric).  
- **Temperature Dependence**: Simulates SNR changes with effective temperature (`T_eff`) and supports cryogenic conditions.  
- **Visualization**: Generates SNR vs. geometry plots and B1-field deviation maps.

## File Structure
- **coilOptimizer.m**: Master script that runs the entire optimization and plots results.  
- **solenoidOptimizer.m**: Class definition for solenoid modeling and SNR computation.  
- **NturnOptimization.m**: Optimizes number of turns and wire diameter for SNR.  
- **RFfield_map.m**: Calculates and visualizes B1 sensitivity maps.  

## Dependencies
- MATLAB R2018 or newer (recommended)  
- No external toolboxes required  

## Usage
1. **Set Sample Parameters**  
   Edit `lsample`, `dsample`, `capthickness`, etc. in `coilOptimizer.m` to match your experiment.  

2. **Run Optimization**  
   Run `coilOptimizer.m`. The script will:  
   - Find the minimum coil length-to-diameter ratio for desired field homogeneity.  
   - Optimize number of turns and wire diameter for maximum SNR.  
   - Output plots of SNR vs. geometry and B-field deviation maps.  

3. **Outputs**  
   - **Console**: Best coil geometry (`n`, `d_wire`, `l_coil`) and estimated SNR.  
   - **Figure 1**: Field inhomogeneity vs. `l_coil/d_coil`.  
   - **Figure 2**: SNR vs. wire diameter for multiple turn counts.  
   - **Figure 3**: B1-field map (sagittal slice).  

## Installation
1. Download and unzip the repository.  
2. Add the folder to your MATLAB path.  
3. Open `coilOptimizer.m`.  
4. Run and view results.  

## Citation
If you use this tool, please cite:  
- Bibek Dhakal, Benjamin Hardy, Adam Anderson et al. *Technical developments for high-resolution magnetic resonance microscopy in a horizontal bore magnet*. Research Square (2025).  
- Minard, K. R., & Wind, R. A. (2001). *Solenoidal microcoil design. Part II: Optimizing winding parameters for maximum signal-to-noise performance.* *Concepts in Magnetic Resonance*, 13(3), 190–210.

## Author
**Bibek Dhakal**  
Ph.D. Candidate, Department of Physics  
Vanderbilt University Institute of Imaging Science (VUIIS)  
📧 bibek.dhakal@vanderbilt.edu  

Special thanks to **Ben M. Hardy** (Remcom Inc., VUIIS alumnus) for the initial version of the code.
