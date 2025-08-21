# microSolenoidGeometryOptimizer-v1.2

## About
The MATLAB class file **solenoidOptimizer.m** optimizes the geometry of a micro-solenoid to provide the maximum SNR for Magnetic resonance microscopy expeirments. It takes various input of a micro-solenoid coil geometry such as: number of turns, coil diameter, wire diameter, conductivity of the wire,  and length of the lead wires.  The other inputs are the sample length, sample diameter, temperature of the coil, temperature of the sample, and operating frequency of the NMR spectrometer. The **NturnOptimization.m** utilizes the class file **solenoidOptimizer.m** to simulate the SNR of micro-solenoids with different number of turns and wire diameter while keeping the rest of the geometry parameter constant. The **RFfield_map.m** generates the B1 sensitivity map of the solenoid coil. The **coilOptimizer.m** utilizes the all of the above matlab files to calculate required coil length, generate optimal number of turns versus wire diameter plot by iterative numerical simulation, and lastly generating the B1 map. 
The final result is similar to Figure 12 in the Minard and Wind (2001) paper. The MATLAB class file **solenoidOptimizer.m** can be edited if a different circuit design apart from Minard and Wind's equivalent circuit is being used. Lastly, the current state of **solenoidOptimizer.m** performs the net loss calculation assuming a conductive sample. For simulating the SNR of the micro-solenoids with non-conductive sample, the **Rnmr** parameter should be changed to **Rnmr = (self.Rcoil + self.Rleads+self.Rcap)** include just the coil resistance, lead resistance, and the capacitive losses. Also, for simulating the effect of temperature while calculating SNR, **SNRTemp** shoule be used which uses **Rnmr = RnmrTemp** need to be used in the solenoidOptimizer. 

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
- **coilOptimizer.m**:               Master script: runs entire optimization, plots results.
- **solenoidOptimizer.m**:           Class definition for solenoid coil modeling and SNR computation.
- **NturnOptimization.m**:           Function to optimize number of turns and wire diameter for SNR.
- **RFfield_map.m**:                 Function to calculate and visualize Bxy field uniformity (sagittal slice).

## Dependencies

- MATLAB R2018 or newer (recommended)
- No external toolboxes required

## How to Use

1. **Edit Sample Parameters**  
   Modify `lsample`, `dsample`, `capthickness`, etc. in `coilOptimizer.m` to reflect your sample setup.

2. **Run Optimization**  
   Execute `coilOptimizer.m`. It will:
   - Find the minimum coil length-to-diameter ratio for desired field homogeneity.
   - Optimize the number of turns and wire diameter for maximum SNR.
   - Plot SNR vs wire diameter and field deviation maps.

3. **Interpret Outputs**
   - **Console Output**: Displays best coil geometry (`n`, `d_wire`, `l_coil`) and SNR.
   - **Figure 1**: Plot of field inhomogeneity vs. `l_coil/d_coil`.
   - **Figure 2**: SNR vs wire diameter for multiple `n`.
   - **Figure 3**: B-field deviation map in the sagittal plane (YZ plane).

## Author

Bibek Dhakal  
Ph.D. Candidate, Department of Physics  
Vanderbilt University Institute of Imaging Science (VUIIS)  
Email: bibek.dhakal@vanderbilt.edu

Special thanks to Ben M. Hardy (Remcom Inc. and VUIIS + Vanderbilt Alumnus) for the initial version of the solenoidOptimizer code.

## Citation

If you use this tool in your research, please cite:
> Bibek Dhakal, Benjamin Hardy, Adam Anderson et al. Technical developments for high resolution magnetic resonance microscopy in a horizontal bore magnet, 10 May 2025, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-6550180/v1].
> Minard, K. R., & Wind, R. A. (2001). *Solenoidal microcoil design. Part II: Optimizing winding parameters for maximum signal-to-noise performance.* Concepts in Magnetic Resonance, 13(3), 190–210.


## Installation:
1. Download and unzip the **microSolenoidGeometryOptimizer** repository
2. Store it in a folder
3. Add the path of the folder to matlab path
4. open the **coilOptimizer.m** file in the Matlab
5. run and enjoy!


