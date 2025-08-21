# microSolenoidGeometryOptimizer-v1.2

## About
The MATLAB class file **solenoidOptimizer.m** optimizes the geometry of a micro-solenoid to provide the maximum SNR for Magnetic resonance microscopy expeirments. It takes various input of a micro-solenoid coil geometry such as: number of turns, coil diameter, wire diameter, conductivity of the wire,  and length of the lead wires.  The other inputs are the sample length, sample diameter, temperature of the coil, temperature of the sample, and operating frequency of the NMR spectrometer. The **NturnOptimization.m** utilizes the class file **solenoidOptimizer.m** to simulate the SNR of micro-solenoids with different number of turns and wire diameter while keeping the rest of the geometry parameter constant. The **RFfield_map.m** generates the B1 sensitivity map of the solenoid coil. The **coilOptimizer.m** utilizes the all of the above matlab files to calculate required coil length, generate optimal number of turns versus wire diameter plot by iterative numerical simulation, and lastly generating the B1 map. 
The final result is similar to Figure 12 in the Minard and Wind (2001) paper. The MATLAB class file **solenoidOptimizer.m** can be edited if a different circuit design apart from Minard and Wind's equivalent circuit is being used. Lastly, the current state of **solenoidOptimizer.m** performs the net loss calculation assuming a conductive sample. For simulating the SNR of the micro-solenoids with non-conductive sample, the **Rnmr** parameter should be changed to include just the coil resistance, lead resistance, and the capacitive losses. 

## Authors
**Bibek Dhakal**  
*PhD Candidate, Department of Physics, Vanderbilt University*  
*Trainee, Vanderbilt University Institute of Imaging Science (VUIIS)*  
Email: bibek.dhakal@vanderbilt.edu

**Benjamin Hardy**  
*Remcom Inc.*  
*VUIIS alumnus*

## Installation:
1. Download and unzip the **microSolenoidGeometryOptimizer** repository
2. Store it in a folder
3. Add the path of the folder to matlab path
4. open the **coilOptimizer.m** file in the Matlab
5. run and enjoy!

## References: 
[1] Minard, K. R. & Wind, R. A. Solenoidal microcoil design. Part II: Optimizing winding parameters for maximum signal-to-noise performance. Concepts Magn. Reson. 13, 190–210 (2001).
