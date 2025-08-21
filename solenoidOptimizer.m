classdef solenoidOptimizer
    % solenoidOptimizerobj Class
    %
    % Description:
    %   This class models a solenoid coil and calculates its properties, including
    %   resistance, inductance, capacitance, magnetic field, and NMR Signal-to-Noise Ratio (SNR).
    %   It considers factors like coil geometry, material properties, operating frequency,
    %   and sample characteristics.
    %
    % Inputs (Constructor):
    %   n             : Number of turns in the coil.
    %   dcoil         : Outer diameter of the coil (meters).
    %   lcoil         : Length of the coil (meters).
    %   dwire         : Diameter of the wire used for the coil (meters).
    %   leadwirelength: combined length of the two lead wires (meters).
    %   f             : Operating frequency (Hz).
    %   rho           : Resistivity of the wire material (Ohm-meters).
    %   alpha         : Diameter of the sample (meters).
    %   beta          : Length of the sample (meters).
    %   sigma         : Conductivity of the sample (Siemens/meter).
    %   tc            : Temperature of the coil (Kelvin).
    %   ts            : Temperature of the sample (Kelvin).
    %
    % Outputs (Properties):
    %   n             : Number of turns in the coil.
    %   dcoil         : Outer diameter of the coil (meters).
    %   dwire         : Diameter of the wire used for the coil (meters).
    %   lcoil         : Length of the coil (meters).
    %   leadwirelength: Length of the lead wires (meters).
    %   s             : Spacing between turns (meters).
    %   u0            : Permeability of free space (H/m).
    %   rho           : Resistivity of the wire material (Ohm-meters).
    %   ur            : Relative permeability (unitless).
    %   f             : Operating frequency (Hz).
    %   tc            : Temperature of the coil (Kelvin).
    %   sigma         : Conductivity of the sample (Siemens/meter).
    %   alpha         : Diameter of the sample (meters).
    %   beta          : Length of the sample (meters).
    %   ts            : Temperature of the sample (Kelvin).
    %   Bxy           : Magnetic field strength (Tesla).
    %   en            : Enhancement factor (unitless).
    %   Rcoil         : Resistance of the coil (Ohms).
    %   Rleads        : Resistance of the lead wires (Ohms).
    %   Lcoil         : Inductance of the coil (Henry).
    %   Rcap          : Resistance due to capacitance (Ohms).
    %   RsampleDielectric: Resistance due to dielectric losses in the sample (Ohms).
    %   RsampleMagnetic: Resistance due to magnetic losses in the sample (Ohms).
    %   Rnmr          : Total resistance (Ohms).
    %   RnmrTemp      : Temperature-adjusted total resistance (Ohms).
    %   SNR           : Signal-to-Noise Ratio (unitless).
    %   SNRtemp       : Temperature-adjusted Signal-to-Noise Ratio (unitless).
    %
    % Methods:
    %   solenoidOptimizerobj: Constructor to initialize the class object.
    %   calculateEnhancement: Calculates the enhancement factor.
    %   calculateRcoil: Calculates the coil resistance.
    %   calculateRleads: Calculates the lead wire resistance.
    %   calculateInductance: Calculates the coil inductance.
    %   calculateCapacitance: Calculates the resistance due to capacitance.
    %   calcualteRsampleMagnetic: Calculates the magnetic sample resistance.
    %   calculateRsampleDielectric: Calculates the dielectric sample resistance.
    %   calculateSNR: Calculates the Signal-to-Noise Ratio (SNR).

    % References: .	Minard, K. R. & Wind, R. A. Solenoidal microcoil design. Part II: 
    % Optimizing winding parameters for maximum signal-to-noise performance. Concepts Magn. Reson. 13, 190–210 (2001).
%
% bibek.dhakal@vanderbilt.edu
% Vanderbilt Univeristy Institute of Imaging Science (VUIIS)
% Date last updated: March 26, 2025
% Date created: 22 January, 2024
%
% Acknowledgement: Thank you to Ben M. Hardy (Remcom Inc. and VUIIS Alumni) 
% who wrote the initial version of this code.

    properties
        % Coil parameters
        n, dcoil, dwire, lcoil, leadwirelength, s, % Number of turns, coil diameter, wire diameter, coil length, lead wire length, spacing
        
        % Physical constants:
        % Permeability of free space, resistivity, relative permeability, frequency
        u0, rho, ur, f, 
        
        % Temperature and sample parameters:
        % Coil temperature, conductivity, sample diameter, sample length, sample temperature
        tc, sigma, alpha, beta, ts, 
        
        % Calculated values:
        % Magnetic field, enhancement factor, coil resistance, lead resistance, coil inductance
        Bxy, en, Rcoil, Rleads, Lcoil, 
        
        % Capacitor resistance, dielectric sample resistance, magnetic sample resistance, 
        % total resistance, temp-adjusted total resistance, SNR, temp-adjusted SNR
        Rcap, RsampleDielectric, RsampleMagnetic, Rnmr, RnmrTemp, SNR, SNRtemp 
    end
    
    methods
        % Permeability of free space, resistivity, relative permeability, frequency
        function self = solenoidOptimizer(n, dcoil, lcoil, dwire,leadwirelength, f, rho, alpha, beta, sigma, tc, ts)
           % Assign input parameters to object properties
            self.n = n;
            self.dcoil = dcoil;
            self.dwire = dwire;
            self.lcoil = lcoil;
            self.leadwirelength = leadwirelength;
            
            % Calculate derived parameters
            self.s = self.lcoil / self.n; % Calculate spacing between turns
            self.u0 = 4 * pi * 1e-7; % Permeability of free space
            self.rho = rho; % Resistivity of wire
            self.ur = 1; % Relative permeability (assuming air core)
            self.f = f; % Frequency
            self.tc = tc; % Coil temperature
            self.sigma = sigma; % Conductivity of sample
            self.alpha = alpha; % Sample diameter
            self.beta = beta; % Sample length
            self.ts = ts; % Sample temperature
            
            % Calculate magnetic field
            self.Bxy = (self.n * self.u0) / (self.dcoil * sqrt(1 + (self.lcoil / self.dcoil)^2));
            
            % Calculate other dependent properties
            self.en = self.calculateEnhancement();
            self.Lcoil = self.calculateInductance();
            self.Rcoil = self.calculateRcoil();
            self.Rleads = self.calculateRleads();
            self.RsampleMagnetic = self.calcualteRsampleMagnetic();
            self.Rcap = self.calculateRcapacitance();
            self.RsampleDielectric = self.calculateRsampleDielectric();
            
            % Calculate SNR
            [self.Rnmr, self.RnmrTemp, self.SNR, self.SNRtemp] = self.calculateSNR(); % Calculate and assign SNR
        end
        
        %Enhancement factor calculation
        function en = calculateEnhancement(self)
            dos = self.dwire / self.s;
            dos = min(max(round(dos, 4, 'significant'), 0.1), 1);
            coil_ratio = round(self.lcoil / self.dcoil, 4, 'significant');
            
            crs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 1000000];
            dss = flip(0.1:0.1:1);
            [X, Y] = meshgrid(crs, dss);
            [Xq, Yq] = meshgrid(0:0.01:100, 1:-0.01:0.1);
            Xq = round(Xq,2);
            Yq = round(Yq, 2);
            
            ef_table = [5.31, 5.45, 5.65, 5.80, 5.80, 5.55, 4.10, 3.54, 3.31, 3.20, 3.23, 3.41; ...
            3.73, 3.84, 3.99, 4.11, 4.17, 4.10, 3.36, 3.05, 2.92, 2.90, 2.93, 3.11;...
                2.74, 2.83, 2.97, 3.10, 3.20, 3.17, 2.74, 2.60, 2.60, 2.62, 2.65, 2.81;...
                    2.12, 2.20, 2.28, 2.38, 2.44, 2.47, 2.32, 2.27, 2.29, 2.34, 2.37, 2.51;...
                        1.74, 1.77, 1.83, 1.89, 1.92, 1.94, 1.98, 2.01, 2.03, 2.08, 2.10, 2.22;...
                            1.44, 1.48, 1.54, 1.60, 1.64, 1.67, 1.74, 1.78, 1.80, 1.81, 1.83, 1.93;...
                                1.26, 1.29, 1.33, 1.38, 1.42, 1.45, 1.50, 1.54, 1.56, 1.57, 1.58, 1.65;...
                                    1.16, 1.19, 1.21, 1.22, 1.23, 1.24, 1.28, 1.32, 1.34, 1.34, 1.35, 1.40;...
                                        1.07, 1.08, 1.08, 1.10, 1.10, 1.10, 1.13, 1.15, 1.16, 1.16, 1.17, 1.19;...
                                            1.02, 1.02, 1.03, 1.03, 1.03, 1.03, 1.04, 1.04, 1.04, 1.04, 1.04, 1.05]; % Enhancement factor table (same as in original code)
            ef = griddata(X, Y, ef_table, Xq, Yq, 'linear');
    
            coil_ratio = round(coil_ratio,2);
            dos = round(dos,2);
    
    
            [row1,col1] = find(Xq==coil_ratio);
            index_coil_ratio = col1(1);
            
            [row2,col2] = find(Yq==dos);
            index_dos = row2(1);
    
            ef_needed = ef(index_dos,index_coil_ratio); %tricky, huh! :)
            en = 1 + (ef_needed - 1) * (1 - (1 / self.n));
        end
        
        %Inductanance Calculation
        function Lcoil = calculateInductance(self)
            %Lead inductance
            Llead = 460 * self.leadwirelength * log10(4 * (self.leadwirelength / self.dwire) - 0.75);
            %coil inductance
            J = 2.33 * log10(self.dwire / self.s) + 0.515;
            K_v = [0.001, 0.008, 0.16, 0.19, 0.22, 0.23, 0.24, 0.25, 0.26, 0.265, 0.27, 0.275, 0.28];

            if (self.n > 13)%from figure 3 of the paper
                K = .3;
            else
                K = K_v(self.n);
            end
            
            L = (9850 * self.dcoil * self.n^2 / (4.5 + 10 * (self.lcoil / self.dcoil))) - 628 * self.dcoil * self.n * (J + K);
            Lcoil = (L + Llead) * 1e-9;
        end

        %Coil Resistance Calculation
        function Rcoil = calculateRcoil(self)
            l = self.lcoil + self.n * self.dcoil * pi;
            rs = ((l / self.dwire) * sqrt(self.ur * self.u0 * self.rho * self.f / pi));
            Rcoil = rs * self.en;
        end
        
        %Coil lead resitance calculation
        function Rleads = calculateRleads(self)
            Rleads = (self.leadwirelength / self.dwire) * sqrt(self.ur * self.u0 * self.rho * self.f / pi);
        end
        


        %calculate Rcap
        function Rcap = calculateRcapacitance(self)
            f0 = self.f/1e06; %in Mhz
            Q = 5.05*(10^9)*((f0)^-2.35); %eqn 13
            L = self.Lcoil;
            %need to estimate capacitance!
            Ctune = 1/(L*(2*pi*self.f)^2);
            %Ctune = 3.5*1e-12;%Farad from the paper Capacitive loss section
            Rcap =  1/(self.f*2*pi*Q*Ctune);
        end

        %calculate RsampleMagnetic
        function RsampleMagnetic = calcualteRsampleMagnetic(self)
            w0 = 2*pi*self.f;
            RsampleMagnetic = (pi*(w0^2)*(self.u0^2)*(self.n^2)*(self.alpha^4)*self.beta*self.sigma)/(128*((self.dcoil^2) + (self.beta^2)));
        end

        %calculate RsampleDielectric
        function RsampleDielectric = calculateRsampleDielectric(self)
            w0 = 2*pi*self.f; %Scanner's frequency
            enot = 8.85e-12; % F/m; permittivity of free space
            e0 = 78.32; %for water
            einf = 5.30; %for water at infinity
            tau = 8.27e-12; % pico second; from paper reference(40)
            epsPrime = einf + (e0 - einf)/(1+(w0*tau)^2);
            epsDoublePrime = ((e0-einf)*w0*tau)/(1+(w0*tau)^2) + self.sigma/(w0*enot);
            
            %!!!careful!!!
            fd = .948; % in van heterens work he reports .25 for larger coils 
            % !!!!!!!!careful!!!!
   
            
            H = .1126*(self.lcoil/self.dcoil) + .08 + (.27/(sqrt(self.lcoil/self.dcoil)));
            Cstray = 100*H*self.dcoil* 1e-12; % change to Farads
            Cprime = .5*Cstray;
            
            C1 = Cprime*(1/(1-fd)); % equation 23
            C2 = epsPrime*Cprime*(1/fd); % equation 23
            
            Rd = epsPrime/(w0*C2*epsDoublePrime);
            Yreal = (Rd*(w0*C1)^2)/(1+(((C1+C2)^2)*((w0*Rd)^2))); % equation 20
            Re = Yreal * self.Lcoil^2 * w0^2; % equation 19
            RsampleDielectric = Re;
        end

        %Total Reistance
        function [Rnmr, RnmrTemp, SNR, SNRtemp] = calculateSNR(self)

            %calculate net loss
            %self.Rnmr = (self.Rcoil + self.Rleads+self.Rcap);%for non-conductive sample
            Rnmr = (self.Rcoil+ self.Rleads+ self.Rcap + self.RsampleDielectric+ self.RsampleMagnetic); %conductive sample
    
            %self.RnmrTemp = (self.Rcoil + self.Rcap + self.Rleads + self.RsampleMagnetic);
            RnmrTemp = ((self.Rcoil + self.Rcap + self.Rleads)*self.tc + (self.RsampleMagnetic + self.RsampleDielectric)*self.ts);
        
               
            %calculate the final SNR
            SNR = (self.Bxy*(2*pi*self.f)^2)/sqrt(Rnmr);
            
            %calculate the SNR addressing temperature
            SNRtemp = (self.Bxy*(2*pi*self.f)^2)/sqrt(RnmrTemp);
            
        end
        
    end
end



