clc;
clear;
close all;

ts = 289; % 20 degrees C in K
tc = 293;
sigma = 1; % S/m 
rho = 1.72e-8; % Ohm·m for copper

% Micro-Solenoid Coil properties
n = [2; 4; 8; 12];
d_inner = 1e-3; % Inner diameter in meters
lcoil = 2.24 * d_inner; % Coil length
leadwirelength = 4e-3;%in 4mm in meters


% Sample size
alpha = 1e-3; % Sample diameter
beta = 2.26e-3; % Sample length

% Frequency
f = 750e6; % Hz

SNR = zeros(length(n), 20); % Pre-allocate SNR storage

% Loop over coil turns
for i = 1:length(n)
    % Define a range of wire diameters
    dwire = linspace(0.1 * lcoil / n(i), lcoil / n(i), 20); 
    
    % Loop over wire diameters
    for k = 1:length(dwire)
        dcoil = d_inner + dwire(k); % Outer coil diameter
        
        % Create an instance of SolenoidCoil class
        coil = solenoidcoil(n(i), dcoil, lcoil, dwire(k),leadwirelength, f, rho, alpha, beta, sigma, tc, ts);
        
        % Store SNR value from the class property
        SNR(i, k) = coil.SNR;
        
        % Display current values
        fprintf('n = %d, dwire = %.2f µm\n', n(i), dwire(k) * 1e6);
    end
end

%% **Plot Results**
figure;
legend_entries = cell(length(n), 1); % Prepare legend entries

for i = 1:length(n)
    dwire = linspace(0.1 * lcoil / n(i), lcoil / n(i), 20); % Wire diameters
    
    %plot
    plot(dwire / 1e-6, SNR(i, :), '-.*', 'LineWidth', 2);
    grid on;
    hold on;
    legend_entries{i} = ['n = ' num2str(n(i))]; % Store legend text
end

legend(legend_entries);
xlabel('d_{wire} (\mum)');
ylabel('SNR');
title('SNR vs Wire Diameter for Different Coil Turns');
