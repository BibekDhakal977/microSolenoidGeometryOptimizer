clc;
clear;
close all;

addpath(pwd);

%---------------------------------------
% Sample / coil dimensions
%---------------------------------------
lsample = 6e-3;         % length of sample = 6 mm
dsample = 2.46e-3;         % diameter of sample = 2.73 mm
capthickness = 0.27e-3; % thickness of capillary = 0.27 mm
dcoil_inner = dsample + 2 * capthickness;  
lsamp_dcoil_r = lsample / dcoil_inner;
percentage_deviation_desired = 15;%limit of field inhomogeneity desired inside sample in percentage

%---------------------------------------
% Sweep of coil length-to-diameter ratio
%---------------------------------------
lcoil_dcoil_r = linspace(0.1, 10, 1000);

%---------------------------------------
% Compute Delta B/B ratio vs. lcoil/dcoil
%---------------------------------------
deltaBx_ratio_max = zeros(size(lcoil_dcoil_r));  % pre-allocate
for ii = 1:length(lcoil_dcoil_r)
    deltaBx_ratio_max(ii) = deltaBxyCalc(lcoil_dcoil_r(ii), lsamp_dcoil_r);
end

%---------------------------------------
% Find minimum lcoil/dcoil that gives % inhomogeneity
%---------------------------------------
[~, idx] = min(abs(deltaBx_ratio_max - percentage_deviation_desired));
minimum_lcoil_dcoil_r = lcoil_dcoil_r(idx);

%---------------------------------------
% Plot
%---------------------------------------
figure;
plot(lcoil_dcoil_r, deltaBx_ratio_max, ...
    '-k', 'LineWidth', 2, 'DisplayName','$\Delta B_{1}^{-} / B_{1}^{-}$');
hold on;

plot(minimum_lcoil_dcoil_r, deltaBx_ratio_max(idx), ...
    '^r', 'LineWidth', 2, ...
    'DisplayName', sprintf('Minimum ratio = %.2f', minimum_lcoil_dcoil_r));

% Create legend now, so it only includes the first two plot objects
legend('Location','best','FontSize',10,'FontWeight','bold','Interpreter','latex');

% Turn off auto-updating the legend
legend('AutoUpdate','off'); 

% Draw your projection lines, which will be excluded
xPt = minimum_lcoil_dcoil_r;
yPt = deltaBx_ratio_max(idx);
plot([xPt xPt], [0 yPt], '--k');
plot([0 xPt], [yPt yPt], '--k');

xlabel('l_{coil}/d_{coil}','FontSize',12,'FontWeight','bold');
ylabel('$\Delta B_{1}^{-} / B_{1}^{-} \, (\%)$','FontSize',12,'FontWeight','bold','Interpreter','latex');
title(sprintf('l_{sample}/d_{coil} = %d',round(lsamp_dcoil_r,2)),'FontSize',12,'FontWeight','bold');
ylim([0, 25]);
grid on;
hold off;

% Make the axes lines thicker and turn on minor ticks and minor grids
ax = gca; 
ax.LineWidth = 1.5; 
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

%---------------------------------------
% Compute and display required coil length
%---------------------------------------
lcoil = minimum_lcoil_dcoil_r * dcoil_inner;
disp(['The required coil length is: ' num2str(lcoil*1e3) ' mm']);


%----------------------------------------------
% Compute the best number of turns of the coil
%----------------------------------------------
Tsample = 298; % 25 degrees C in K
Tcoil = 298; % 25 degrees C in K
sample_conductivity = 1; % S/m for conductive sample
wire_resistivity = 1.72e-8; % Ohm·m for copper

% Micro-Solenoid Coil properties
leadwirelength = 16e-3;%in 4mm in meters

%Operating frequency of the B0 field
operating_frequency = 650e6;%Hz

% Micro-Solenoid Coil properties
n = 2:1:8;

% Calculate the best number of turns for the micro-solenoid
[best_n, best_dwire,peakValsAll,peakDwireAll] = NturnOptimization(n,dcoil_inner, lcoil, leadwirelength, dsample, lsample, Tsample, Tcoil, sample_conductivity, wire_resistivity, operating_frequency);
peakValsAllnorm = peakValsAll/min(peakValsAll(:));
%---------------------------------------
%% Plot the B-field map
%---------------------------------------
%Best coil in terms of SNR
dcoil = dcoil_inner + best_dwire;
RFfield_map(best_n, dcoil, lcoil, lsample,dsample)    

%% local functions:
%---------------------------------------
% Local function for calculating dBxy/Bxy(max)
%---------------------------------------
function deltaBx_ratio_max = deltaBxyCalc(lcoil_dcoil_r, lsamp_dcoil_r)
    term1 = 1 + (1 / lcoil_dcoil_r)^2;
    term2 = (lcoil_dcoil_r + lsamp_dcoil_r) / sqrt(1 + (lcoil_dcoil_r + lsamp_dcoil_r)^2);
    term3 = (lcoil_dcoil_r - lsamp_dcoil_r) / sqrt(1 + (lcoil_dcoil_r - lsamp_dcoil_r)^2);
    deltaBx_ratio_max = 1 - 0.5 * sqrt(term1) * (term2 + term3);
    deltaBx_ratio_max = deltaBx_ratio_max*100;%percentage
end


