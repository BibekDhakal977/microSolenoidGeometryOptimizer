function [best_n, best_dwire,peakValsAll,peakDwireAll] = NturnOptimization(n, dcoil_inner, lcoil, leadwirelength, dsample, lsample, Tsample, Tcoil, sample_conductivity, wire_resistivity, operating_frequency)
   

    
    SNR = zeros(length(n), 20); % Pre-allocate SNR storage
    % Loop over coil turns
    for i = 1:length(n)
        % Define a range of wire diameters
        dwire = linspace(0.1 * lcoil / n(i), lcoil / n(i), 20); 
        
        % Loop over wire diameters
        for k = 1:length(dwire)
            dcoil = dcoil_inner + dwire(k); % Outer coil diameter
            
            % Create an instance of solenoidOptimizer class
            coil = solenoidOptimizer(n(i), dcoil, lcoil, dwire(k),leadwirelength, operating_frequency, wire_resistivity, dsample, lsample, sample_conductivity, Tcoil, Tsample);
            
            % Store SNR value from the class property
            SNR(i, k) = coil.SNR;
            
            % Display current values
            %fprintf('n = %d, dwire = %.2f µm\n', n(i), dwire(k) * 1e6);
        end
    end
    
    %%
    %Normalizing SNR
    SNR = SNR/min(SNR(:));
    
    %% Find the global maximum SNR across all n and dwire
    [maxVal, maxIdx] = max(SNR(:));          % max SNR and its index in the 2D array
    [minVal, minIdx] = min(SNR(:));
    [iMax, kMax] = ind2sub(size(SNR), maxIdx); % convert linear index to row,col
    
    % Reconstruct the dwire range for the row iMax
    dwireRange = linspace(0.1 * lcoil / n(iMax), lcoil / n(iMax), 20);
    
    % Get the best n
    best_n     = n(iMax);
    best_dwire = dwireRange(kMax);
    
    fprintf('\nGlobal max SNR = %.3f at n = %d and dwire = %.2f um\n', ...
        maxVal, best_n, best_dwire*1e6);
    %% **Plot Results**
    figure;
    legend_entries = cell(length(n), 1); % Prepare legend entries
    
    for i = 1:length(n)
        dwire = linspace(0.1 * lcoil / n(i), lcoil / n(i), 20); % Wire diameters
        
        %plot
        plot(dwire / 1e-6, SNR(i, :), '-', 'LineWidth', 2);
        grid on;
        hold on;
        legend_entries{i} = ['n = ' num2str(n(i))]; % Store legend text
    end
    ylim([0.9*minVal 1.1*maxVal]);
    plot(best_dwire / 1e-6,maxVal,'^r','LineWidth',2);
    legend(legend_entries,'FontSize',10,'FontWeight','bold');
    
    % Turn off auto-updating the legend
    legend('AutoUpdate','off'); 
    
    % Draw your projection lines, which will be excluded
    xPt = best_dwire/1e-6;
    yPt = maxVal;
    plot([xPt xPt], [0.9*min(SNR(:)) yPt], '--k');
    plot([0 xPt], [yPt yPt], '--k');
    
    % Make the axes lines thicker and turn on minor ticks and minor grids
    ax = gca; 
    ax.LineWidth = 1.5; 
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    
    
    xlabel('d_{wire} (\mum)','FontSize',12,'FontWeight','bold');
    ylabel('SNR_{Normalized}','FontSize',12,'FontWeight','bold');
    title('SNR vs Wire Diameter for Different Coil Turns','FontSize',12,'FontWeight','bold');


    
    %% Find local maxima in each row

for iRow = 1:length(n)
    % Extract the 1D SNR data for row iRow
    SNR_row = SNR(iRow, :);

    % Find local maxima in that row
    [peakVals, peakLocs] = findpeaks(SNR_row);

    % Reconstruct the dwire range for this n
    dwireRange = linspace(0.1*lcoil/n(iRow), lcoil/n(iRow), 20);

    if ~isempty(peakLocs)
        % The wire diameters at local maxima
        localDwirePeaks = dwireRange(peakLocs);

        % Store them in cell arrays
        peakValsAll(iRow) = peakVals;
        peakDwireAll(iRow) = localDwirePeaks;

        fprintf('\nLocal maxima for n=%d:\n', n(iRow));
        for p = 1:length(peakLocs)
            fprintf('  SNR=%.3f at dwire=%.3f µm (index %d)\n',...
                peakVals(p), localDwirePeaks(p)*1e6, peakLocs(p));
        end
    else
        % No local maxima
        peakValsAll(iRow) = [];       % empty
        peakDwireAll(iRow) = [];
        fprintf('\nNo local maxima found for n=%d\n', n(iRow));
    end
end


end