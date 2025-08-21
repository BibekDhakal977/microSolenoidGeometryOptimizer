function RFfield_map(nTurns, d_coil, l_coil, lsample,dsample)
    % ----------------------------------------------------
    %  Coil geometry similar to your code:
    %   - coil axis = y
    %   - we compute the B_xy = sqrt(Bx^2 + By^2) at x=0
    %   - then show fractional deviation from coil center
    % ----------------------------------------------------
    
    %% Coil Parameters
    mu0    = 4*pi*1e-7;   % Permeability of free space (H/m)
    I      = 1.0;         % Current in Amps

    % Helix discretization
    Nseg   = 3000;                         
    phiVals = linspace(0, 2*pi*nTurns, Nseg);
    dphi   = phiVals(2) - phiVals(1);
    r_coil = d_coil/2;

    % Helical coordinates (coil axis = y)
    Xc = r_coil .* cos(phiVals);
    Yc = (l_coil/(2*pi*nTurns))*phiVals - (l_coil/2);
    Zc = r_coil .* sin(phiVals);

    % Differentiate to get dL for each coil segment
    dXc = gradient(Xc, dphi);
    dYc = gradient(Yc, dphi);
    dZc = gradient(Zc, dphi);

    %% Set Up Grid for YZ plane (x=0)

    Ny = 101;  % number of points in y
    Nz = 101;  % number of points in z

    % We'll define our "y range" around +/- lsample
    % and z range around +/- (d_coil/3), or as you prefer
    yvals = linspace(-lsample, lsample, Ny);
    %zvals = linspace(-d_coil/3, d_coil/3, Nz);
    zvals = linspace(-(0.8*dsample)/2, (0.8*dsample)/2, Nz);

    % Preallocate the B map in that plane
    Bmap_YZ = zeros(Ny, Nz);

    %% Compute Bmap in the YZ plane (x=0)
    totalIterations = Ny*Nz;
    w = waitbar(0,'Computing field map in YZ plane...');
    iteration = 0;

    for iy = 1:Ny
        for iz = 1:Nz
            iteration = iteration + 1;
            if mod(iteration,1000)==0
                waitbar(iteration/totalIterations, w);
            end

            xPt = 0.0;              % fixed at x=0
            yPt = yvals(iy);
            zPt = zvals(iz);

            % Biot-Savart sum
            Bx_sum = 0; 
            By_sum = 0; 
            Bz_sum = 0;
            for k=1:Nseg
                rx = xPt - Xc(k);
                ry = yPt - Yc(k);
                rz = zPt - Zc(k);
                r3 = (rx^2 + ry^2 + rz^2)^(3/2);
                if r3 < 1e-20, continue; end

                dL = [dXc(k); dYc(k); dZc(k)];
                r_vec = [rx; ry; rz];
                cross_prod = cross(dL, r_vec);
                Bx_sum = Bx_sum + cross_prod(1)/r3;
                By_sum = By_sum + cross_prod(2)/r3;
                Bz_sum = Bz_sum + cross_prod(3)/r3;
            end
            constFactor = mu0/(4*pi) * I * dphi;
            Bx = constFactor * Bx_sum;
            By = constFactor * By_sum;
            Bz = constFactor * Bz_sum;

            % B_xy = sqrt(Bx^2 + By^2)
            Bmap_YZ(iy, iz) = sqrt(Bx^2 + By^2);
        end
    end
    close(w);

    %% Fractional Deviation from the coil center
    % The coil center is (y=0, z=0), so we check Bmap at (0,0).
    B0 = interp2(zvals, yvals, Bmap_YZ, 0, 0);
    fracDev_YZ = (Bmap_YZ - B0) / B0;
    fracDev_YZ = fracDev_YZ*100;

    %% Plot the fractional deviation in the YZ plane
    figure;
    [Ygrid, Zgrid] = meshgrid( yvals*1e3, zvals*1e3 );  % in mm
    contourf(Ygrid, Zgrid, fracDev_YZ', 100, 'LineColor','none');
    colormap('jet'); 
    cb = colorbar;

    % Colorbar label text
    cb.Label.String = '$\Delta B_{1}^{-} / B_{1}^{-} \, (\%)$';
    cb.Label.FontSize = 14;       % size for the label
    cb.Label.FontWeight = 'bold'; % weight for the label
    cb.Label.Interpreter = 'latex'; % enable LaTeX for colorbar label

    % Alternatively, set the font for colorbar tick labels, too:
    cb.FontSize = 12;     % tick label font size
    cb.FontWeight = 'bold';

    %xlabel('y (mm)');
    %ylabel('z (mm)');    
    xlabel('Axial direction (mm)','FontSize',15,'FontWeight','bold');
    ylabel('Radial direction (mm)','FontSize',15,'FontWeight','bold');
    %title('Fractional Deviation in YZ Plane (x=0)');
    title('Field homogeneity in Saggital Plane','FontSize',12,'FontWeight','bold');

    %% Right after the contourf(...) and labeling:
    
    hold on;
    
    % Convert coil endpoints to mm (since you used *1e3 in the meshgrid):
    coilStart_mm = -0.5 * l_coil * 1e3;
    coilEnd_mm   = +0.5 * l_coil * 1e3;
    
    % We'll define the vertical extents in the plot for the lines.
    % Since your Zgrid is from zvals(1)*1e3 to zvals(end)*1e3:
    zMin_mm = min(Zgrid(:));  % in mm
    zMax_mm = max(Zgrid(:));  % in mm
    
    % Draw dashed lines at coil boundaries:
    plot([coilStart_mm coilStart_mm], [zMin_mm zMax_mm], 'k--','LineWidth',2);
    plot([coilEnd_mm coilEnd_mm], [zMin_mm zMax_mm], 'k--','LineWidth',2);
    
    % Optional text labels:
    text(coilStart_mm, zMax_mm, 'Coil Start', ...
         'Color','k','FontWeight','bold', ...
         'HorizontalAlignment','center','VerticalAlignment','bottom');
    text(coilEnd_mm, zMax_mm, 'Coil End', ...
         'Color','k','FontWeight','bold', ...
         'HorizontalAlignment','center','VerticalAlignment','bottom');


    %plot sample as well
     hold on;
    
    % Convert coil endpoints to mm (since you used *1e3 in the meshgrid):
    sampleStart_mm = -0.5 * lsample * 1e3;
    sampleEnd_mm   = +0.5 * lsample * 1e3;
    
    % We'll define the vertical extents in the plot for the lines.
    % Since your Zgrid is from zvals(1)*1e3 to zvals(end)*1e3:
    zMin_mm = min(Zgrid(:));  % in mm
    zMax_mm = max(Zgrid(:));  % in mm
    
    % Draw dashed lines at coil boundaries:
    plot([sampleStart_mm sampleStart_mm], [zMin_mm zMax_mm], 'b--','LineWidth',2);
    plot([sampleEnd_mm sampleEnd_mm], [zMin_mm zMax_mm], 'b--','LineWidth',2);
    
    % Optional text labels:
    text(sampleStart_mm, zMax_mm, 'Sample Start', ...
         'Color','b','FontWeight','bold', ...
         'HorizontalAlignment','center','VerticalAlignment','bottom');
    text(sampleEnd_mm, zMax_mm, 'Sample End', ...
         'Color','b','FontWeight','bold', ...
         'HorizontalAlignment','center','VerticalAlignment','bottom');



end