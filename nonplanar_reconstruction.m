function nonplanar_reconstruction()
% Non-planar reconstruction procedure described in "Wave-Based 
% Non-Line-of-Sight Imaging using Fast f-k Migration" by David B. Lindell, 
% Matthew O'Toole, and Gordon Wetzstein. 
%
% This function loads the measurements and calibration files and displays
% the reconstructed scene with and without the non-planar correction. A 3D
% visualization of the scanned surface is also displayed. The full
% wave-based non-planar correction can be applied (default), or a faster,
% approximate histogram-shifting-based approach can be used (adjust
% parameter settings accordingly). To speed up the reconstruction, the
% volumes are processed at low resolution by default.
%
% Reconstructed volumes at full resolution are also provided in the 
% 'nonplanar' directory for the filtered backprojection approach as well 
% for f-k migration applied to the uncorrected data, with the
% naive histogram-shifting correction, or with the full wave-based 
% non-planar correction.

    % load data files
    addpath('util/');
    load('nonplanar/tof.mat');
    load('nonplanar/meas_150sec.mat');
    load('nonplanar/calibration.mat');

    % Parameters
    use_naive = 0; % toggle using wave correction or naive histogram shifting
    use_low_resolution = 1; % process at lower resolution to reduce computation time

    % Constants
    bin_resolution = 32e-12; 
    c              = 3e8; 

    % realign measurements so that time zero corresponds to arrival of direct
    % component to the surface
    rect_data = meas;
    for ii = 1:size(rect_data, 1)
       for jj = 1:size(rect_data,2 )
           rect_data(ii, jj, :) = circshift(rect_data(ii, jj, :), -floor(tofgrid(ii, jj) / (bin_resolution*1e12)));
       end
    end  
    rect_data = sqrt(rect_data(:, :, 1:512));

    if use_low_resolution
        rect_data = imresize3(rect_data, [32, 32, 512]);
        z_offset = imresize(z_offset, [32, 32]);
    end

    % set dimension sizes and space/fourier domain indices
    Nx = size(rect_data, 1);
    Ny = size(rect_data, 2);
    Nz = size(rect_data, 3);
    Nx_iter = Nx;
    Ny_iter = Ny;
    Nz_iter = Nz;

    spectra = zeros(size(rect_data));
    meas_rect = zeros(size(rect_data));
    [N, ~, M] = size(rect_data);
    range = M.*c.*bin_resolution;

    x = linspace(-width/2, width/2, Nx);
    y = linspace(-width/2, width/2, Ny);
    z = linspace(0, range/2, Nz);
    [X, Y] = ndgrid(x, y);

    kx = fftshift(fftfreq(Nx, x(2) - x(1)));
    ky = fftshift(fftfreq(Ny, y(2) - y(1)));
    kz = fftshift(fftfreq(Nz, z(2) - z(1)));
    [KX, KY, KZ] = ndgrid(kx, ky, kz);
    vec = @(x) x(:);

    % run the non-planar correction
    for ii = 1:Nx_iter
        progress(ii, Nx_iter);
        for jj = 1:Ny_iter
            if use_naive
                meas_rect(ii, jj, :) = circshift(rect_data(ii, jj, :), 2*floor(z_offset(ii, jj) / (bin_resolution*c)));
            else
                tmp_vol = zeros(size(rect_data));
                tmp_vol(ii, jj, :) = rect_data(ii, jj, :);
                tmp_spectra = fftshift(fftn(tmp_vol));
                phase_mask = exp(2*pi*1j * (sqrt(KZ.^2 - KX.^2 - KY.^2)) .* z_offset(ii, jj));
                spectra = spectra + tmp_spectra.*phase_mask;
            end
        end

        if ~use_naive
            meas_rect = real( rifft3(ifftshift(spectra)));
        end
    end
    
    % run reconstruction
    algorithm = 2;
    wall_size = 1;

    fprintf('\nReconstructing corrected volume\n');
    cnlos_reconstruction(meas_rect.^2, [], wall_size, algorithm);
    fprintf('\nReconstructing volume without correction\n');
    cnlos_reconstruction(rect_data.^2, [], wall_size, algorithm);

    fprintf('\nVisualizing non-planar scan surface\n');
    load('nonplanar/scan_points.mat');

    figure;
    scan_grid(3, :) = -(scan_grid(3, :) - max(scan_grid(3, :)));
    scan_grid_pc = pointCloud(scan_grid');

    pcshow(scan_grid_pc, 'VerticalAxis', 'Y', 'MarkerSize', 14);
    set(gca, 'zlim', [0.0, 0.11]);
    set(gca, 'ztick', [0.0 0.1], 'zticklabel', {'0', '0.1'});
    xlabel('x'); ylabel('y'); zlabel('z');
    camup([0 1 0]);
    campos([-1, 1, 1])
    set(gca, 'clim', [0, 0.1]);
    title('Scan Surface');
end

function result = fftfreq(n, d)
    val = 1.0 / (n * d);
    N = floor((n-1) / 2) + 1;
    p1 = 0:N-1;
    p2 = -floor(n/2):-1;
    result = [p1 p2] .* val;
end

function spectra = rifft3(spectra)
    % do ifft in x
    spectra = ifft(spectra, [], 1);
    % ifft in y
    spectra = ifft(spectra, [], 2);

    % real ifft in z (enforce conj symmetry, then do fft)
    N = size(spectra, 3);
    spectra(:, :, 2:N/2) = conj(flip(spectra(:, :, N/2+2:end), 3));
    spectra = ifft(spectra, [], 3);
end