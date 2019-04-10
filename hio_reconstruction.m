function hio_reconstruction()
% Phase retrieval algorithm described in "Wave-Based  Non-Line-of-Sight 
% Imaging using Fast f-k Migration" by David B. Lindell, Matthew O'Toole, 
% and Gordon Wetzstein. 
% 
% This function runs the hybrid input-output (HIO) phase retrieval procedure
% to reproduce the result shown in the supplementary material
%
% Reconstruction results are also provided for f-k migration without phase
% and after 50 iterations of HIO in the 'hio_init.mat', and
% 'hio_50.mat' files.
%
% The reconstruction runs at lower spatial resolution by default to 
% decrease the amount of computation time. Change this by adjusting the
% 'use_low_resolution' flag below.

    addpath('util/');
    use_low_resolution = 1; % run at lower resolution to decrease computation time

    load('phase_retrieval/meas_180min.mat');
    meas = sqrt(meas512);
    if use_low_resolution
        meas = imresize3(meas, [32, 32, 512]);
    end

    N = size(meas,1);        % Spatial resolution of data
    M = size(meas,3);        % Temporal resolution of data
    meas = permute(meas,[3 2 1]);

    % constants
    width = 1;
    bin_resolution = 32e-12; % Native bin resolution for SPAD is 4 ps
    c              = 3e8;   % Speed of light (meters per second)
    range = M.*c.*bin_resolution; % Maximum range for histogram

    % falloff correction factor
    [z,y,x] = ndgrid(0:M-1,0:N-1,0:N-1);
    meas = meas .* (z+1);

    % compute initial reconstruction
    geom = fkback3(meas, width, range);
    geom_init = geom;
    mag = meas;

    N_iters = 50;
    figure;
    for ii = 1:N_iters
        progress(ii, N_iters);   
        % recover geometry and enforce zero phase
        grad = geom - fkback3(meas, width, range);
        geom = geom - 0.005 * grad;

        % keep mag and add phase from measurements
        phase = angle(fkforward3(geom, width, range));
        meas = mag .* exp(1j .* phase);

        if mod(ii, 10) == 0
            subplot(231);
            imagesc(squeeze(max(abs(geom_init(1:end-5, :, :)).^2, [], 1)).');
            axis square;
            title('Initial Reconstruction (xy)');

            subplot(232);
            imagesc(squeeze(max(abs(geom_init(1:end-5, :, :)).^2, [], 2)).');
            axis square;
            title('(yz)');

            subplot(233);
            imagesc(squeeze(max(abs(geom_init(1:end-5, :, :)).^2, [],3)).');
            axis square;
            title('(xz)');

            subplot(234);
            imagesc(squeeze(max(abs(geom(1:end-5, :, :)).^2, [], 1)).');
            axis square;
            title(sprintf('HIO Reconstruction (xy) (%d Iters)', ii));

            subplot(235);
            imagesc(squeeze(max(abs(geom(1:end-5, :, :)).^2, [], 2)).');
            axis square;
            title('(yz)');

            subplot(236);
            imagesc(squeeze(max(abs(geom(1:end-5, :, :)).^2, [], 3)).');
            axis square;
            title('(xz)');
            drawnow;
        end
end
