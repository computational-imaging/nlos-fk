%% Run this demo script to step through the included reconstruction procedures

% First, run FBP, LCT, and f-k migration reconstructions for one of the
% captured datasets

% Optionally replace the below filenames with files from other scenes:
% bike, discoball, dragon, outdoor, resolution, statue, teaser
load('statue/tof.mat');
load('statue/meas_10min.mat');

% resize to low resolution to reduce memory requirements
measlr = imresize3(meas, [64, 64, 2048]); % y, x, t
tofgridlr = imresize(tofgrid, [64, 64]); 
wall_size = 2; % scanned area is 2 m x 2 m

% run FBP
fprintf('\nRunning FBP\n');
algorithm = 0;
fbp = cnlos_reconstruction(measlr, tofgridlr, wall_size, algorithm);

% run LCT
fprintf('\nRunning LCT\n');
algorithm = 1;
lct = cnlos_reconstruction(measlr, tofgridlr, wall_size, algorithm);

% run f-k migration
fprintf('\nRunning f-k migration\n');
algorithm = 2;
fk = cnlos_reconstruction(measlr, tofgridlr, wall_size, algorithm);

%% Reconstruct a frame from the interactive results
% perform reconstruction and visualization for 32 x 32 resolution scene
% captured at 4 Hz
fprintf('\nReconstructing interactive results\n');

% this processes the frame, but it could also be directly loaded from the
% preprocessed version in 'interactive/fk_32.mat'
load('interactive/meas_32.mat');
load('interactive/tof_32.mat');
v = VideoReader('interactive/video_32.mov');

video_idx = 1400;
frame_idx = 73;
frame = squeeze(meas(frame_idx, :, :, :));
crop = 1024; % crop measurements to 1024 bins
cnlos_reconstruction(frame, tofgrid, wall_size, algorithm, crop);
figure; imshow(read(v, video_idx),'InitialMagnification', 20); 
title('Photo of interactive capture');
clear v;

%% Run non-planar reconstruction
% reconstruction of a scene from a scanned non-planar surface
fprintf('\nRunning non-planar reconstructions\n');
nonplanar_reconstruction;

%% Run phase retrieval
% run the phase retrieval algorithm on one of the captured datasets
fprintf('\nRunning phase-retrieval reconstruction\n');
hio_reconstruction;