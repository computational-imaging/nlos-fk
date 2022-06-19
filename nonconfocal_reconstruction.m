function nonconfocal_reconstruction(varargin)
% Non-confocal reconstruction procedure described in "Wave-Based 
% Non-Line-of-Sight Imaging using Fast f-k Migration" by David B. Lindell, 
% Matthew O'Toole, and Gordon Wetzstein. 
% 
% The simulated data for this procedure can be downloaded from the Zaragoza
% NLOS dataset at https://graphics.unizar.es/nlos_dataset.html
% The direct download link to the dataset is 
% https://drive.google.com/uc?export=download&id=1HZrlcSnme0qDLWanOlHgK2OgwLBeU6JX
%
% We also provide example code for the nonconfocal dataset of Klein et al.
% This can be downloaded at https://nlos.cs.uni-bonn.de/challenges/Geometry or
% we provide a direct link to an example dataset at
% https://drive.google.com/uc?export=download&id=1DISUkCfiA_NSBRrF7w8cBXduyH_m_EQa
%
% This function loads simulated nonconfocal measurements, calculates the
% normal moveout correction, which results in emulated confocal
% measurements, and then processes the data using f-k migration. 

    fprintf('Running non-confocal processing...\n');

    if length(varargin) == 0
        method = 'zaragoza';
    else
        method = varargin{1};
    end

    % Zaragoza dataset (exhaustively scanned 2D laser + 2D camera measurements)
    if strcmp(method, 'zaragoza')
       
        fname = './Z_l[0.00,-0.50,0.00]_r[1.57,0.00,3.14]_v[0.81,0.01,0.81]_s[16]_l[16]_gs[1.00].hdf5';
        if ~exist(fname, 'file')
           error(['Download the simualted non-confocal dataset from the Zaragoza NLOS webpage' ...
                  ' at https://graphics.unizar.es/nlos_dataset.html, or use the direct download' ...
                  ' link https://drive.google.com/uc?export=download&id=1HZrlcSnme0qDLWanOlHgK2OgwLBeU6JX']);
        end

        fprintf('\tProcessing Zaragoza dataset\n')
        fprintf('\tLoading data\n');
        laserGridPositions = h5read(fname, '/laserGridPositions');
        laserGridPositions = reshape(laserGridPositions, [], 3);

        laserPosition = h5read(fname, '/laserPosition');

        cameraPosition = h5read(fname, '/cameraPosition');

        % simulated data is laser_x, laser_y, cam_x, cam_y, bounces, time
        % collapse to laser_pos, cam_pos, time
        data = h5read(fname, '/data');
        dims = size(data);
        data = reshape(data, dims(1)*dims(2), dims(3)*dims(4), dims(5), dims(6));
        data = squeeze(sum(data, 3));

        deltaT = h5read(fname, '/deltaT');

        fprintf('\tReparamaterize + normal moveout correction\n');
        meas = reparameterize(laserGridPositions, laserPosition, cameraPosition, data, deltaT);

        fprintf('\tf-k migration\n');
        alg = 2; % f-k migration
        crop = 1024; % crop measurements to 1024 bins
        wall_size = 1; % simulated scanned area is 1 x 1 m
        tofgrid = [];
        bin_resolution = double(deltaT / 3e8);
        cnlos_reconstruction(meas, tofgrid, wall_size, alg, crop, bin_resolution);

    % Klein dataset (single laser position, grid of 2D camera positions)
    elseif strcmp(method, 'klein')
        fname = 'letterk.mat';
        if ~exist(fname, 'file')
           error(['Download the simualted non-confocal dataset from the NLOS Benchmark webpage' ...
                  ' at https://nlos.cs.uni-bonn.de/challenges/Geometry, or use the direct download' ...
                  ' link https://drive.google.com/uc?export=download&id=1DISUkCfiA_NSBRrF7w8cBXduyH_m_EQa']);
        end

        fprintf('\tProcessing Klein dataset\n')
        fprintf('\tLoading data\n');
        load(fname)

        zero_append = tMin / tDelta;
        data = cat(3, zeros(size(data, 1), size(data, 2), zero_append), data);

        fprintf('\tReparameterize + normal moveout correction\n');
        R = 256;
        data = data(:, :, 1:1024);
        data = imresize3(data, [R, R, 1024]);
        [Y, X] = meshgrid(linspace(-0.256, 0.256, size(data, 1)), linspace(-0.256, 0.256, size(data, 2)));

        x = X(:);
        y = Y(:);

        mx = X(:) / 2;
        my = Y(:) / 2;

        mx = mx(:);
        my = my(:);

        midpoints = [mx my];
        midpoints = round(midpoints, 4);
        [midpoints] = unique(midpoints, 'rows');

        hx = abs(X(:)) / 2;
        hy = abs(Y(:)) / 2;

        hx = hx(:);
        hy = hy(:);

        offsets = [hx hy];
        offsets = round(offsets, 4);
        offsets = unique(offsets, 'rows');

        t = tDelta * (0:(size(data,3)-1));
        v = 1;

        meas = single(zeros(size(midpoints, 1), size(data,3)));
        num_added = zeros(size(offsets, 1), 1);
        m_added = zeros(size(midpoints, 1), 1);
        data_flat = reshape(data, size(data,1)*size(data,2), []);
        for ii = 1:length(y)
            % find closest midpoint and offset
            m = [x(ii) / 2, y(ii) / 2];
            h = abs([x(ii) / 2, y(ii) / 2]);
            [~, m_idx] = min(sqrt((m(1) - midpoints(:, 1)).^2 + (m(2) - midpoints(:, 2)).^2));  
            [~, h_idx] = min(sqrt((h(1) - offsets(:, 1)).^2 + (h(2) - offsets(:, 2)).^2));  

            % assign nmo-corrected measurement to new array
            offset_x = offsets(h_idx, 1);
            offset_y = offsets(h_idx, 2);
            val = squeeze(data_flat(ii, :));
            meas(m_idx, :) = meas(m_idx, :) + NMO(t, v, offset_x, offset_y, val);
            m_added(m_idx) = m_added(m_idx) + 1;
        end

        for ii = 1:size(m_added, 1)
                meas(ii, :) = meas(ii, :) / m_added(ii); 
        end

        meas(isnan(meas)) = 0;
        m = int32(sqrt(size(midpoints, 1)));
        meas = reshape(meas, m, m, []);

        fprintf('\tf-k migration\n');
        alg = 2; % f-k migration
        crop = 1024; % crop measurements to 1024 bins
        wall_size = 0.256; 
        tofgrid = [];
        bin_resolution = double(tDelta / 3e8);
        cnlos_reconstruction(meas, tofgrid, wall_size, alg, crop, bin_resolution);
    end
end

function meas = reparameterize(laserGridPositions, laserPosition, cameraPosition, data, deltaT)

    % correct for laser/camera positions
    x = laserGridPositions(:, 1);
    y = laserGridPositions(:, 3);
    z = laserGridPositions(:, 2);

    xc = cameraPosition(1);
    yc = cameraPosition(3);
    zc = cameraPosition(2);

    xl = laserPosition(1);
    yl = laserPosition(3);
    zl = laserPosition(2);

    dist = sqrt((x-xl).^2 + (y-yl).^2 + (z-zl).^2) + sqrt((x'-xc).^2 + (y'-yc).^2 + (z'-zc).^2);
    for ii = 1:length(y)
        for jj = 1:length(x)
            data(ii, jj, :) = circshift(data(ii, jj, :), -floor(dist(ii,jj) / (deltaT)));
        end
    end

    min_x = min(x);
    max_x = max(x);
    min_y = min(y);
    max_y = max(y);
    [X, Y] = meshgrid(linspace(min_x, max_x, sqrt(size(x, 1))), linspace(min_y, max_y, sqrt(size(y, 1))));

    mx = (X(:) + X(:)') / 2;
    my = (Y(:) + Y(:)') / 2;

    mx = mx(:);
    my = my(:);

    midpoints = [mx my];
    midpoints = round(midpoints, 3);
    [midpoints] = unique(midpoints, 'rows');

    hx = abs(X(:) - X(:)') / 2;
    hy = abs(Y(:) - Y(:)') / 2;

    hx = hx(:);
    hy = hy(:);

    offsets = [hx hy];
    offsets = round(offsets, 3);
    offsets = unique(offsets, 'rows');

    t = deltaT * (0:(size(data,3)-1));
    v = 1;

    meas = single(zeros(size(midpoints, 1), size(offsets, 1), size(data,3)));
    m_added = zeros(size(midpoints, 1), 1);
    for ii = 1:length(y)
    %     disp(ii);
        for jj = 1:length(x)
            % find closest midpoint and offset
            m = [(x(ii) + x(jj)) / 2, (y(ii) + y(jj)) / 2];
            h = abs([(x(ii) - x(jj)) / 2, (y(ii) - y(jj)) / 2]);
            [~, m_idx] = min(sqrt((m(1) - midpoints(:, 1)).^2 + (m(2) - midpoints(:, 2)).^2));  
            [~, h_idx] = min(sqrt((h(1) - offsets(:, 1)).^2 + (h(2) - offsets(:, 2)).^2));  

            % assign nmo-corrected measurement to new array
            offset_x = offsets(h_idx, 1);
            offset_y = offsets(h_idx, 2);
            val = squeeze(data(ii, jj, :));
            meas(m_idx, h_idx, :) = squeeze(meas(m_idx, h_idx, :))' + NMO(t, v, offset_x, offset_y, val);
            m_added(m_idx) = m_added(m_idx) + 1;
        end
    end

    meas = squeeze(sum(meas, 2));
    meas = reshape(meas, 31, 31, []);
    m_added = reshape(m_added, 31, 31);
    for ii = 1:size(m_added, 1)
        for jj = 1:size(m_added, 2)
            meas(ii, jj, :) = meas(ii, jj, :) / m_added(ii, jj); 
        end
    end

    meas(isnan(meas)) = 0;
    m = int32(sqrt(size(midpoints, 1)));
    meas = reshape(meas, m, m, []);

end

% perform normal moveout correction
function nmo = NMO(t, v, offset_x, offset_y, meas)
    tq = sqrt(t.^2 + 4*offset_x.^2 / v.^2 + 4*offset_y.^2 / v.^2);
    tq(isnan(tq)) = 0;
    nmo = interp1(t, meas, tq, 'linear', 0);
end
