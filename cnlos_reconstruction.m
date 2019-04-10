function result = cnlos_reconstruction(meas, tofgrid, wall_size, alg, crop)
% Reconstruction procedures for "Confocal Non-Line-of-Sight Imaging Based on the Light Cone Transform"
% by Matthew O'Toole, David B. Lindell, and Gordon Wetzstein.
% and for "Wave-Based Non-Line-of-Sight Imaging using Fast f-k Migration"
% by David B. Lindell, Matthew O'Toole, and Gordon Wetzstein.
% 
% Usage/Parameters: 
% cnlos_reconstruction(meas, tofgrid, wall_size, alg)
%     meas: measurements provided from the data files for each scene
%     tofgrid: time-of-flight delays in picoseconds provided by the associated
%          datafiles for each scene. Pass in an empty array if data is
%          already corrected.
%     wall_size: size of scanned wall along one dimension (assumes 
%         identical horizontal, vertical scanned dimensions)
%     alg == 0 : FBP
%     alg == 1 : LCT
%     alg == 2 : f-k migration
%     crop: optional argument, bin index to crop measurements after correcting
%         for time-of-flight delays

    % Constants
    bin_resolution = 32e-12; % Native bin resolution for SPAD is 4 ps
    c              = 3e8;    % Speed of light (meters per second)
    width = wall_size / 2;
    if ~exist('crop', 'var') % bin index to crop measurements after aligning
        crop = 512;          % so that direct component is at t=0
    end
    
    % Parameters
    isdiffuse  = 1; % Toggle diffuse reflection (LCT only)
    isbackproj = 0; % Toggle backprojection vs LCT (LCT only)
    if alg == 0
       isbackproj = 1; 
    end
    snr = 1e-1; % SNR value (LCT only)

    % adjust so that t=0 is when light reaches the scan surface
    if ~isempty(tofgrid)
        for ii = 1:size(meas, 1)
            for jj = 1:size(meas,2 )
                meas(ii, jj, :) = circshift(meas(ii, jj, :), -floor(tofgrid(ii, jj) / (bin_resolution*1e12)));
            end
        end  
    end
    meas = meas(:, :, 1:crop);
    
    N = size(meas,1);        % Spatial resolution of data
    M = size(meas,3);        % Temporal resolution of data
    range = M.*c.*bin_resolution; % Maximum range for histogram
    
    % Permute data dimensions
    data = permute(meas,[3 2 1]);

    % Define volume representing voxel distance from wall
    grid_z = repmat(linspace(0,1,M)',[1 N N]);

    display('Inverting...');
    tic;

    if (alg == 2) % f-k migration
        [z,y,x] = ndgrid(-M:M-1,-N:N-1,-N:N-1);
        z = z./M; y = y./N; x = x./N;

        % Step 0: Pad data
        data = data .* grid_z.^2;
        data = sqrt(data);
        tdata = zeros(2.*M,2.*N,2.*N);
        tdata(1:end./2,1:end./2,1:end./2) = data;

        % Step 1: FFT
        tdata = fftshift(fftn(tdata));

        % Step 2: Stolt trick
        tvol = interpn(z,y,x,tdata,sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(x.^2+y.^2)+z.^2)),y,x,'linear',0);
        tvol = tvol.*(z > 0);
        tvol = tvol.*abs(z)./max(sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(x.^2+y.^2)+z.^2)),1e-6);

        % Step 3: IFFT
        tvol = ifftn(ifftshift(tvol));
        tvol = abs(tvol).^2;    
        vol = abs(tvol(1:end./2,1:end./2,1:end./2));

    else % LCT     
        % Define NLOS blur kernel 
        psf = definePsf(N,M,width./range);

        % Compute inverse filter of NLOS blur kernel
        fpsf = fftn(psf);
        if isbackproj
            invpsf = conj(fpsf);
        else
            invpsf = conj(fpsf) ./ (abs(fpsf).^2 + 1./snr);
        end
        
        % Define transform operators
        [mtx,mtxi] = resamplingOperator(M);
        mtx = full(mtx);
        mtxi = full(mtxi);

        % Step 1: Scale radiometric component
        if (isdiffuse)
            data = data.*(grid_z.^4);
        else
            data = data.*(grid_z.^2);
        end

        % Step 2: Resample time axis and pad result
        tdata = zeros(2.*M,2.*N,2.*N);
        tdata(1:end./2,1:end./2,1:end./2)  = reshape(mtx*data(:,:),[M N N]);

        % Step 3: Convolve with inverse filter and unpad result
        tvol = ifftn(fftn(tdata).*invpsf);
        tvol = tvol(1:end./2,1:end./2,1:end./2);

        % Step 4: Resample depth axis and clamp results
        vol  = reshape(mtxi*tvol(:,:),[M N N]);
        vol  = max(real(vol),0);
        
        if isbackproj
           % apply laplacian of gaussian filter
           vol = filterLaplacian(vol);
           
           % normalize, apply optional threshold           
           vol = max(vol./max(vol(:)),0);
        end
    end

    display('... done.');
    time_elapsed = toc;

    display(sprintf(['Reconstructed volume of size %d x %d x %d '...
        'in %f seconds'], size(vol,3),size(vol,2),size(vol,1),time_elapsed));

    tic_z = linspace(0,range./2,size(vol,1));
    tic_y = linspace(width,-width,size(vol,2));
    tic_x = linspace(width,-width,size(vol,3));

    % clip artifacts at boundary, rearrange for visualization
    vol(end-10:end, :, :) = 0;
    vol = permute(vol, [1, 3, 2]);
    result = permute(vol, [2, 3, 1]);
    vol = flip(vol, 2);
    vol = flip(vol, 3);

    % View result
    figure

    subplot(1,3,1);
    imagesc(tic_x,tic_y,squeeze(max(vol,[],1)));
    title('Front view');
    set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
    set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
    xlabel('x (m)');
    ylabel('y (m)');
    colormap('gray');
    axis square;

    subplot(1,3,2);
    imagesc(tic_x,tic_z,squeeze(max(vol,[],2)));
    title('Top view');
    set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
    set(gca,'YTick',linspace(min(tic_z),max(tic_z),3));
    xlabel('x (m)');
    ylabel('z (m)');
    colormap('gray');
    axis square;

    subplot(1,3,3);
    imagesc(tic_z,tic_y,squeeze(max(vol,[],3))')
    title('Side view');
    set(gca,'XTick',linspace(min(tic_z),max(tic_z),3));
    set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
    xlabel('z (m)');
    ylabel('y (m)');
    colormap('gray');
    axis square;

    drawnow;

end

function psf = definePsf(U,V,slope)
    % Local function to compute NLOS blur kernel
    x = linspace(-1,1,2.*U);
    y = linspace(-1,1,2.*U);
    z = linspace(0,2,2.*V);
    [grid_z,grid_y,grid_x] = ndgrid(z,y,x);

    % Define PSF
    psf = abs(((4.*slope).^2).*(grid_x.^2 + grid_y.^2) - grid_z);
    psf = double(psf == repmat(min(psf,[],1),[2.*V 1 1]));
    psf = psf./sum(psf(:,U,U));
    psf = psf./norm(psf(:));
    psf = circshift(psf,[0 U U]);
end

function [mtx,mtxi] = resamplingOperator(M)
 % Local function that defines resampling operators
     mtx = sparse([],[],[],M.^2,M,M.^2);
     
     x = 1:M.^2;
     mtx(sub2ind(size(mtx),x,ceil(sqrt(x)))) = 1;
     mtx  = spdiags(1./sqrt(x)',0,M.^2,M.^2)*mtx;
     mtxi = mtx';
     
     K = log(M)./log(2);
     for k = 1:round(K)
         mtx  = 0.5.*(mtx(1:2:end,:)  + mtx(2:2:end,:));
         mtxi = 0.5.*(mtxi(:,1:2:end) + mtxi(:,2:2:end));
     end
end

function out = filterLaplacian(vol)
    w = fspecial3('log',5, 1);
    out = convn(vol, w, 'same');
end