function geom = fkback3(meas, width, range)
    N = size(meas,2);        % Spatial resolution of data
    M = size(meas,1);        % Temporal resolution of data
    data = meas;
    
    [z,y,x] = ndgrid(-M:M-1,-N:N-1,-N:N-1);
    z = z./M; x = x./N; y = y./N;
    
    % Step 0: Pad data
    tdata = zeros(2.*M,2.*N,2.*N);
    tdata(1:end/2,1:end/2,1:end/2) = data;

    % Step 1: FFT
    tdata = fftshift(fftn(tdata));

    % Step 2: Stolt trick
    tvol = interpn(z,y,x,tdata,sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(x.^2+y.^2)+z.^2)),y,x,'linear',0);
    tvol = tvol.*abs(z)./max(sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(x.^2+y.^2)+z.^2)),1e-8);
    tvol = tvol.*(z > 0);
    geom = ifftn(ifftshift(tvol));
    geom = geom(1:end/2,1:end/2,1:end/2);
end