function meas = fkforward3(geom, width, range)
    N = size(geom,2);        % Spatial resolution of data
    M = size(geom,1);        % Temporal resolution of data
    data = geom;

    [z,y,x] = ndgrid(-M:M-1,-N:N-1,-N:N-1);
    z = z./M; x = x./N; y = y./N;

    % Step 0: Pad data
    tdata = zeros(2.*M,2.*N,2.*N);
    tdata(1:end./2,1:end./2,1:end./2) = data;
    
    % Forward
    tvol = fftshift(fftn(tdata));
    tdata = interpn(z,y,x,tvol,sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(x.^2+y.^2)-z.^2)),y,x,'linear',0);
    tdata = tdata.*(abs((((N.*range)./(M.*4*width)).^2).*(x.^2+y.^2))<abs(z.^2));
    tdata = tdata.*(z > 0);
    tdata = ifftn(ifftshift(tdata));
    meas = tdata(1:end/2,1:end/2,1:end/2);
end
