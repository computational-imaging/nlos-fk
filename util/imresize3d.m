function out = imresize3d(vol, Nx, Ny, Nz)
    x = linspace(0, 1, size(vol, 1));
    y = linspace(0, 1, size(vol, 2));
    z = linspace(0, 1, size(vol, 3));
    [X, Y, Z] = ndgrid(x, y, z);
    F = griddedInterpolant(X, Y, Z, vol, 'nearest', 'none');
    
    xq = linspace(0, 1, Nx);
    yq = linspace(0, 1, Ny);
    zq = linspace(0, 1, Nz);
    [Xq, Yq, Zq] = ndgrid(xq, yq, zq);
    out = F(Xq, Yq, Zq);
end
