function [A, ells] = get_spherical_function(dims, order, safe)
%SPH_HARM_BASIS  Real spherical harmonics basis on a disk (mapped to hemisphere).
%
%   [A, ells] = sph_harm_basis(dims, order, safe)
%
%   Inputs:
%       dims  - [Nx, Ny] image size
%       order - maximum spherical harmonic degree (l)
%       safe  - if true, excludes all axis-symmetric terms (m=0, l>0)
%
%   Outputs:
%       A     - (Nx*Ny) x Nterms basis matrix
%       ells  - vector of spherical harmonic degrees l for each column of A
%
%   Notes:
%     - Uses a disk->hemisphere mapping: r in [0,1] maps to theta in [0,pi/2].
%     - Returns real harmonics:
%          cos(mφ), sin(mφ) versions instead of complex exp(imφ).
%     - To suppress lobes, you can penalize higher l with a spectral weight,
%       e.g. R = diag((ells.*(ells+1)).^p), p=1 or 2.

    if nargin < 3, safe = false; end

    % Coordinates
    [yy,xx] = meshgrid(linspace(-1,1,dims(1)), linspace(-1,1,dims(2)));
    x = xx(:); y = yy(:);

    % Polar coords on disk
    [phi,rho] = cart2pol(x,y);
    rho(rho > 1) = 1;                  % clamp outside disk
    theta = acos(sqrt(1-rho.^2));      % map radius -> polar angle
    ct = cos(theta);

    cols = {};
    ells = [];

    for l = 0:order
        % Associated Legendre values at cos(theta)
        Pl_all = legendre(l, ct, 'sch');
        for m = -l:l
            if safe && (m==0 && l>0)
                continue  % drop all axis-symmetric terms except constant
            end

            Pm = Pl_all(abs(m)+1,:).'; 
            if m < 0
                Y = sqrt(2) * Pm .* sin(abs(m)*phi);
            elseif m > 0
                Y = sqrt(2) * Pm .* cos(m*phi);
            else
                Y = Pm;
            end

            cols{end+1} = Y; 
            ells(end+1,1) = l;
        end
    end

    A = cat(2, cols{:});
end