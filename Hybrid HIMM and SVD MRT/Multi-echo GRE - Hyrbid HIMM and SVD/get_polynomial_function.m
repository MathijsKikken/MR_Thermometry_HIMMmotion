function A = get_polynomial_function(dims, order, safe)
% Safe 2D monomial basis up to "order".
% Excludes pure even powers (x^(2k), y^0) and (x^0, y^(2k)) when safe=true.
    if nargin < 3, safe = false; end
    [yc,xc] = meshgrid(linspace(-0.5,0.5,dims(1)), linspace(-0.5,0.5,dims(2)));
    xc = xc(:); yc = yc(:);

    cols = {};
    for yp = 0:order
        for xp = 0:(order-yp)
            if safe && ( (xp>0 && yp==0 && mod(xp,2)==0) || (yp>0 && xp==0 && mod(yp,2)==0) )
                continue; % exclude x^(2k) or y^(2k)
            end
            cols{end+1} = (xc.^xp).*(yc.^yp); 
        end
    end
    A = cat(2, cols{:});
end