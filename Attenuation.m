function [ Xr, SNR ] = Attenuation(X, fs, r, hr, r0)
    
    f = linspace(0, fs/2, size(X,1));
    L1 = 3.2 * 10^(-7) * f.*f;
    h = hr * 10 ^ (4.6151 - 6.8346);
    muO = 0.00105;
    muN = 0.0002;
    c = 343;
    fro = 24 + 40400 * h * ((0.02 + h)/(0.391 + h));
    frn = 9 + 280 * h;
    L2 = 2*8686*muO*(f/c).*(2*(f/fro)./(1+(f/fro).^2));
    L3 = 2*8686*muN*(f/c).*(2*(f/frn)./(1+(f/frn).^2));
    L0 = (L1+L2+L3)*((r-r0)/1000);
    L4 = 20*log(r/r0)/log(10);
    L0 = L0 + L4;
    L0 = 10.^(L0./20);
    Xr = X./repmat(L0.', 1, size(X,2), size(X,3));
    SNR = mag2db(rssq(reshape(Xr, 1, [], 1))/rssq(reshape(X, 1, [], 1)));

end