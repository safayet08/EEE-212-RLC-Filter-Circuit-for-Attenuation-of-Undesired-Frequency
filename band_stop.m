function band_stop(WC1, WC2) 

    L = 2.7 * 10 ^ -3;
    BW = abs(WC2 - WC1);
    R = BW * L;
    W0 = sqrt(WC1 * WC2);
    disp("Quality Factor: ");
    Q = W0 / BW
    C = 1 / (W0^2 * L);
    
    disp("Values of components"); 
    R, L, C
    
    w = linspace(0, WC2 + BW, 10000);
    TF =@(w, R, L, C) (1 - w.*w.*L*C) ./ sqrt((1 - w .* w .* L*C) + w .* w.* (R*C)^ 2); %ACROSS LC
    
    r = abs(TF(w, R, L, C));
    
    plot(w, r);
end