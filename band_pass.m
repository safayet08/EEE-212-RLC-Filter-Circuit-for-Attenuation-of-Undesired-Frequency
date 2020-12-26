
function band_pass(WC1, WC2) 
    L = 2.7 * 10 ^ -3;
    BW = abs(WC2 - WC1);
    R = BW * L;
    W0 = sqrt(WC1 * WC2);
    disp("Quality Factor: ");
    Q = W0 / BW
    C = 1 / (W0^2 * L);
    
    disp("Values of components"); 
    R, L, C
    
    w = linspace(max(0, -BW + WC1), WC2 + 2*BW, 10000);
    TF =@(w, R, L, C) (1i*R*C .* w) ./ (1i*R*C .* w - w .* w .* L*C + 1); %Across R
    
    r = abs(TF(w, R, L, C));
    
    figure;
    im = imread('BAND_PASS.jpg');
    imshow(im);
    
    figure;
    DB = 20 .* log10(r);
    plot(w, DB, 'LineWidth', 1.5); hold on; grid on;
    
    [p, idx0] = findpeaks(DB);
    disp("CENTER ANGULAR FREQ: "); disp(w(idx0));  
    plot(w(idx0), DB(idx0), 'ro', 'LineWidth', 3); hold on; %PEAK POINT PLOTTING
    
    idxC = cutoff(DB, DB(idx0(1)) - 3, 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    disp("CUTOFF ANGULAR FREQ / HALF POWER ANGULAR FREQ: ");
    disp(w(idxC(1:end)));
    plot(w(idxC(1:end)), DB(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    xlabel("---------------> Angular Frequency(\omega)")
    ylabel("Gain in DB")
    title("Plot of gain in DB Across R")
    
    figure;
    plot(w, r, 'LineWidth', 1.5); hold on; grid on; %PLOTTING THE TRANSFER FUNCTION
    title("Plot of the transfer function");
    xlabel("---------------> Angular Frequency(\omega)")
    
    [p, idx0] = findpeaks(r);
    disp("CENTER ANGULAR FREQ: "); disp(w(idx0));  
    plot(w(idx0), r(idx0), 'ro', 'LineWidth', 3); hold on; %PEAK POINT PLOTTING
    
    idxC = cutoff(r, r(idx0(1)) / sqrt(2), 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    disp("CUTOFF ANGULAR FREQ / HALF POWER ANGULAR FREQ: ");
    disp(w(idxC(1:end)));
    
    plot(w(idxC(1:end)), r(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    band_pass_simulation(R, L, C, TF, WC1, WC2);
end