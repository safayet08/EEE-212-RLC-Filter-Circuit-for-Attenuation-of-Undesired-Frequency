function low_pass(WC) %RL circuit

    L = 2.7 * 10 ^ -3;
    R = WC * L;
    
    disp("Values of components"); 
    R, L
    
    w = linspace(0, 9 * WC, 10000);
    TF =@(w, R, L) (R) ./ (1i*L.* w + R); %Across L
    r = abs(TF(w, R, L));
    
    figure;
    im = imread('RL_HIGH_PASS.jpg');
    imshow(im);
    
        
    figure;
    DB = 20 .* log10(r);
    plot(w, DB, 'LineWidth', 1.5); hold on; grid on;
    
    title("Plot of gain in DB");
    xlabel("---------------> Angular Frequency(\omega)");
    ylabel("---------------> Gain in DB Across L");
    peak = max(DB);
    idxC = cutoff(DB, peak - 3, 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    plot(w(idxC(1:end)), DB(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    

    
    figure;
    plot(w, r, 'LineWidth', 1.5); hold on; %PLOTTING THE TRANSFER FUNCTION
    grid on;
    title("Plot of the transfer function");
    xlabel("---------------> Angular Frequency(\omega)")
    ylabel("---------------> The Transfer Function Across L")

    idxC = cutoff(r, 1 / sqrt(2), 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    disp("CUTOFF FREQ / HALF POWER FREQ: ");
    disp(w(idxC(1:end)));
    
    plot(w(idxC(1:end)), r(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    hold on;
    low_pass_simulation(R, L, TF, WC);
end