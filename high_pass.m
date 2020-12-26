function high_pass(WC) %RL circuit

    C = 10 ^ -6;
    R = 1 / (WC * C);
    
    disp("Values of components"); 
    R, C
    
    w = linspace(0, 9 * WC, 10000);
    TF =@(w, R, C) (1i * R * C .* w) ./ ((1i * R * C .* w) + 1) ; %Across R
    r = abs(TF(w, R, C));
    
    figure;
    im = imread('RC_HIGH_PASS.jpg');
    imshow(im);
    
    figure;
    DB = 20 .* log10(r);
    plot(w, DB, 'LineWidth', 1.5); hold on; grid on;
    
    title("Plot of gain in DB");
    xlabel("---------------> Angular Frequency(\omega)");
    ylabel("---------------> Gain in DB Across R");
    peak = max(DB);
    idxC = cutoff(DB, peak - 3, 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    plot(w(idxC(1:end)), DB(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    figure;
    plot(w, r, 'LineWidth', 1.5); hold on; %PLOTTING THE TRANSFER FUNCTION
    grid on;
    title("Plot of the transfer function");
    xlabel("---------------> Angular Frequency(\omega)")
    ylabel("---------------> The Transfer Function Across R")
    
    idxC = cutoff(r, r(length(r)) / sqrt(2), 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    disp("CUTOFF FREQ / HALF POWER FREQ: ");
    disp(w(idxC(1:end)));
    
    plot(w(idxC(1:end)), r(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    hold on;
    
    high_pass_simulation(R, C, TF, WC);
end