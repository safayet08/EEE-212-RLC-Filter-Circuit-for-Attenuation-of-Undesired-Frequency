function band_pass_simulation(R, L, C, TF, WC1, WC2) 
    sinusoid =@(w, t, phi) sin(w * t + phi);
    
%   For a series of wt, 2wt, 3wt, 4wt, ....., nwt; set base here
    base = 500;
    upto = 10;
    
    t = linspace(0, 2.5 * (2 * pi) / (base), 10000);
    
    total_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        v_out = TF(w, R, L, C);
        response = abs(v_out) * sinusoid(w, t, 0);
        total_response = total_response + response;
    end
    
    pass_band_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        if w >= WC1 && w <= WC2
            v_out = TF(w, R, L, C);
            response = abs(v_out) * sinusoid(w, t, 0);
            pass_band_response = pass_band_response + response;
        end
    end
    
    
    figure;
    plot(t, total_response, 'LineWidth', 1.5); hold on; grid on;
    plot(t, pass_band_response, 'LineWidth', 1.5); hold on;
    legend("Total Response", "Pass Band Response");
    xlabel("---------------> Time (t)");
    ylabel("---------------> Voltage across R (V_{out})");
end