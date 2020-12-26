function low_pass_simulation(R, L, TF, WC) 
    sinusoid =@(w, t, phi) sin(w * t + phi);
    
%   For a series of wt, 2wt, 3wt, 4wt, ....., nwt; set base here
    base = 1000;
    upto = 5;
    
    t = linspace(0, 6 * (2 * pi) / (base), 10000);
    
    total_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        v_out = TF(w, R, L);
        response = abs(v_out) * sinusoid(w, t, 0);
        total_response = total_response + response;
    end
    
    pass_band_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        if w <= WC
            w
            v_out = TF(w, R, L);
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