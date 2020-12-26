function [idx] = cutoff(r, value, tol)
    v = abs(r - (value));
    d = v .* ones(size(r));
    minDev = min(d);
    eps = tol;
    idx = find(abs(d - minDev) <= eps);
end