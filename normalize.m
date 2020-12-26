function pos = normalize(x, idx) 
    pos = [];
    pos = idx(1);
    for i = 2: 1 : length(idx)
        if idx(i) - idx(i - 1) > 1
            pos = [pos idx(i)];
        end
    end
end