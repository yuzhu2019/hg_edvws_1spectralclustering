function hgw = comp_hgw(parameter_list, mode, delta)
    if mode == 'c'
        func_w_e = @w_c;
    elseif mode == 's'
        func_w_e = @w_s;
    else % mode == 'l'
        func_w_e = @w_l;
    end
    n_e = length(parameter_list);
    hgw = zeros(1, n_e);
    b_max = 0;
    for i = 1:n_e
        pi = parameter_list{i};
        %
        if mod(i, 30) == 0
            fprintf('%d %d %d\n', i, n_e, length(pi));
        end
        %
        t = sum(pi);
        s = subset_sum_closest(pi, t/2, length(pi));
        hgw(i) = func_w_e(s, t, delta);
        b_max = max(b_max, min(s, t-s)/t);
    end
    fprintf('b_max %.5f\n', b_max);
end

function res = w_c(s, t, ~)
    res = s * (t - s);
end

function res = w_s(s, t, ~)
    res = min(s, t - s);
end

function res = w_l(s, t, delta)
    res = min([s, t - s, delta * t]);
end
