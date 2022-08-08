function func = hg_expansion(mode)
    if strcmp(mode, 'c')
        func = @clique_expansion;
    elseif strcmp(mode, 's')
        func = @star_expansion;
    else 
        func = @lawler_expansion;
    end
end

function A = clique_expansion(incidence_list, parameter_list, n_v, ~)
    n_e = length(incidence_list);
    A = zeros(n_v);
    for e_i = 1:n_e
        nodes = incidence_list{e_i};
        EDVWs = parameter_list{e_i};
        esize = length(nodes);
        for i = 1:esize
            vi = nodes(i);
            for j = i+1:esize
                vj = nodes(j);
                A(vi, vj) = A(vi, vj) + EDVWs(i) * EDVWs(j);
                A(vj, vi) = A(vi, vj);
            end
        end
    end
end

function A = star_expansion(incidence_list, parameter_list, n_v, ~)
    n_e = length(incidence_list);
    A = zeros(n_v + n_e);
    for e_i = 1:n_e
        nodes = incidence_list{e_i};
        EDVWs = parameter_list{e_i};
        v_e = n_v + e_i;
        for i = 1:length(nodes)
            v = nodes(i);
            A(v, v_e) = EDVWs(i);
            A(v_e, v) = A(v, v_e);
        end
    end
end

function A = lawler_expansion(incidence_list, parameter_list, n_v, delta)
    n_e = length(incidence_list);
    A = zeros(n_v + 2 * n_e);
    for e_i = 1:n_e
        nodes = incidence_list{e_i};
        EDVWs = parameter_list{e_i};
        e1 = n_v + e_i * 2 - 1;
        e2 = n_v + e_i * 2;
        A(e1, e2) = delta * sum(EDVWs);
        for i = 1:length(nodes)
            v = nodes(i);
            A(v, e1) = EDVWs(i);
            A(e2, v) = EDVWs(i);
        end
    end
end
