function [L, W_triu, W_tril, ix, jx, n, m] = find_W(W)
    % W: weighted adjacency matrix of a digraph
    n = size(W, 1); % number of nodes in the reduced graph
    W2 = W + W';
    L = 4 * max(sum(W2.^2)); % upper bound on Lipschitz constant
    W2 = triu(W2, 1);
    [ix, jx] = find(W2); % 1-index
    m = length(ix);
    W_triu = zeros(1, m);
    W_tril = zeros(1, m);
    for cnt = 1:m
        u = ix(cnt);
        v = jx(cnt);
        W_triu(cnt) = W(u, v);
        W_tril(cnt) = W(v, u);
    end
    ix = ix - 1; % change to 0-index
    jx = jx - 1;
    ix = ix'; % column vector -> row vector
    jx = jx';
end
