function [B, edge_weights, iy, jy, n, M] = find_W_pdhg(W)
    % W: weighted adjacency matrix of a digraph
    n = size(W, 1); % number of nodes in the reduced graph
    M = nnz(W); % number of (directed) edges in the reduced graph
    
    [iy, jy, edge_weights] = find(W);
    smatB = sparse([1:M,1:M], [iy;jy], [edge_weights;-edge_weights]);
    B = sqrt(eigs(smatB'*smatB, 1));
    iy = iy'; % column vector -> row vector
    jy = jy'; 
    edge_weights = edge_weights';
    iy = iy - 1; % change to 0-index
    jy = jy - 1;
    
%     matB = zeros(M, n); % infeasible for hugh matrix
%     edge_weights = zeros(1, M);
%     iy = zeros(1, M);
%     jy = zeros(1, M);
%     cnt = 0;
%     for u = 1:n
%         for v = 1:n
%             if W(u, v) > 0 % directed edge u -> v
%                 cnt = cnt + 1;
%                 matB(cnt, u) = W(u, v);
%                 matB(cnt, v) = -W(u, v);
%                 edge_weights(cnt) = W(u, v);
%                 iy(cnt) = u - 1; % change to 0-index
%                 jy(cnt) = v - 1;
%             end
%         end
%     end
%     B = norm(matB);
end
