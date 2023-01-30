% clear; clc; close all;
function [] = test_rcv1_pdhg(para_p, delta)

addpath('subfunc');
load('data/rcv1/data_C14_C23.mat', 'R', 'y'); % R: |E|x|V|, y: 1x|V|
random_seed = 0;
rng(random_seed);

%% balance
% idx0 = find(y == 0);
% idx1 = find(y == 1);
% num_ones = sum(y);
% idx0 = idx0(1:num_ones);
% idx = [idx0, idx1];
% R = R(:, idx);
% y = y(idx);

%%
mode_str = 'l';
mode_num = 2;
scale_B = 0.9;

folder = 'data_pdhg/rcv1_C14_C23/rng0_ww_out3_in2_out1_in52/';
file = strcat(folder, mode_str, '_B', string(scale_B), '_p', string(para_p), '_delta', string(delta), '/');
mkdir(file)
file_log = strcat(file, 'log.txt');
file_log_id = fopen(file_log, 'w');
fclose(file_log_id);
diary(file_log);

[n_e, n_v, ~, ~, R, incidence_list, parameter_list, mu] = hg_para(R, para_p, mode_str, delta, 'w', 'text');

all_err_ncc = zeros(13, 2);
all_labels = zeros(13, n_v);
all_eigvec = zeros(13, n_v);
random_start = zeros(10, n_v);

%% rw-based clique
tic 
f_rwc = hg_rw_laplacian('c');
[~, eigvec_rwc, eigval_rwc] = f_rwc(R, 'std', 2);
vmin_rwc = eigvec_rwc(:, 2)';
fprintf('%f %f\n', eigval_rwc);

[labels_rwc, NCut_rwc] = general_optthreshold(incidence_list, parameter_list, mu, n_v, n_e, mode_num, delta, vmin_rwc);
err_rwc = comp_err(labels_rwc, y);

all_err_ncc(1, 1) = err_rwc;
all_err_ncc(1, 2) = NCut_rwc;
all_labels(1, :) = labels_rwc;
all_eigvec(1, :) = vmin_rwc;

fprintf('baseline %f %f\n', NCut_rwc, err_rwc);
toc

%% submodular
f_sm = hg_expansion(mode_str);
W = f_sm(incidence_list, parameter_list, n_v, delta);
[B, edge_weights, iy, jy, n, M] = find_W_pdhg(W);

B = B/scale_B;

dec_outloop = 1e-3;
err_inloop = 1e-2;
maxiter_outer = 10;
maxiter_inner_bound = 100000;
maxiter_inner = 100;

%%
tic
warmstart = vmin_rwc;

[labels, NCut, vmin] = reducible_hypergraph_partition_pdhg(incidence_list, parameter_list, mu, ...
    n_v, n_e, mode_num, delta, dec_outloop, err_inloop, warmstart, ...
    B, edge_weights, iy, jy, n, M, maxiter_inner_bound, maxiter_inner, maxiter_outer);

err = comp_err(labels, y);

all_err_ncc(2, 1) = err;
all_err_ncc(2, 2) = NCut;
all_labels(2, :) = labels;
all_eigvec(2, :) = vmin;

fprintf('proposed - 1 %f %f\n', NCut, err);
toc

%%
% tic
% warmstart = labels_rwc;
% 
% [labels, NCut, vmin] = reducible_hypergraph_partition_pdhg(incidence_list, parameter_list, mu, ...
%     n_v, n_e, mode_num, delta, dec_outloop, err_inloop, warmstart, ...
%     B, edge_weights, iy, jy, n, M, maxiter_inner_bound, maxiter_inner, maxiter_outer);
% 
% err = comp_err(labels, y);
% 
% all_err_ncc(3, 1) = err;
% all_err_ncc(3, 2) = NCut;
% all_labels(3, :) = labels;
% all_eigvec(3, :) = vmin;
% 
% fprintf('proposed - 2 %f %f\n', NCut, err);
% toc

%%
% for i = 1:10
%     tic
%     warmstart = randn(1, n_v);
%     random_start(i, :) = warmstart;
% 
%     [labels, NCut, vmin] = reducible_hypergraph_partition_pdhg(incidence_list, parameter_list, mu, ...
%         n_v, n_e, mode_num, delta, dec_outloop, err_inloop, warmstart, ...
%         B, edge_weights, iy, jy, n, M, maxiter_inner_bound, maxiter_inner, maxiter_outer);
% 
%     err = comp_err(labels, y);
% 
%     all_err_ncc(i+3, 1) = err;
%     all_err_ncc(i+3, 2) = NCut;
%     all_labels(i+3, :) = labels;
%     all_eigvec(i+3, :) = vmin;
%     
%     fprintf('proposed - random %f %f\n', NCut, err);
%     toc
% end

%% save results

diary off

save(strcat(file, 'all_err_ncc.mat'), 'all_err_ncc');
save(strcat(file, 'all_labels.mat'), 'all_labels');
save(strcat(file, 'all_eigvec.mat'), 'all_eigvec');
save(strcat(file, 'random_start.mat'), 'random_start');

end
