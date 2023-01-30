% gunzip('https://archive.ics.uci.edu/ml/machine-learning-databases/covtype/covtype.data.gz');
% load covtype.data % covtype 581012 x 55, class = [1,2,3,4,5,6,7]
% clear; clc; close all;
function [] = test_covtype_pdhg(para_p, delta)

addpath('subfunc');
load('data/covtype/covtype.mat', 'covtype');
random_seed = 0;
rng(random_seed);

covertype_cell = cell(1, 7);
cnt = zeros(1, 7);
for i = 1:7
    idx = find(covtype(:, end) == i);
    covertype_cell{i} = covtype(idx, 1:end-1);
    % num of samples in each class: 211840, 283301, 35754, 2747, 9493, 17367, 20510
    cnt(i) = length(idx);
end

%%
careclass = [4, 5];
data = [covertype_cell{careclass(1)}; covertype_cell{careclass(2)}];
N_0 = cnt(careclass(1));
N_1 = cnt(careclass(2));
N = N_0 + N_1;
y = zeros(1, N); 
y(N_0+1:end) = 1;

% for test
% idx_sample = randsample(N, round(N*0.01));
% idx1 = randsample(N_0, 500);
% idx2 = randsample(N_1, 500) + N_0;
% idx_sample = [idx1; idx2];
% data = data(idx_sample, :);
% y = y(idx_sample);

n_v = length(y);

%%
% para_p = 0.4; 
% delta = 0.3162;
mode_str = 'l';
mode_num = 2;
scale_B = 0.9;

% folder = 'data_pdhg/covtype/type45_rng0_bins20_ww_out3_in2_out1_in52/';
% file = strcat(folder, mode_str, '_B', string(scale_B), '_p', string(para_p), '_delta', string(delta), '/');
% mkdir(file)
% file_log = strcat(file, 'log.txt');
% file_log_id = fopen(file_log, 'w');
% fclose(file_log_id);
% diary(file_log);

%%
n_bins = 20;
p = linspace(0, 1, n_bins+1);
quat_flag = 10; % the first 10 features are numerical, the rest are categorical
quantile_th = zeros(quat_flag, n_bins+1);

n_e = 0;
R = zeros(n_e, n_v);
for i = 1:quat_flag
    quantile_th(i, :) = quantile(data(:,i), p);
    quantile_th(i, 1) = quantile_th(i, 1) - 1;
    for j = 1:n_bins
        temp_list = find(data(:,i) <= quantile_th(i,j+1) & data(:,i) > quantile_th(i,j)); 
        if length(temp_list) > 1 && length(temp_list) < n_v % consider |e| >= 2
            n_e = n_e + 1;
            center = median(data(temp_list,i));
            dist = abs(data(temp_list,i)' - center);
            if all(dist == 0)
                R(n_e, temp_list) = exp(-para_p * dist);
            else
                R(n_e, temp_list) = exp(-para_p * dist / max(dist)); 
            end
        end
    end
end

% kappa: hyperedge weights in EDVWs-based hypergraphs
% R: EDVWs in EDVWs-based hypergraphs
% mu: vertex weights
[n_e, n_v, ~, ~, R, incidence_list, parameter_list, mu] = hg_para(R, para_p, mode_str, delta, 'w', 'covtype');

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
tic
warmstart = labels_rwc;

[labels, NCut, vmin] = reducible_hypergraph_partition_pdhg(incidence_list, parameter_list, mu, ...
    n_v, n_e, mode_num, delta, dec_outloop, err_inloop, warmstart, ...
    B, edge_weights, iy, jy, n, M, maxiter_inner_bound, maxiter_inner, maxiter_outer);

err = comp_err(labels, y);

all_err_ncc(3, 1) = err;
all_err_ncc(3, 2) = NCut;
all_labels(3, :) = labels;
all_eigvec(3, :) = vmin;

fprintf('proposed - 2 %f %f\n', NCut, err);
toc

%%
for i = 1:10
    tic
    warmstart = randn(1, n_v);
    random_start(i, :) = warmstart;

    [labels, NCut, vmin] = reducible_hypergraph_partition_pdhg(incidence_list, parameter_list, mu, ...
        n_v, n_e, mode_num, delta, dec_outloop, err_inloop, warmstart, ...
        B, edge_weights, iy, jy, n, M, maxiter_inner_bound, maxiter_inner, maxiter_outer);

    err = comp_err(labels, y);

    all_err_ncc(i+3, 1) = err;
    all_err_ncc(i+3, 2) = NCut;
    all_labels(i+3, :) = labels;
    all_eigvec(i+3, :) = vmin;
    
    fprintf('proposed - random %f %f\n', NCut, err);
    toc
end

%% save results

diary off

save(strcat(file, 'all_err_ncc.mat'), 'all_err_ncc');
save(strcat(file, 'all_labels.mat'), 'all_labels');
save(strcat(file, 'all_eigvec.mat'), 'all_eigvec');
save(strcat(file, 'random_start.mat'), 'random_start');

end
