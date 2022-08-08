function [] = test_rcv1(para_p, delta)

addpath('subfunc');
load('data/rcv1/data.mat'); % R: |E|x|V|, y: 1x|V|
random_seed = 0;
rng(random_seed);

%%
mode_str = 'l';
mode_num = 2;

folder = 'data/rcv1/rng0_ww_out3_in2_out1_in52/';
file = strcat(folder, mode_str, '_p', string(para_p), '_delta', string(delta), '/');
mkdir(file)
file_log = strcat(file, 'log.txt');
file_log_id = fopen(file_log, 'w');
fclose(file_log_id);
diary(file_log);

[n_e, n_v, card, kappa, R, incidence_list, parameter_list, mu] = hg_para(R, para_p, mode_str, delta, 'w', 'text');

%% rw-based clique
f_rwc = hg_rw_laplacian('c');
tic
[~, eigvec_rwc, eigval_rwc] = f_rwc(R, 'std', 2);
toc
vmin_rwc = eigvec_rwc(:, 2)';
fprintf('%f %f\n', eigval_rwc);

[labels_rwc, NCut_rwc] = general_optthreshold(incidence_list, parameter_list, mu, n_v, n_e, mode_num, delta, vmin_rwc);
err_rwc = comp_err(labels_rwc, y);

%% submodular
tic
f_sm = hg_expansion(mode_str);
W = f_sm(incidence_list, parameter_list, n_v, delta);
[L, W_triu, W_tril, ix, jx, n, m] = find_W(W);

dec_outloop = 1e-3;
err_inloop = 1e-2;
warmstart = vmin_rwc;
maxiter_outer = 10;
maxiter_inner_bound = 100000;
maxiter_inner = 100;

[labels, NCut, vmin] = reducible_hypergraph_partition(incidence_list, parameter_list, mu, ...
    n_v, n_e, mode_num, delta, dec_outloop, err_inloop, warmstart, ...
    L, W_triu, W_tril, ix, jx, n, m, maxiter_inner_bound, maxiter_inner, maxiter_outer);
toc

err = comp_err(labels, y);

fprintf('baseline %f %f\n', NCut_rwc, err_rwc);
fprintf('proposed %f %f\n', NCut, err);

diary off
%% save results
save(strcat(file, 'vmin_rwc.mat'), 'vmin_rwc');
save(strcat(file, 'labels_rwc.mat'), 'labels_rwc');
save(strcat(file, 'vmin.mat'), 'vmin');
save(strcat(file, 'labels.mat'), 'labels');

save(strcat(file, 'NCut_rwc.mat'), 'NCut_rwc');
save(strcat(file, 'err_rwc.mat'), 'err_rwc');
save(strcat(file, 'NCut.mat'), 'NCut');
save(strcat(file, 'err.mat'), 'err');

end
