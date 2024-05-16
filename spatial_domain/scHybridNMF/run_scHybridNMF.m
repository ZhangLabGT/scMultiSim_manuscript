addpath('scHybridNMF-main/'); 

expr = readmatrix("expression_1.txt");
locs = readmatrix("meta_1.txt");

expr = expr';
locs = locs(:, 2:3)';

% A: m x n
% L: 2 x n

[H_A, W_A, nmf_ids, j_ids] = scHybridNMF(expr, locs, 3, 0.2, 0.01, 500);

[~, C] = max(H_A);

filename = 'C.csv';
csvwrite(filename, C.');
