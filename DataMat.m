% DataMat.m
% Generate input data.mat for SMODEL

clear; clc;

% -----------------------------
% File paths
% -----------------------------
adt_file   = 'adt.csv';
rna_file   = 'RNA.csv';
label_file = 'label.csv';%If available
pool_file  = 'results.csv';
spot_file  = 'Spot_location_data.csv';

out_mat = 'data.mat';

% Optional feature truncation (set Inf to keep all)
n_adt_features = Inf;
n_rna_features = Inf;

% -----------------------------
% Load data
% -----------------------------
adt = table2array(readtable(adt_file, 'ReadRowNames', true));
rna = table2array(readtable(rna_file, 'ReadRowNames', true));

n_cells = min(size(adt,1), size(rna,1));

adt_data = adt(1:n_cells, 1:min(size(adt,2), n_adt_features));
rna_data = rna(1:n_cells, 1:min(size(rna,2), n_rna_features));

X = {adt_data', rna_data'};

% Other inputs
Spot     = readmatrix(spot_file);
Y        = readmatrix(label_file);
poolsget = readmatrix(pool_file);

% -----------------------------
% Save
% -----------------------------
save(out_mat, 'X', 'Spot', 'Y', 'poolsget');
