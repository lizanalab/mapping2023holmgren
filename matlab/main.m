clear all;
close all;

% Add helper functions
addpath('HelperFunctions');
addpath('../GenLouvain');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE CHROMOSOME DATA
chr1 = 10; %choose the first chromosome
chr2 = 10; %choose the second chromosome (for interchromsomal Hi-C) or repeat the first one (for intrachromosomel Hi-C)

chromosomes = chr1:-1:chr2;
chr_list = [int2str(chromosomes(1))];

for chr = chromosomes(2:end)
    chr_list = [chr_list '-' int2str(chr)];
end

% data path
data_path = ['./matlab_data/' 'chr_' ...
    chr_list '.data.normKR.mat'];
% path to the information about chromosome's length
info_path = ['./raw_data/' 'chr_' ...
    chr_list '.data.info'];

% Load chrosomoes' Hi-C matrix KR normalized
data = load(data_path);
data = data.hicmap;
data = double(data);
% Load chromosomes' length
info = load(info_path);

% SAVE CHR Hi-C map to the same folder where I will have community assignments
%save(['./output/chr' num2str(chr1) '_dna_00per.mat'], 'data')
%fprintf('Data has been loaded for chr%s\n', num2str(chr1))

%thresh = 5e-3;
%data(data < thresh) = 0;

% PLOT the HI-C map
FIG = 1;
figure(FIG);
imshow(log(data), [min(min(log(data))) max(max(log(data)))]);
colormap(hot);
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenLouvain community partition

% we set three parameters for the null-model:
% ---> alpha - the exponent of the decay of contact frequencies between node i and node j
% ---> gamma - the resolution paramater
% ---> distance matrix D

%alpha = 1.08;
%alphas = 0.75;
alphas = 1.27;

k = 1;
chromosome = info(k, 3);

data_part = data(1:info(k, 1), 1:info(k, 2)); % gets a slice of the Hi-C matrix accordingly to the INFO datafile

%%
%imshow(log(data_part), [min(min(log(data_part))) max(max(log(data_part)))]);
%colormap(hot);

%data_part = data_part - diag(diag(data_part)); % to remove self-loops

% to calculate the distance matrix between node i and node j;
% Dij = |i-j| -- this should be the way how D matrix is calculated
D = toeplitz(1-1:info(k,1)-1); % It makes the D matrix to have zeros on the main diagonal, all following diagonals are +1

[num_nodes, ~] = size(data_part);
num_iters = 100;

gammas = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2];

for alpha = alphas
    Qs = zeros(size(gammas, 2), num_iters);
    
    for i = 1:size(gammas, 2)
        gamma = gammas(i);
    
        % Calculare the modularity matrix
        [B, twom] = fractalglobulepoly_f(data_part, D, gamma, alpha);
        S = zeros(num_nodes + 1, num_iters);
    
        tic;
        for seed = 1:num_iters
            [S_p, Q, n_it] = iterated_genlouvain(B,10000,0,0,'moverandw',[],[],seed);
            fprintf('Partition %u: Q=%.2d gamma=%.2d (%u iters)\n', seed, Q, gamma, n_it)
    
            S(:, seed) = [Q / twom; S_p];
        end
        toc;
    
        S = sortrows(S.',1,'descend').';
        Qs(i,:) = S(1,:);
        S(1,:) = [];
        filename = sprintf("./output/A1_chr%u_gamma%s_alpha%s_partitions.csv", chr1, num2str(gamma*100), num2str(alpha*100));
        writematrix(S, filename, "Delimiter", "space");
        S = [[0:length(data_part)-1]' S];
        filename = sprintf("./output/A1_chr%u_gamma%s_alpha%s.csv", chr1, num2str(gamma*100), num2str(alpha*100));
        writematrix(S, filename, "Delimiter", "space");
    end
    
    writematrix(Qs, sprintf("./output/A1_chr%u_alpha%s_modularity.csv", chr1, num2str(alpha*100)), "Delimiter", "space");
end