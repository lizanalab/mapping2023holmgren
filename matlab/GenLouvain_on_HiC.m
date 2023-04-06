clear all;
close all;

% Add helper functions
addpath('HelperFunctions');
addpath('GL_output');


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
data_path = ['C:\Users\Dolores_B\Documents\PhD\HiC\python\matlab_data\' 'chr_' ...
  chr_list '.data.normKR.mat'];
% path to the information about chromosome's length
info_path = ['C:\Users\Dolores_B\Documents\PhD\HiC\python\raw_data\' 'chr_' ...
  chr_list '.data.info'];

% Load chrosomoes' Hi-C matrix KR normalized
data = load(data_path);
data = data.hicmap;
data = double(data);
% Load chromosomes' length
info = load(info_path);

% SAVE CHR Hi-C map to the same folder where I will have community assignments
save(['C:\Users\Dolores_B\Documents\PhD\HiC\matlab\GL_output\chr' num2str(chr1) '_dna_00per.mat'], 'data')
fprintf('Data has been loaded for chr%s\n', num2str(chr1))

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


%for alpha = [0.75, 0.85, 0.9, 0.95, 1, 1.05, 1.08, 1.1, 1.15, 1.2, 1.25, 1.3] % VARY the exponent
alpha = 1.08; % FIX one parameter

for gamma = [0.66] % [0.54, 0.6, 0.66, 0.85, 0.87] %[0.4, 0.5, 0.6, 0.7, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2]
   for iterN = linspace(1,5,5) % Repeat GenLouvain algorithm 5 times and save results of each iteration        
        S = []; % to save community assignemnts for nodes

        for k = 1:size(info, 1)          

          chromosome = info(k, 3);
          
          if (k == 1) % if gL will be done on a single chromosome (OUR CASE)
            data_part = data(1:info(k, 1), 1:info(k, 2)); % gets a slice of the Hi-C matrix accordingly to the INFO datafile
            
            % to remove self-loops (?)
            % data_part = data_part - diag(diag(data_part));

            % to calculate the distance matrix between node i and node j;
            % Dij = |i-j| -- this should be the way how D matrix is calculated
            D = toeplitz(1-1:info(k,1)-1); % It makes the D matrix to have zeros on the main diagonal, all following diagonals are +1
          else % (NEVER HAVE BEEN IN THIS CASE)
            data_part = data((info(k-1, 1) + 1):info(k, 1), (info(k-1, 2) + 1):info(k, 2));
            D = toeplitz(1:(info(k,1) - info(k-1, 1)));
          end

          % Calculare the modularity matrix
          B = fractalglobulepoly_f(data_part, D, gamma, alpha);
          
          % Get communities with genlouvain
          [S_p, Q] = iterated_genlouvain(B);
          fprintf('Modularity. Partition %s', num2str(iterN))
          fprintf(':   %s\n', num2str(Q))
         
          % Save results into S matrix
          S = [S; S_p ones(length(S_p), 1) * chromosome ones(length(S_p), 1) * Q];          

          gamma_n = num2str(gamma *100);
          save(['C:\Users\Dolores_B\Documents\PhD\HiC\matlab\GL_output\A1_chr' num2str(chr1) '_gamma' gamma_n '_dna00per' '_iter' num2str(iterN) '.mat'], 'S')
  
        end
    end
end