%%
%clear all;
close all;

% CHOOSE CHROMOSOME DATA
chr1 = 10; %choose the first chromosome
chr2 = 10; %choose the second chromosome (for interchromsomal Hi-C) or repeat the first one (for intrachromosomel Hi-C)

chromosomes = chr1:-1:chr2;
chr_list = [int2str(chromosomes(1))];

for chr = chromosomes(2:end)
    chr_list = [chr_list '-' int2str(chr)];
end

% data path
data_path = ['./matlab_data/' 'chr_' chr_list '.data.normKR.mat'];
% path to the information about chromosome's length
info_path = ['./raw_data/' 'chr_' chr_list '.data.info'];

% Load chrosomoes' Hi-C matrix KR normalized
data = load(data_path);
data = data.hicmap;
data = double(data);
% Load chromosomes' length
info = load(info_path);
%save(['./output/chr' num2str(chr1) '_dna_00per.mat'], 'data')

k = 1;

data = data(1:info(k, 1), 1:info(k, 2)); % gets a slice of the Hi-C matrix accordingly to the INFO datafile
save(['./output/chr' num2str(chr1) '.mat'], 'data');


% PLOT the HI-C map
FIG = 1;
figure(FIG);
imshow(log(data), [min(min(log(data))) max(max(log(data)))]);
axis("on", "image")
%xlabel("Position (kb)")
%ylabel("Position (kb)");
fontsize(gca, 22, "points");
%colormap(hot);
hold on;

%%
%gamma = 55;


opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["module", "node"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
significant = readtable(sprintf("../output/A1_chr%u_gamma%u_significant_0.csv", chr1, gamma), opts);
significant = table2array(significant);
insignificant = readtable(sprintf("../output/A1_chr%u_gamma%u_insignificant_0.csv", chr1, gamma), opts);
insignificant = table2array(insignificant);
clear opts

%%

c = [
  "#3182bd"
  "#e6550d"
  "#31a354"
  "#756bb1"
  "#9e9ac8"
  "#bcbddc"
  "#dadaeb"
  "#636363"
  "#969696"
  "#bdbdbd"
  "#d9d9d9"
];

%%

for k = 1:2
    if k == 1
        node_set = significant;
        line_style = "-";
        alpha = 0.5;
    else
        node_set = insignificant;
        line_style = ":";
        alpha = 0;
    end

    modules = unique(node_set(:,1));

    sizes = zeros(length(modules), 2);

    for j = 1:length(modules)
        module = modules(j);
        module_size = length(node_set(node_set(:,1) == module, 2));
        sizes(j,:) = [module module_size];
    end
    sizes = sortrows(sizes, 2, 'descend');

    modules = sizes(:,1);
    
    for j = 1:length(modules)
        module = modules(j);
        nodes = node_set(node_set(:,1) == module, 2);
    
        if j < length(c)
            color = c(j);
        else
            color = c(end);
        end

        hex_color = [hex2rgb(eraseBetween(color, 1, 1)) alpha];
    
        % find consequtive nodes
        nodes = sort(nodes);
        start = nodes(1);
        prev = start;
        for i = 2:length(nodes)
            node = nodes(i);
            if node ~= prev + 1 || node == nodes(end)
                rectangle("Position", [start start prev-start prev-start], ...
                    "FaceColor", hex_color, ...
                    "EdgeColor", color, ...
                    "LineWidth", 3,  ...
                    "LineStyle", line_style)
                start = node;
                prev = start;
            else
                prev = node;
            end
        end
    end
end

saveas(gca, sprintf("../results/chr%u_gamma%u.png", chr1, gamma));