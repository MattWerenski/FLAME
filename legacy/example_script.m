addpath code

%% Required external dependencies: see README.txt for more information
%addpath ../libsvm-3.24/matlab % LIBSVM package (for cross_validation.m)
%addpath ../L-BFGS-B-C/Matlab % L-BFGS package (only if svd_approx = false)

%% Example parameters
org = 'yeast';      % use human or yeast data
onttype = 'mf'; % which type of annotations to use
                    %   options: {bp, mf, cc} for human or yeast GO,
                    %            {level1, level2, level3} for yeast MIPS
ontsize = [101 300];       % consider terms in a specific size range (*human GO only*)
                    %   examples: [11 30], [31 100], [101 300]
nperm = 5;          % number of cross-validation trials
svd_approx = false;  % use SVD approximation for Mashup
                    %   recommended: true for human, false for yeast
ndim = 500;         % number of dimensions
                    %   recommended: 800 for human, 500 for yeast

%% Construct network file paths 
string_nets = {'neighborhood', 'fusion', 'cooccurence', 'coexpression', ...
               'experimental', 'database'};
network_files = cell(1, length(string_nets));
for i = 1:length(string_nets)
  network_files{i} = sprintf('data/networks/%s/%s_string_%s_adjacency.txt', ...
                             org, org, string_nets{i});
end

%% Load gene list
gene_file = sprintf('data/networks/%s/%s_string_genes.txt', org, org);
genes = textread(gene_file, '%s');
ngene = length(genes);

%% Mashup integration
fprintf('[Mashup]\n');
%x = mashup(network_files, ngene, ndim, svd_approx);

%% Load known annotations
fprintf('[Loading annotations]\n');
if strcmp(org, 'human')
  anno = load_go(onttype, genes, ontsize, true);
elseif strcmp(org, 'yeast')
  anno = load_mips(onttype, genes);
end
fprintf('Number of functional labels: %d\n', size(anno, 1)); 

%% Function prediction
%fprintf('[Function prediction]\n');
%[acc, f1, aupr] = cross_validation(x, anno, nperm);

%% Output summary
%fprintf('[Performance]\n');
%fprintf('Accuracy: %f (stdev = %f)\n', mean(acc), std(acc, 1));
%fprintf('F1: %f (stdev = %f)\n', mean(f1), std(f1, 1));
%fprintf('AUPR: %f (stdev = %f)\n', mean(aupr), std(aupr, 1));
