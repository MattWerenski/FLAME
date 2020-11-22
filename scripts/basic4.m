addpath mtba-nix/mtba-nix/mtba-nix/mtba;
addpath libsvm-3.24/matlab; % LIBSVM package (for cross_validation.m)
addpath L-BFGS-B-C/Matlab; % L-BFGS package (only if svd_approx = false)
addpath code;
addpath code/embed;
%% Example parameters

% use human or yeast data
options.org = 'human';

% which type of annotations to use
% options: {bp, mf, cc} for human GO,
%          {level1, level2, level3} for yeast MIPS
options.onttype = 'mf'; 

% consider terms in a specific size range (GO only)
% examples: [11 30], [31 100], [101 300]
options.ontsize = [31 100];       



% Number of bi-clusters to create (-1 to not bi-cluster)
options.num_clusters = 8; 



% use SVD approximation for Mashup
% recommended: true for human, false for yeast
options.embedding.svd_approx = true;

% whether to stack the matrics really tall for svd or use the
% log transofrm and sum - needs implementation of dummy nodes
% options.embedding.svd_full = false;

% number of dimensions
% recommended: 800 for human, 500 for yeast
options.embedding.ndim = 800;

% the weight of the edges connecting dummy nodes to true nodes
options.embedding.mustlink_penalty = 1; 

% the weight of the edges connecting dummy nodes to dummy nodes
options.embedding.cannotlink_penalty = 32; 


% when using go, whether or not to append the extra link matrix
% generated from the labels
options.walk.use_go_link = false;

% chance that the random walk restarts itself
options.walk.restart_prob = 0.5;



% portion of the labelled vertices to go in the 
% testing set. (1 - test_fraction) is used to train
options.test_fraction = 0.2; 
                    
                    
%% Logs the options so we can see the parameters used in log file later
log_options(options);

 
%% Construct network file paths
string_nets = {'neighborhood', 'fusion', 'cooccurence', 'coexpression', ...
               'experimental', 'database'};
network_files = cell(1, length(string_nets));
for i = 1:length(string_nets)
     network_files{i} = sprintf('data/networks/%s/%s_string_%s_adjacency.txt', ...
         options.org, options.org, string_nets{i}); 
end

%% Load gene list
fprintf('[Loading annotations]\n');
[genes, ngene, anno] = load_anno(options);
fprintf('Number of functional labels: %d\n', size(anno, 1));

%% Generate training and testing sets
fprintf('Acquiring test filter using %d testing fraction\n', options.test_fraction);
[train_filt, test_filt, ntrain, ntest, ...
    train_labels, test_labels] = create_train_test(anno, options);

%% SMashup integration
fprintf('[SMashup]\n');

fprintf('[Performing Biclustering]')
[gene_clusters, label_clusters] = bicluster(anno, train_filt, options);

%% Performs the specified variant of RWR
fprintf('[Performing RWR step]\n');
walks = compute_rwr(network_files, ngene, train_filt, options);

fprintf('[Performing embedding step]\n');
x = compute_embedding(walks, gene_clusters, options);

%% Use the embedding with SVMs
fprintf('[Perfoming our version]\n');
run_svm(x, anno, test_filt);

%% Performs the base Mashup for comparison
fprintf('[Performing base version]\n');
if options.walk.use_go_link
    x_base = svd_embed(walks(1:end-1,:,:), options.embedding.ndim);
else
    x_base = svd_embed(walks, options.embedding.ndim);
end
run_svm(x_base, anno, test_filt);
