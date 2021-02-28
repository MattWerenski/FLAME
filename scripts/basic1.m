addpath mtba-nix/mtba-nix/mtba-nix/mtba;
addpath libsvm-3.24/matlab; % LIBSVM package (for cross_validation.m)
addpath L-BFGS-B-C/Matlab; % L-BFGS package (only if svd_approx = false)
addpath code;
addpath code/embed;
%% Example parameters

% use human or yeast data
options.org = 'yeast';

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
options.embedding.ndim = 500;

% the weight of the edges connecting dummy nodes to true nodes
options.embedding.mustlink_penalty = 1; 

% the weight of the edges connecting dummy nodes to dummy nodes
options.embedding.cannotlink_penalty = 50; 

% whether to add 1/ngene^2 strength constraint between all genes
options.embedding.use_unsupervised = false;

% when using go, whether or not to append the extra link matrix
% generated from the labels
options.walk.use_go_link = false;

% when using go and the link matrix, what fraction of links to use
options.walk.go_link_fraction = 1.0;

% chance that the random walk restarts itself
options.walk.restart_prob = 0.5;

% number of folds for k-fold cross validation
% if set to 1 or less then a single experiment is run
options.kfolds = 5;

% if options.kfolds is set to 1 or less then this is
% the portion of the labelled vertices to go in the 
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

if options.kfolds <= 1
    %% Generate training and testing sets
    fprintf('Acquiring test filter using %d testing fraction\n', options.test_fraction);
    folds = create_kfolds(anno, options);
else
    folds = create_kfolds(anno, options);
end

%% Performs the specified variant of RWR

fprintf('[Performing RWR step]\n');
if isfile(sprintf('data/walks/%s.mat', options.org))
    % loading this file automatically adds walks to the workspace.
    load(sprintf('data/walks/%s.mat', options.org));
else
    walks = compute_rwr(network_files, ngene, -1, options);
    %save('walks.mat', 'walks', '-v7.3');
end

if options.walk.use_go_link
    x_base = svd_embed(walks(1:end-1,:,:), options.embedding.ndim);
else
    x_base = svd_embed(walks, options.embedding.ndim);
end

for i = 1:length(folds)
    options.embedding.mustlink_penalty = options.embedding.mustlink_penalty * 5;
    options.embedding.cannotlink_panalty = options.embedding.cannotlink_penalty * 5;
    
    train_filt = folds(i).train_filt;
    test_filt = folds(i).test_filt;


    %% SMashup integration
%    fprintf('[SMashup] Fold %d / %d \n', i, options.kfolds);

%    fprintf('[Performing Biclustering]')
    [gene_clusters, label_clusters] = bicluster(anno, train_filt, options);

%    fprintf('[Performing embedding step]\n');
    x = compute_embedding(walks, gene_clusters, options);
    
    fprintf('[Perfoming our version]\n');
    [dist_mat, knn] = compute_knn(x,20);
    acc = majority_voting(anno, test_filt, train_filt, knn, dist_mat);
    fprintf('Acc: %f \n', acc);

    fprintf('[Performing base version]\n');
    [dist_mat, knn] = compute_knn(x_base,20);
    acc = majority_voting(anno, test_filt, train_filt, knn, dist_mat);
    fprintf('Acc: %f \n', acc);

end
