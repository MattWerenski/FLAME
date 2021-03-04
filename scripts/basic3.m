clear;clc;
addpath(genpath('code'));
addpath(genpath('data'));
%addpath(genpath('libsvm/matlab'));
%addpath(genpath('mtba'));

%% Example parameters

% use human or yeast data
options.org = 'yeast';

% which type of annotations to use
% options: {bp, mf, cc} for human GO,
%          {level1, level2, level3} for yeast MIPS
options.onttype = 'bp'; 

% consider terms in a specific size range (GO only)
% examples: [11 30], [31 100], [101 300]
options.ontsize = [31 100];       

% number of kNN
k=5;

% use 1-p percent of the training labels to form extra graph
extra_graph_filt = 0.0;

% Number of bi-clusters to create (-1 to not bi-cluster)
options.num_clusters = 4; 



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
options.embedding.cannotlink_penalty = 64; 

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

% compute the network walks for both mu and f.
walks = compute_rwr(network_files, ngene, -1, options);

% mu embedding is shared across all folds so we take it before cross-val
x_mu = svd_embed(walks, options.embedding.ndim);

% mashup results
acc_mu = zeros(length(folds), 1);
f1_mu = zeros(length(folds), 1);
auc_mu = zeros(length(folds), 1);

% flame results
acc_f = zeros(length(folds), 1);
f1_f = zeros(length(folds), 1);
auc_f = zeros(length(folds), 1);

weighted = true;
fprintf('weighted: true \n');

for i = 1:length(folds)
    fprintf('Fold %d / %d \n', i, options.kfolds);
    fprintf('Using k = %f \n', k);
    
    train_filt = folds(i).train_filt;
    test_filt = folds(i).test_filt;
     
    training_labels = anno.*(train_filt.');
    [l1, l2] = size(training_labels);
    
    rf = rand(l1,l2) > extra_graph_filt; % use 1-p fraction of the labels
    filtered_labels = (rf .* training_labels) > 0;

    % copy the original walk matrices
    walks2 = walks;

    % get the additional label matrix
    link_mat = mustlink(filtered_labels);
    restart_prob=0.5;
    c_walk = constraint_walk(link_mat,restart_prob);
    % add it to the network walks 
    walks2(end+1,:,:) = c_walk;
    
    % compute the embedding with added network
    x_f = svd_embed(walks2, options.embedding.ndim);
     
    
    %% Performs Mashup for comparison
    fprintf('[Performing mu version]\n');
    [dist_mat,knn] = compute_knn_labelled(x_mu, k, train_filt);
    [acc, f1, auc] = matrix_majority_voting(anno, test_filt,train_filt, knn, dist_mat, weighted);
    acc_mu(i) = acc;
    f1_mu(i) = f1;
    auc_mu(i) = auc;
    
    %% Perform flame version
    fprintf('[Perfoming f version]\n');
    [dist_mat,knn] = compute_knn_labelled(x_f, k, train_filt);
    [acc, f1, auc] = matrix_majority_voting(anno, test_filt,train_filt, knn, dist_mat, weighted);
    acc_f(i) = acc;
    f1_f(i) = f1;
    auc_f(i) = auc;
end


fprintf('[mu accuracy mean = %f ]\n', mean(acc_mu));
fprintf('[mu accuracy std = %f ]\n', std(acc_mu));
fprintf('[mu f1 mean = %f ]\n', mean(f1_mu));
fprintf('[mu f1 std = %f ]\n', std(f1_mu));
fprintf('[mu auc mean = %f ]\n', mean(auc_mu));
fprintf('[mu auc std = %f ]\n', std(auc_mu));

fprintf('[f accuracy mean = %f ]\n', mean(acc_f));
fprintf('[f accuracy std = %f ]\n', std(acc_f));
fprintf('[f f1 mean = %f ]\n', mean(f1_f));
fprintf('[f f1 std = %f ]\n', std(f1_f));
fprintf('[f auc mean = %f ]\n', mean(auc_f));
fprintf('[f auc std = %f ]\n', std(auc_f));
 


    
     
