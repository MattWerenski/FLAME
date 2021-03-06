clear;clc;
addpath(genpath('code'));
addpath(genpath('data'));
addpath(genpath('libsvm/matlab'));
addpath(genpath('mtba'));

%% Example parameters

% use human or yeast data
options.org = 'yeast';

% which type of annotations to use
% options: {bp, mf, cc} for human GO,
%          {level1, level2, level3} for yeast MIPS
options.onttype = 'bp'; 

% consider terms in a specific size range (GO only)
% examples: [11 30], [31 100], [101 300]
options.ontsize = [101 300]; 
train_ontsize1 = [11 30]; % used to load all labels in embedding
train_ontsize2 = [31 100];

% number of kNN
fprintf('number of k nearest neighbors: \n ');
k1=10
k2=1

% use 1-p percent of the training labels to form extra graph
%extra_graph_filt = 0.5;

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
[genes, ngene, anno, golevels] = load_anno(options);

options.ontsize = train_ontsize1;
[~, ~, anno_train, golevels1] = load_anno(options);
options.ontsize = train_ontsize2;
[~, ~, anno_train2, golevels2] = load_anno(options);

train_anno = vertcat(anno,anno_train,anno_train2);
train_golevels = vertcat( golevels, golevels1, golevels2);
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
    %exit
    walks = compute_rwr(network_files, ngene, -1, options);
    %save('walks.mat', 'walks', '-v7.3');
end

% Create the standard Mashup embedding
x_mash = svd_embed(walks, options.embedding.ndim);
% and get the pairwise distance for that embedding
dist_mat_mash = squareform(pdist(x_mash'));


acc_flame = zeros(length(folds), 1);
acc_mash = zeros(length(folds), 1);

for i = 1:length(folds)
    
    walks2 = walks;
    
    train_filt = folds(i).train_filt;
    test_filt = folds(i).test_filt;
    
    training_labels = train_anno.*(train_filt.');
    [l1, l2] = size(training_labels);
    
    rand_filt = rand(l1,l2);
    rf = rand_filt > 0; % use 1-p fraction of the labels
    filtered_labels = (rf .* training_labels) > 0;
    
    
    link_mat = mustlink_weighted(filtered_labels, train_golevels);
    
    
    restart_prob=0.5;
    c_walk = constraint_walk(link_mat,restart_prob);
    walks2(end+1,:,:) = link_mat;
    
    x_flame = svd_embed(walks2, options.embedding.ndim);
    
    %% SMashup integration
    fprintf('[SMashup] Fold %d / %d \n', i, options.kfolds);
    
    %% Performs the base Mashup for comparison
    fprintf('[Performing base version]\n');
    [acc, predictions] = majority_voting_k1_k2(dist_mat_mash,k1,k2, anno, test_filt,train_filt);
    acc_mash(i) = acc;

    %% Perform supervised version
    fprintf('[Perfoming supervised version]\n');

    dist_mat_flame = squareform(pdist(x_flame'));
    acc_flame(i) = majority_voting_backup(dist_mat_flame,k1, predictions, anno, test_filt,train_filt);
 
end

fprintf('[Mashup accuracy mean = %f ]\n', mean(acc_mash));
fprintf('[Mashup accuracy std = %f ]\n', std(acc_mash));

fprintf('[Flame accuracy mean = %f ]\n', mean(acc_flame));
fprintf('[Flame accuracy std = %f ]\n', std(acc_flame));

    
     
