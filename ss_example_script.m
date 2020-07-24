addpath code

%% Required external dependencies: see README.txt for more information
addpath ../libsvm-3.24/matlab % LIBSVM package (for cross_validation.m)
addpath ../L-BFGS-B-C/Matlab % L-BFGS package (only if svd_approx = false)

%% Example parameters
test_fraction = 0.5 % portion of the labelled vertices to go in the 
                    % testing set. (1 - test_fraction) is used to train
org = 'yeast';      % use human or yeast data
onttype = 'level1'; % which type of annotations to use
                    %   options: {bp, mf, cc} for human GO,
                    %            {level1, level2, level3} for yeast MIPS
ontsize = [];       % consider terms in a specific size range (*human GO only*)
                    %   examples: [11 30], [31 100], [101 300]
svd_approx = false;  % use SVD approximation for Mashup
                    %   recommended: true for human, false for yeast
ndim = 500;         % number of dimensions
                    %   recommended: 800 for human, 500 for yeast
restart_prob = 0.5; % chance that the random walk restarts itself
walk_mode = 'unsupervised'; % determines what type of RWR to perform
                    % 'unsupervised' is the default in mashup
                    % 'semi-supervised' uses "class teleportation" allowing escape
                    % 'supervised' uses class teleportation without escape
teleport_prob= 0.5; % only uses in semi-supervised walk mode. The chance that
                    % the chance that a walk entering a labeled vertex jumps to
                    % another vertex with that label instead of doing the normal walk
embedding_mode = 'unsupervised'; % determines what type of embedding to do
                    % Only applies when svd_approx is set to false. Options are
                    % 'unsupervised' performs the standard optimization of mashup
                    % 'supervised' introduces penalties using the labels
mustlink_penalty = 0.1; % in supervised embedding the amount of penalty placed
                    % on the mustlink constraints
cannotlink_penalty = 0.01; % in supervised embedding the amount of penalty placed
                    % on the cannot link constraints



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



%% Load known annotations
fprintf('[Loading annotations]\n');
if strcmp(org, 'human')
  anno = load_go(onttype, genes, ontsize, true);
elseif strcmp(org, 'yeast')
  anno = load_mips(onttype, genes);
end
fprintf('Number of functional labels: %d\n', size(anno, 1));



%% Generate the training and testing sets
test_frac = 0.5;
fprintf('Acquiring test filter using %d testing fraction\n', test_frac);
[~, test_filt] = cv_partition(anno, test_frac); 
% filters proteins with no labels
label_filt = (sum(anno) > 0).'; 
% removes the unlabeled proteins from the train and test sets
train_filt =(~test_filt) & label_filt;
test_filt = test_filt & label_filt;

ntest = sum(test_filt);
ntrain = sum(train_filt);

% applies the filters as masks to the annotations
training_labels = anno.*(train_filt.');
test_labels = anno.*(train_filt.');



%% SMashup integration
fprintf('[SMashup]\n');

%% Performs the specified variant of RWR
fprintf('[Performing RWR step]\n');
if strcmp(walk_mode, 'unsupervised')
  walks = unsupervised_rwr(network_files, ngene, restart_prob);
elseif strcmp(walk_mode, 'semi-supervised')
  walks = semisupervised_rwr(network_files, ngene, restart_prob, ...
    teleport_fraction, training_labels);
else
  walks = supervised_rwr(network_files, ngene, restart_prob, ...
    1.0, training_labels);
end

fprintf('[Performing embedding step]\n');
if svd_approx
  x = svd_embed(walks, ndim);
else
  if strcmp(embedding_mode, 'unsupervised')
    x = unsupervised_embed(walks, ndim);
  else
    x = supervised_embed(walks, ndim, mustlink_penalty, cannotlink_penalty);
  end
end



%% Split annotations 
function [ntest, test_filt] = cv_partition(anno, test_frac)
    ng = size(anno, 2);
    ntest = floor(ng * test_frac);
    test_ind = randperm(ng, ntest);
    test_filt = false(ng, 1);
    test_filt(test_ind) = true;
    ntest = length(test_ind);
end