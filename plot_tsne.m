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
options.onttype = 'mf'; 

% consider terms in a specific size range (GO only)
% examples: [11 30], [31 100], [101 300]
options.ontsize = [11 31];       



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
options.embedding.mustlink_penalty = 32; 

% the weight of the edges connecting dummy nodes to dummy nodes
options.embedding.cannotlink_penalty = 128; 


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

%options.onttype = 'bp'; 
%[genes1, ngene1, anno1] = load_anno(options);

%anno = [anno; anno1];

fprintf('Number of functional labels: %d\n', size(anno, 1));

%% Generate training and testing sets
fprintf('Acquiring test filter using %d testing fraction\n', options.test_fraction);
[train_filt, test_filt, ntrain, ntest, ...
    train_labels, test_labels] = create_train_test(anno, options);

%% SMashup integration
fprintf('[SMashup]\n');

size(anno)
fprintf(class(anno));

fprintf('[Performing Biclustering]')
[gene_clusters, label_clusters] = bicluster(anno, train_filt, options);
%% Performs the specified variant of RWR
fprintf('[Performing RWR step]\n');
walks = compute_rwr(network_files, ngene, train_filt, options);

fprintf('[Performing embedding step]\n');
[x] = compute_embedding(walks, gene_clusters, options);

%% Use the embedding with SVMs
fprintf('[Perfoming our version]\n');
%run_svm(x, anno, test_filt);

%% Performs the base Mashup for comparison
fprintf('[Performing base version]\n');
if options.walk.use_go_link
    x_base = svd_embed(walks(1:end-1,:,:), options.embedding.ndim);
else
    x_base = svd_embed(walks, options.embedding.ndim);
end
%run_svm(x_base, anno, test_filt);

%% plot

[least_anno, most_anno]=calc_anno_frequency(anno);

anno_one = sum(anno,1)~=0;

least_anno_label = least_anno(anno_one);
most_anno_label = most_anno(anno_one);
[iid,jid,s] = find(label_clusters);
least_anno_cluster = zeros(size(least_anno_label));
most_anno_cluster = zeros(size(most_anno_label));


for i = 1:size(label_clusters,1)
    least_anno_cluster(least_anno_label == iid(i)) = jid(i);
    most_anno_cluster(most_anno_label == iid(i)) = jid(i); 
end



%% plot by clusters
Y = tsne(x_base');
Y = Y(anno_one',:);
%gscatter(Y(:,1),Y(:,2),least_anno(anno_one));
scatter(Y(:,1),Y(:,2),35,least_anno_cluster, 'filled');
figure(2)
scatter(Y(:,1),Y(:,2),35,most_anno_cluster, 'filled');
scatter(Y2(:,1),Y2(:,2))
title('yeast\_mf\_2D', 'FontSize', 20)
savefig('fig1.fig')


xreal = real(x);
Y2=tsne(xreal');
Y2 = Y2(anno_one',:);
figure(3)
scatter(Y2(:,1),Y2(:,2),35,least_anno_cluster, 'filled');
figure(4)
scatter(Y2(:,1),Y2(:,2),35,most_anno_cluster, 'filled');


%% plot by labels

Y = tsne(x_base');
Y = Y(anno_one',:);
%gscatter(Y(:,1),Y(:,2),least_anno(anno_one));
scatter(Y(:,1),Y(:,2),35,least_anno(anno_one), 'filled');
scatter(Y(:,1),Y(:,2),35,most_anno(anno_one), 'filled');
%gscatter(Y(:,1),Y(:,2))
title('yeast\_mf\_2D', 'FontSize', 20)
savefig('fig1.fig')

figure(2);
xreal = real(x);
Y2=tsne(xreal');
Y2 = Y2(anno_one',:);
scatter(Y2(:,1),Y2(:,2),35,least_anno(anno_one), 'filled');
scatter(Y2(:,1),Y2(:,2),35,most_anno(anno_one), 'filled');


% figure(3);
% xabs = abs(x);
% Y3 = tsne(xabs');
% Y3 = Y3(anno_one',:);
% scatter(Y2(:,1),Y2(:,2),35,most_anno(anno_one), 'filled');

%gscatter(Y3(:,1),Y3(:,2))

figure(4)
Y4 = tsne(xreal','NumDimensions',3);
Y4 = Y4(anno_one',:);
scatter3(Y4(:,1),Y4(:,2),Y4(:,3),35,least_anno(anno_one),'filled')
scatter3(Y4(:,1),Y4(:,2),Y4(:,3),35,most_anno(anno_one),'filled')
savefig('fig4.fig')


scatter3(Y4(:,1),Y4(:,2),Y4(:,3))

figure(5)
Y5 = tsne(x_base','NumDimensions',3);
Y5 = Y5(anno_one',:);
colormap default
scatter3(Y5(:,1),Y5(:,2),Y5(:,3),35,least_anno(anno_one),'filled')
scatter3(Y5(:,1),Y5(:,2),Y5(:,3),35,most_anno(anno_one),'filled')
savefig('fig5.fig')
