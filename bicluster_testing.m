addpath code;
addpath mtba-nix/mtba-nix/mtba-nix/mtba;

% use human or yeast data
options.org = 'yeast';

% which type of annotations to use
% options: {bp, mf, cc} for human GO,
%          {level1, level2, level3} for yeast MIPS
options.onttype = 'bp'; 

% consider terms in a specific size range (GO only)
% examples: [11 30], [31 100], [101 300]
options.ontsize = [31 100];       
                    
options.num_clusters = 12;
%% Load gene list
fprintf('[Loading annotations]\n');
[genes, ngene, anno] = load_anno(options);
fprintf('Number of functional labels: %d\n', size(anno, 1));

%train_filt = sum(anno, 1) > 0;

%fprintf('[Performing Biclustering]')

%num_clusters = options.num_clusters;

%clusters_data = spectralCoClustering(anno(:,train_filt), num_clusters, 1);

%imagesc(anno(:,train_filt));