function x = svd_cluster_embed(walks, gene_clusters, options)

    embedding = options.embedding;
    ndim = embedding.ndim;
    cl_penalty = embedding.cannotlink_penalty;
    ml_penalty = embedding.mustlink_penalty;
    num_clusters = options.num_clusters;
    

    [nnetworks, ngene, ~] = size(walks);
    
    RR_sum = zeros(ngene);
    % add the networks together without dummy nodes
    for i = 1:nnetworks
        W = squeeze(walks(i,:,:));
        R = log(W + 1/ngene); % smoothing
        RR_sum = RR_sum + R * R';
    end
    clear R
    
    % normalize to have mean 0
    RR_Sum = RR_sum - mean(RR_sum(:));
    % normalize to have std 1
    RR_Sum = RR_Sum / std(RR_sum(:));
    
    
    nnodes = ngene + num_clusters;
    augmented = zeros(nnodes);
    augmented(1:ngene, 1:ngene) = RR_Sum;
    % adds the cannot link between dummy nodes
    augmented(ngene+1:nnodes,ngene+1:nnodes) = -cl_penalty + (cl_penalty * eye(num_clusters));
    % adds the must links from dummy nodes to genes
    augmented(ngene+1:nnodes,1:ngene) = ml_penalty * gene_clusters;
    augmented(1:ngene,ngene+1:nnodes) = ml_penalty * gene_clusters';
        
    [V, d] = eigs(augmented, ndim);
    x = diag(sqrt(sqrt(diag(d)))) * V';
    x = x(1:ngene, 1:ngene);
end

