function walks = compute_rwr(network_files, ngene, train_filt, options)
    addpath code/random_walk;
    addpath code/constraints;
    
    restart_prob = options.walk.restart_prob;
    use_go_link = options.walk.use_go_link;
    org = options.org;
    onttype = options.onttype;

    walks = unsupervised_rwr(network_files, ngene, restart_prob);
    
    if use_go_link
        fprintf('[Adding Constraint Matrix]')
        filepath = sprintf('data/annotations/%s/%s-total-index.txt', org, onttype);
        sparse_link = go_link_matrix(filepath, ngene, train_filt);
        
        dense_link = full(sparse_link);
        link_markov = markov_mat(dense_link);
        
        link_walk = rwr(link_markov, restart_prob);
        link_walk = reshape(link_walk, 1,ngene,ngene); % add an extra index for appending
        walks = cat(1, walks, link_walk); % tacks the extra walk onto the end
    end 
end

