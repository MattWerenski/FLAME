function x = svd_supervised_embed(walks, ndim, labels, ... 
    cl_penalty, ml_penalty)

    [nnetworks, ngenes, ~] = size(walks);
    Q = zeros(ngenes*nnetworks, ngenes);
    
    % probably a better way to do this
    for i=1:nnetworks
        Q((i-1)*ngenes+1 : i*ngenes, :) = squeeze(walks(i,:,:)) / nnetworks;
    end
    
    % create the must-link and cannot-link matrices (0s and 1s)
    ml = mustlink(labels);
    cl = cannotlink(ml, labels);
    % reweight them using 
    n_ml = sum(sum(ml))/2;
    n_cl = sum(sum(cl))/2;
    ml = (ml_penalty/(n_ml*n_ml)) * ml;
    cl = (cl_penalty/(n_cl*n_cl)) * cl;
    
    % creates the combined constraint matrix
    S = ml - cl;
    nz = nnz(S);
    % Checking how dense supervised matrix is
    fprintf('  Supervised penalty has %d non-zero (%f percent)\n', nz, 100 * nz / (ngenes * ngenes));
    LS = diag(sum(S))-S;
    
    [nnetworks, ngene, ~]  = size(walks);
    RR_sum = zeros(ngene);
    for i = 1:nnetworks
        W = walks(i);
        R = log(W + 1/ngene); % smoothing
        RR_sum = RR_sum + R * R';
    end
    clear R Q
    [V, d] = eigs(RR_sum-LS, ndim);
    %x = diag(sqrt(sqrt(diag(d)))) * V';
    x = V';
     

    function cl = cannotlink(ml, labels)
        [n,~] = size(ml);
        cl = ones(n,n);
        % indicator vector for if a vertex has a label
        unlabelled = sum(labels) == 0;
        % remove the ml constraints
        cl = cl - ml;
        % remove the unlabeled constraints
        cl(unlabelled, :) = 0;
        % remove diagonal
        cl = cl - diag(diag(cl));
    end

    function ml = mustlink(labels)
        [c, n] = size(labels);
        % for each class create a clique adj. mat.
        ml = zeros(n,n);
        for i = 1:c
            % get the indicator vector for that class
            indicator = labels(i,:) == 1;
            % fills in the ml matrix for that class
            ml(indicator, indicator) = 1;
        end
        % remove diagonal
        ml = ml - diag(diag(ml));
    end
end
