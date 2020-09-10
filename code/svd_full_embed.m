function x = svd_full_embed(walks, ndim, labels, ... 
    cl_penalty, ml_penalty)

    [nnetworks, ngene, ~] = size(walks);
    
    % create the must-link and cannot-link matrices (0s and 1s)
    ml = mustlink(labels);
    cl = cannotlink(ml, labels);
    % normally we don't need this check but for debugging I have it
    ml = ml_penalty * ml;
    cl = cl_penalty * cl;
    
    % creates the combined constraint matrix
    constraints = ml - cl;
    nz = nnz(constraints);
    % Checking how dense supervised matrix is
    fprintf('  Supervised penalty has %d non-zero (%f percent)\n', nz, 100 * nz / (ngene * ngene));
    L = diag(sum(constraints))-constraints;
    
    % reshapes the walks so that they are stacked into a very tall matrix
    % the first ngene rows correspond to the first network, rows ngene+1 to 2*ngene
    % are for the second network and so on
    Wfull = reshape(permute(W, [3, 1, 2]), nnetworks * ngene, ngene);

    % stacks nnetworks copies of the constraints. The output matrix has
    % nnetworks * ngene rows and ngene columns.
    Lfull = repmat(L,nnetworks,1);

    % One question is whether or not to keep the values between 0 and 1 in 
    % Wfull + LSfull because that provides a rw-interpretation
    [~, S, V] = svd(Wfull + Lfull);

    x = S(1:ndim, :) * V';
end
