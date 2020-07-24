function [w, x, P, fval] = ss_vector_embedding(Q, ndim, maxiter, labels, alpha, beta)
    [nnode ncontext] = size(Q);
    nparam = (nnode + ncontext) * ndim;
    
    % create the must-link and cannot-link matrices (0s and 1s)
    ml = mustlink(labels);
    cl = cannotlink(ml, labels);
    % reweight them using alpha and beta
    n_ml = sum(sum(ml))/2;
    n_cl = sum(sum(cl))/2;
    
    ml = (beta/(n_ml*n_ml)) * ml;
    cl = (alpha/(n_cl*n_cl)) * cl;
    
    % creates the combined constraint matrix
    S = ml - cl;
    
    %% Optimize
    opts = struct('factr', 1e4, 'pgtol', 0, 'm', 5, 'printEvery', 50, 'maxIts', maxiter);

    while true
      %% Initialize vectors
      fprintf('Initializing vectors ... '); tic
      wx = rand(ndim, nnode + ncontext) / 10 - .05;
      fprintf('done. '); toc

      opts.x0 = wx(:);
      [xopt, fval, info] = lbfgsb(@optim_fn, -inf(nparam,1), inf(nparam,1), opts);
      if info.iterations > 10
        break
      end
      fprintf('Premature termination (took %d iter to converge); trying again.\n', info.iterations);
      info
    end
    wx = reshape(xopt, ndim, nnode + ncontext);

    fprintf('Done.\n');
    
    %% Summarize output
    w = wx(:,1:ncontext);
    x = wx(:,ncontext+1:end);
    P = P_fn(w,x);
    fval = obj_fn(P);
    
    function [fval, grad] = optim_fn(wx)
        
        wx = reshape(wx, ndim, nnode + ncontext);
        w = wx(:,1:ncontext);
        x = wx(:,ncontext+1:end);
        
        P = P_fn(w,x);

        % add the supervised penalty
        fval = obj_fn(P,S,x);

        % unsupervised gradietns
        wgrad = x * (P-Q);
        xgrad = w * (P-Q)';
        
        % supervised gradient
        
        grad = [wgrad, xgrad];

        grad = grad(:);
    end

    function P = P_fn(w, x)
        P = exp(x' * w);
        P = bsxfun(@rdivide, P, sum(P));
    end

    function res = obj_fn(P,S,x)
        unsup = zeros(ncontext,1);
        for j = 1:ncontext
            unsup(j) = kldiv(Q(:,j),P(:,j));
        end
        % pairwise distances times the penalty (squared euclidean)
        sup = (squareform(pdist(x')).^2).*S;
        % adds both unsupervised (kl-divergence) and supervised (distance)
        res = sum(unsup) + sum(sum(sup));
    end
   
    function res = kldiv(p,q)
        filt = p > 0;
        res = sum(p(filt) .* log(p(filt) ./ q(filt)));
    end

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
