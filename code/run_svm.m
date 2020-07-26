function [acc, f1, aupr] = run_svm(x, anno, test_filt)
    % Parameters
    gvec = -3:1:0;
    cvec = -2:1:2;
  
    % Scale features
    maxval = max(x, [], 2);
    minval = min(x, [], 2);
    x = bsxfun(@times, bsxfun(@minus, x, minval), 1 ./ (maxval - minval));
  
    % Filter genes with no annotations
    filt = sum(anno) > 0;
    anno = anno(:,filt);
    x = x(:,filt);
  
    [nclass, ngene] = size(anno);
      
    test_filt = test_filt(filt);
    ntest = sum(test_filt);
    ntrain = ngene - ntest;
  
    fprintf('Pregenerating kernels:\n');
    rbfK = cell(length(gvec), 1);
    for i = 1:length(gvec)
        fprintf('%d / %d ... ', i, length(gvec)); tic
        rbfK{i} = rbf_kernel(x', x', 10^gvec(i));
        fprintf('done. '); toc
    end
  
      
    retmax = -inf;
    gmax = 1;
    cmax = 1;
    fprintf('Running nested cross validation...\n');
    for gi = 1:length(gvec)
        for ci = 1:length(cvec)
            tt = tic;

            Ktrain = rbfK{gmax}(~test_filt,~test_filt);
            Ktest = rbfK{gmax}(test_filt,~test_filt);
  
            class_score = zeros(ntest, nclass);
            %parfor s = 1:nclass
            for s = 1:nclass
                Ytrain = full(double(anno(s,~test_filt)') * 2 - 1);
                Ytest = full(double(anno(s,test_filt)') * 2 - 1);
  
                model = svmtrain(Ytrain, [(1:size(Ktrain,1))', Ktrain], ['-t 4 -b 1 -q -c ', num2str(10^cvec(ci))]);
                posind = find(model.Label > 0);
                if ~isempty(posind)
                    [~, ~, dec] = svmpredict(Ytest, [(1:size(Ktest,1))', Ktest], model, '-q');
                    class_score(:,s) = dec(:,posind);
                end
            end
            % TODO make it just the passed indices
            [~, ~, cv_result] = evaluate_performance(class_score, anno(:,test_filt)');

            if retmax < cv_result
                retmax = cv_result;
                gmax = gi;
                cmax = ci;
            end
  
            fprintf('gi:%d, ci:%d, cv_result:%f, ', gi, ci, cv_result); toc(tt)
        end
    end
  

    % use the best parameters we computed above
    fprintf('Using best paramaters\n')
    Ktrain = rbfK{gmax}(~test_filt,~test_filt);
    Ktest = rbfK{gmax}(test_filt,~test_filt);
  
    class_score = zeros(ntest, nclass);
    %parfor s = 1:nclass
    for s = 1:nclass
        Ytrain = full(double(anno(s,~test_filt)') * 2 - 1);
        Ytest = full(double(anno(s,test_filt)') * 2 - 1);
  
        model = svmtrain(Ytrain, [(1:ntrain)', Ktrain], ['-t 4 -b 1 -q -c ', num2str(10^cvec(cmax))]);
        posind = find(model.Label > 0);
        if ~isempty(posind)
            [~, ~, dec] = svmpredict(Ytest, [(1:ntest)', Ktest], model, '-q');
            class_score(:,s) = dec(:,posind);
        end
    end
  

    [acc, f1, aupr] = evaluate_performance(class_score, anno(:,test_filt)');
    fprintf('[Trial #%d] acc: %f, f1: %f, aupr: %f\n', 1, acc, f1, aupr);
end
  