function [accuracy] = majority_voting(X, anno, test_filt,train_filt, knn, dist_mat)
%----------------------------------------------------
% function prediction with majority voting by kNN
%----------------------------------------------------
    filt = sum(anno) > 0;
    %anno = anno(:,filt);
    %X = X(:,filt);
  
    [nclass, ngene] = size(anno);
      
    test_filt = test_filt(filt);
    train_ind = find(train_filt);
    train_filt = train_filt(filt);
    test_ind = find(test_filt);
    train_filt_ind = find(train_filt);
    %dist_mat = dist_mat(filt,filt);
    ntest = sum(test_filt);
    ntrain = ngene - ntest;
    
    
    
    acc = zeros(ntest,1);
    
    for p = 1:ntest
        i = test_ind(p);    % gene id
        
        % collect all function labels from kNN neighbors
        fun = [];
        fun_weight = [];
        
        iknn = knn(:,i);
        
        voting_knn = intersect(iknn, train_ind);
        
        if isempty(voting_knn)==1 % kNN neighbors not labelled thus can't vote
            acc(p) = 10;          % delete these from the accuracy calculation
        else
            for j = 1:length(voting_knn)
                %j_filt = find(train_ind==voting_knn(j));
                j_fun = find(anno(:, voting_knn(j)));
                j_num_fun = length(j_fun);
                if j_num_fun ~= 0
                    %dist = dist_mat(i,knn(i,voting_knn(j)));
                    dist = dist_mat(i,voting_knn(j));
                    j_weight = 1/dist*ones(1,j_num_fun);
                    fun = [fun;j_fun]; %array of all funcitons from knn, including repeat terms.
                    fun_weight = [fun_weight,j_weight];
                end
            end         
        end
        
        % count frequency of all functions voted by kNN neighbors 
        fun_list=unique(fun);  
        fun_freq=zeros(size(fun_list)); % array of weight/frequency of functions from neighbors
        
        for j = 1:length(fun_list)
            fun_freq(j) = sum(fun_weight(find(fun==fun_list(j))));
        end
        
        [~,max_freq_ind] = max(fun_freq);
        pred_fun = fun_list(max_freq_ind);
        true_fun = find(anno(:,i));
        if ~isempty(intersect(pred_fun, true_fun)) 
            acc(p)=1;
        else
            acc(p)=0;
        end 
    end
mean(acc)

end