function [train_filt, test_filt, ntrain, ntest, ...
        train_labels, test_labels] = create_train_test(anno, options)
    
    addpath code/evaluate;
    
    test_fraction = options.test_fraction;
    
    [~, test_filt] = cv_partition(anno, test_fraction); 
    % filters proteins with no labels
    label_filt = (sum(anno) > 0).'; 
    % removes the unlabeled proteins from the train and test sets
    train_filt =(~test_filt) & label_filt;
    test_filt = test_filt & label_filt;

    ntest = sum(test_filt);
    ntrain = sum(train_filt);

    % applies the filters as masks to the annotations
    train_labels = anno.*(train_filt.');
    test_labels = anno.*(test_filt.');
end

