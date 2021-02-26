# FLAME

Functional Label Augmented Multinetwork Embedding

This code builds upon the [Mashup](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5225290/) algorithm for embedding networks. The primary advantage in our approach is the integration of the functional annotations. This allows us to create an embedding that incorportates both network topology and information about what proteins are doing in the cell. This leads to performance gains when trying to infer the functions of other proteins in the networks compared to only using the topology.

## Setup 

To run the basic version of the code which uses a singular value decomposition to embed the networks, and a k-nearest neighbor classifier to infer protein functions, there are no additional dependencies.

Optionally, there are two dependencies that can be used to run Mashup as a baseline for comparison. These are
 - [libsvm](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) which is used to train an SVM classifier.
 - [L-BFGS-B-C.](https://github.com/stephenbeckr/L-BFGS-B-C) which is used to perform an alternative embedding.


## Data

In the future we will release the dataset that we have been working on but in the meantime, you can use the supplementary Mashup data which can be accessed [here.](http://cb.csail.mit.edu/cb/mashup/)
