nnodes = 60;
nclasses = 20;
noise = 0.03;
clusters = gcd(nnodes, nclasses) / 4;
genes_per_label = nnodes / clusters;
labels_per_gene = nclasses / clusters;

base_labels = zeros(nclasses, nnodes);

for i=0:(clusters-1)
    base_labels(i*labels_per_gene + 1 : (i+1)*labels_per_gene, ... 
        i*genes_per_label + 1 : (i+1)*genes_per_label) = 1;
end

perturb_labels = rand(nclasses, nnodes) < noise;

labels = xor(base_labels, perturb_labels);

imagesc(labels)

shuffled_labels = labels(randperm(size(labels,1)),:);


for i=1:1000
   r1 = randi([1,nclasses]);
   r2 = randi([1,nclasses]);
   c1 = randi([1,nnodes]);
   c2 = randi([1,nnodes]);
   temp = shuffled_labels(r1,:);
   shuffled_labels(r1,:) = shuffled_labels(r2,:);
   shuffled_labels(r2,:) = temp;
   temp = shuffled_labels(:,c1);
   shuffled_labels(:,c1) = shuffled_labels(:,c2);
   shuffled_labels(:,c2) = temp;
end

imagesc(shuffled_labels);