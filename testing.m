pos = 0; neg = 0; pass = 0;
for i = 1:length(neg1)
    id1 = neg1(i); id2 = neg2(i);
    clust1 = gene_clusters(:,id1);
    clust2 = gene_clusters(:,id2);
    if sum(clust1) == 0 || sum(clust2) == 0
        pass = pass + 1;
    end
    if dot(clust1, clust2) == 1
        pos = pos + 1;
    else
        neg = neg + 1;
    end
end