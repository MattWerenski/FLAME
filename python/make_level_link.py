# organism being used
org = 'human'

# go tree being used
ont = 'mf'

# path to the list of string genes
string_gene_filepath = f'data/networks/{org}/{org}_string_genes.txt'

# path to the list of go genes
go_gene_filepath     = f'data/annotations/{org}/go_{org}_ref_genes.txt'
# path to the list of go terms
go_term_filepath     = f'data/annotations/{org}/go_{org}_ref_{ont}_terms.txt'
# path to a list of genes-go term pairs
go_adj_filepath      = f'data/annotations/{org}/go_{org}_ref_{ont}_adjacency.txt'

# path to a file that maps go terms to an index
go_term_map_filepath = f'data/annotations/{org}/graph/go_{ont}.map'
# path to the a file containing the structure of GO
go_dag_filepath      = f'data/annotations/{org}/graph/go_{ont}.links'

# lowest level to include in the data sets (h.bp - 5, h.mf - 8, y.mf - 5)
min_level            = 8

# where to output the dataset to
location             = f'data/annotations/{org}'

# whether to create files for each level's gene pairs
create_individual    = False
# whether to create one big file containing all the levels
create_total         = True

# ================= Parsing gene files ====================

# reads in the list of string genes
file_gene_list = open(string_gene_filepath, 'r')
string_genes = []
for gene_line in file_gene_list.readlines():
    gene = gene_line.split()[0]
    string_genes.append(gene)
file_gene_list.close()

# reads the list of GO genes in
file_go_gene_list = open(go_gene_filepath, 'r')
go_genes = []
for gene_line in file_go_gene_list.readlines():
    gene = gene_line.split()[0]
    go_genes.append(gene)
file_go_gene_list.close()


# Creates a map from genes to their string index which is useful
# for converting from GO indices to string indices
# (shifted +1 becasuse matlab is 1-indexed)
idx = 1
string_map = {}
for gene in string_genes:
    string_map[gene] = idx
    idx += 1

# creates a map from GO genes back to their GO indexing
# (This is a bit of a nightmare since the GO and string indexing aren't shared)
idx = 1
gene_map = {}
for gene in go_genes:
    gene_map[gene] = idx
    idx += 1

# Creates a set for checking if a string gene is also a GO gene
filt = set()
for i in range(len(string_genes)):
    gene = string_genes[i]
    if gene in go_genes:
        filt.add(gene)

# ================= Parsing GO term files ====================

# reads in the list of GO terms being used
file_go_term_list = open(go_term_filepath, 'r')
go_terms = []
for term_line in file_go_term_list.readlines():
    term = term_line.split()[0]
    go_terms.append(term)
file_go_term_list.close()


# creates dict with the gene indices as keys and the 
# list of their go term indices as values
file_go_adj = open(go_adj_filepath, 'r')
gene_anno = {}
# pre-paring the dict, we need +1 to cope with matlab's 1-indexing
for i in range(len(go_genes)):
    gene_anno[i+1] = {'orig': []}

for map_line in file_go_adj.readlines():
    [gene_index, term_index] = map_line.split()
    gene_index = int(gene_index)
    term_index = int(term_index)

    # did you know that list + list does concatenation?
    gene_anno[gene_index]['orig'] += [term_index]
file_go_adj.close()


# reads in the mapping from go terms to their indices
# basically instead of working on terms directly, they index them to ints
# (if this is confusing it's because it is.)
file_go_term_map = open(go_term_map_filepath, 'r')

go_term_map = [None] # pads for 1-indexing
for map_line in file_go_term_map.readlines():
    [term, index] = map_line.split()
    index = int(index)
    
    # the index is sequential and this means we've seen it before
    if len(go_term_map) == index + 1:
        go_term_map[index] += [term]
    # this means its the first occurrence
    else:
        go_term_map += [[term]] # 8/10 slick
file_go_term_map.close()

# ================= Creating GO term DAG ====================

# creates a DAG of the GO terms 
file_go_links = open(go_dag_filepath, 'r')
# index_tree uses the go term indices (the same as in go_term_map)
# and fills in a dict the indices as values and the tree structure as values
index_tree = {}
for link_line in file_go_links.readlines():
    [parent_index, child_index] = link_line.split()
    parent_index = int(parent_index)
    child_index = int(child_index)

    if not parent_index in index_tree:
        index_tree[parent_index] = {'parents':[], 'children':[]}
    if not child_index in index_tree:
        index_tree[child_index] = {'parents':[], 'children':[]}
    index_tree[parent_index]['children'].append(child_index)
    index_tree[child_index]['parents'].append(parent_index)
file_go_links.close()

# ================= Computing term levels ====================

# find all the terms which have no parents and are therefore the root terms
# (this should probably be only one term but the list is to be safe)
roots = []
for index in index_tree:
    if len(index_tree[index]['parents']) == 0:
        roots.append(index)


# computing term depths in the DAG
#preparing all nodes with auxillary fields
for index in index_tree:
    index_node = index_tree[index]
    index_node['visited'] = False
    index_node['level'] = 0 # indicates unknown

#preparing root nodes 
queue = roots + []
for index in queue:
    index_node = index_tree[index]
    index_node['level'] = 1
    index_node['visited'] = True

# propagates the levels down the dag using BFS 
while len(queue) > 0:
    parent_index = queue.pop(0)
    parent_node = index_tree[parent_index]
    
    child_indices = parent_node['children']
    for child_index in child_indices:
        child_node = index_tree[child_index]
        if child_node['visited']:
            continue
        else:
            child_node['visited'] = True
            child_node['level'] = parent_node['level'] + 1
            queue.append(child_index)

max_level = 0
for index in index_tree:
    level = index_tree[index]['level']
    max_level = max(level, max_level)

level_counts = [0] * max_level # this is also a thing
for index in index_tree:
    level = index_tree[index]['level']
    level_counts[level - 1] += 1

# ================= Propagating parent terms down the DAG ====================

#prepares each tree node
for index in index_tree:
    index_node = index_tree[index]
    index_node['ancestors'] = {index}

for i in range(max_level):
    # go level by level pushing down ancestors
    cur_level = i + 1
    
    for index in index_tree:
        index_node = index_tree[index]
        
        # only look at nodes on the right level
        if index_node['level'] != cur_level:
            continue
        
        parents = index_node['parents']
        # grabs all the parents' ancestors (list of sets)
        ancestors = list(map(lambda pidx: index_tree[pidx]['ancestors'], parents))
        # unions over all parent ancestors + self (* unpacks the list of sets)
        index_node['ancestors'] = index_node['ancestors'].union(*ancestors) # 9/10 slick

# associates each gene index with all its ancestors
for gene_index in gene_anno:
    gene_node = gene_anno[gene_index]
    # 10/10 slick
    gene_node['comp'] = set().union(*list(map(lambda idx: index_tree[idx]['ancestors'], gene_node['orig'])))

# ================= Associating genes to terms ====================

# prepares each term index
for index in index_tree:
    index_node = index_tree[index]
    index_node['gene_inds'] = []
    
# go over each gene and append it to the indices
for gene_index in gene_anno:
    term_indices = list(gene_anno[gene_index]['comp'])
    list(map(lambda tidx: index_tree[tidx]['gene_inds'].append(gene_index), term_indices))

# ================= Computing best matchings ====================

# used to track the matchings we've seen on any level
total_matchings = set()

# tracks the matchings at every level (but only the highest level)
matchings = []
# can't do multiplication of lists because each entry references the same set
for i in range(max_level):
    matchings.append(set())

# go from max level -> min level in decreasing order
for working_level in range(max_level, min_level-1, -1):
    for index in index_tree:
        index_node = index_tree[index]
        
        if index_node['level'] != working_level:
            continue
        
        # list of genes using the current term 
        # which means they share a term at the working level
        gene_inds = index_node['gene_inds']

        for gene_index1 in gene_inds:
            for gene_index2 in gene_inds:
                # don't count the self-edges
                if gene_index1 == gene_index2:
                    continue
                    
                key = f'{gene_index1}-{gene_index2}'

                # used so we don't overwrite a high level with a low level
                if not key in total_matchings:
                    matchings[working_level - 1].add(key)
                    total_matchings.add(key)

# ================= Creating files ====================

# creates lists of genes for each separate level of GO
if create_individual:
    for level in range(max_level, min_level-1, -1):
        file_index = open(f'{location}/{ont}-level-{level}-index.txt','w')
        file_gene = open(f'{location}/{ont}-level-{level}-gene.txt','w')

        for pair in matchings[level - 1]:
            [index1, index2] = pair.split('-')
            index1 = int(index1)
            index2 = int(index2)
            
            gene1 = go_genes[index1]
            gene2 = go_genes[index2]
            
            if gene1 in filt and gene2 in filt:
                adj1 = string_map[gene1]
                adj2 = string_map[gene2]
                file_index.write(f'{adj1}\t{adj2}\n')
                file_gene.write(f'{gene1}\t{gene2}\n')
            
        file_index.close()
        file_gene.close()

# compiles all the genes that are linked into a single file
for i in range(max_level):
    print(f'level: {i+1} - matches: {len(matchings[i])}')
print('--------')

if create_total:
    file_index = open(f'{location}/{ont}-total-index.txt','w')
    file_gene = open(f'{location}/{ont}-total-gene.txt','w')

    for level in range(min_level - 1, max_level):
        matching = matchings[level]
        print(f'level: {level+1} - matches: {len(matching)}')
        for pair in matching:
            [index1, index2] = pair.split('-')
            index1 = int(index1)
            index2 = int(index2)

            gene1 = go_genes[index1 - 1]
            gene2 = go_genes[index2 - 1]

            if gene1 in filt and gene2 in filt:
                adj1 = string_map[gene1]
                adj2 = string_map[gene2]
                file_index.write(f'{adj1}\t{adj2}\t{level+1}\n')
                file_gene.write(f'{gene1}\t{gene2}\t{level+1}\n')
            
    file_index.close()
    file_gene.close()