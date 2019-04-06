HOW TO USE:
In order to generate the distance matrix, edges and newick files
execute "python3 hw3.py hw3.fna"
edges.txt is a tab delimited file that has 3 columns. The first column is the parent node, second column is the descendent node and the third column is the edge length. Is is shown in preorder.
tree.txt is in NEWICK format with all edge distances with only tips named.
HOW TO BUILD TREE
execute "Rscript hw3-plot-newick.r tree.txt hw3-tip-labels.txt" to build tree graph based on NEWICK form
excutre "Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt" to build tree graph based on edges form