## UPGMA clustering

This is an efficient implementation of a hierarchical clustering method,
UPGMA. UPGMA is used (most commonly) to create phylogenetic trees for
use in biological research.

This implementation allows creation of a tree with any arbitrary distance
function and member objects. While Python is generally not ideal for 
a highly optimized data structure, the flexibility allowed by Python 
and the gains to time complexity motivated me to implement UPGMA
clustering myself.

This method has superior time complexity to other implementations (see
*Gronau I, Moran S (2007) Optimal implementations of UPGMA and other
common clustering algorithms*, from which this implementation was
adapted). Note that one proposed optimization - storing distances in a
heap - was not implemented.
