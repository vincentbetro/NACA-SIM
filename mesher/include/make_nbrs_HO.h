#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Linked_List.h"
#include "List.h"

//this routine takes a list of triangles to make nabors from, number of nodes, number of triangles, expanded tri array (even though we only use 0, 1, 2), nabor array (passed in init to -1 and out full), and node_node hash table

void make_nbrs_HO(Linked_List *mknbr, int nn, int nt, int **tri, int **nbr, Linked_List **Nhash);
