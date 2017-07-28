void tabulate
(const int *x, const int *n, const int *nbin, int *ans) {
    int i, tmp;
    for (i=0; i < *nbin; i++) ans[i]=0L; 
    for (i=0; i < *n; i++) {
        tmp = x[i];
        if( (tmp>0) & (tmp<(*nbin+1L)) )   
        ans[tmp-1L] ++;
    }
}

void phangorn_reorder
(const int *from, const int *to, const int *n, 
 const int *sumNode, int *neworder, int *root) { 
    int i, j, sum=0, k, Nnode, ind, *ord, *csum, *tips, *stack, z=0;  // l, 
    double *parent;
    int m=sumNode[0];
    parent = (double *) R_alloc((*n), sizeof(double));
    tips = (int *) R_alloc(m, sizeof(int));
    ord = (int *) R_alloc((*n), sizeof(int));
    csum = (int *) R_alloc( (m+1), sizeof(int));
    stack = (int *) R_alloc(m, sizeof(int));
    for(j=0;j<(*n);j++) parent[j] = (double)from[j];
   
    for(j=0;j<(*n);j++) ord[j] = j;
    for(j=0;j<m;j++) tips[j] = 0;
        
    rsort_with_index(parent, ord, *n);
    tabulate(from, n, sumNode, tips);
    csum[0]=0;
    for(i=0;i<(*sumNode);i++){
        sum+=tips[i];                 
        csum[i+1] = sum;                        
    }      
    k = (*n)-1;
    Nnode = 0;
    stack[0] = *root;
    
    while(z > -1){
        j=stack[z];          
        if(tips[j]>0){   
            for(i=csum[j];i<csum[j+1];i++){
                ind = ord[i];                     
                neworder[k] = ind + 1;        
                stack[z] = to[ind]-1;
                k -=1;
                z++;                            
            }         
            Nnode += 1; 
            }
        z--;       
    }                
    root[0]=Nnode;     
}


/* reorder_phylo.c       2012-09-03 */
/* Copyright 2008-2012 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

static int iii;

void ape_foo_reorder(int node, int n_tips, int n_nodes, int *parent, int *child, int *neworder, int *L, int *pos)
{
	int i = node - n_tips - 1, j, k;

/* 'i' is the C index corresponding to 'node' */

	for (j = 0; j < pos[i]; j++) {
		k = L[i + n_nodes * j];
		neworder[iii++] = k + 1;
		if (child[k] > n_tips) /* is it an internal edge? */
			ape_foo_reorder(child[k], n_tips, n_nodes, parent, child, neworder, L, pos);
	}
}

void ape_bar_reorder(int node, int n_tips, int n_nodes, int *parent, int *child, int *neworder, int *L, int *pos)
{
	int i = node - n_tips - 1, j, k;

	for (j = pos[i] - 1; j >= 0; j--)
		neworder[iii--] = L[i + n_nodes * j] + 1;

	for (j = 0; j < pos[i]; j++) {
		k = child[L[i + n_nodes * j]];
		if (k > n_tips)
			ape_bar_reorder(k, n_tips, n_nodes, parent, child, neworder, L, pos);
	}
}

static R_NativePrimitiveArgType ape_neworder_phylo_t[] = {
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
void ape_neworder_phylo(int *n_tips, int *parent, int *child, int *n_edges, int *neworder, int *order)
/* n_tips: nb of tips
   n_nodes: nb of nodes
   n_edges: nb of edges */
{
	int i, j, k, *L, *pos, n_nodes = *n_edges - *n_tips + 1, max_node = *n_tips - n_nodes + 1;

/* max_node is the largest value that a node degree can be */

/* L is a one-dimensional array storing, for each node, the C indices of the rows of
   the edge matrix where the node is ancestor (i.e., present in the 1st
   column). It is used in the same way as a matrix (which is actually
   a vector) is used in R as a 2-d structure. */

	L = (int*)R_alloc(n_nodes * max_node, sizeof(int));

/* pos gives the position for each 'row' of L, that is the number of elements
   which have already been stored for that 'row'. */

	pos = (int*)R_alloc(n_nodes, sizeof(int));
	memset(pos, 0, n_nodes * sizeof(int));

/* we now go down along the edge matrix */

	for (i = 0; i < *n_edges; i++) {
		k = parent[i] - *n_tips - 1; /* k is the 'row' index in L corresponding to node parent[i] */
		j = pos[k]; /* the current 'column' position corresponding to k */
		pos[k]++; /* increment in case the same node is found in another row of the edge matrix */
		L[k + n_nodes * j] = i;
	}

/* L is now ready: we can start the recursive calls. */
/* We start with the root 'n_tips + 1': its index will be changed into
   the corresponding C index inside the recursive function. */

	switch(*order) {
	case 1 : iii = 0;
		ape_foo_reorder(*n_tips + 1, *n_tips, n_nodes, parent, child, neworder, L, pos);
		break;
	case 2 : iii = *n_edges - 1;
		ape_bar_reorder(*n_tips + 1, *n_tips, n_nodes, parent, child, neworder, L, pos);
		break;
	}
}

#define DO_NODE_PRUNING\
    /* go back down in `edge' to set `neworder' */\
    for (j = 0; j <= i; j++) {\
        /* if find the edge where `node' is */\
        /* the descendant, make as ready */\
        if (edge2[j] == node) ready[j] = 1;\
	if (edge1[j] != node) continue;\
	neworder[nextI] = j + 1;\
	ready[j] = 0; /* mark the edge as done */\
	nextI++;\
    }

static R_NativePrimitiveArgType ape_neworder_pruningwise_t[] = {
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
void ape_neworder_pruningwise(int *ntip, int *nnode, int *edge1,
			  int *edge2, int *nedge, int *neworder)
{
    int *ready, *Ndegr, i, j, node, nextI, n_tips;

    nextI = *ntip +  *nnode;
    Ndegr = (int*)R_alloc(nextI, sizeof(int));
    memset(Ndegr, 0, nextI*sizeof(int));
    for (i = 0; i < *nedge; i++) (Ndegr[edge1[i] - 1])++;

    ready = (int*)R_alloc(*nedge, sizeof(int));

    /* `ready' indicates whether an edge is ready to be */
    /* collected; only the terminal edges are initially ready */
    for (i = 0; i < *nedge; i++)
	    ready[i] = (edge2[i] <= *ntip) ? 1 : 0;

    /* `n_tips' counts the number of times a node has been seen. */
    /* This algo will work if the tree is in cladewise order, */
    /* so that the nodes of "cherries" will be contiguous in `edge'. */
    n_tips = 0;
    nextI = 0;
    while (nextI < *nedge - Ndegr[*ntip]) {
        for (i = 0; i < *nedge; i++) {
            if (!ready[i]) continue;
	    if (!n_tips) {
	        /* if found an edge ready, initialize `node' and start counting */
	        node = edge1[i];
		n_tips = 1;
	    } else { /* else counting has already started */
	        if (edge1[i] == node) n_tips++;
		else {
		    /* if the node has changed we checked that all edges */
		    /* from `node' have been found */
		    if (n_tips == Ndegr[node - 1]) {
		        DO_NODE_PRUNING
		    }
		    /* in all cases reset `n_tips' and `node' and carry on */
		    node = edge1[i];
		    n_tips = 1;
		}
	    } /* go to the next edge */
	    /* if at the end of `edge', check that we can't do a node */
	    if (n_tips == Ndegr[node - 1]) {
	        DO_NODE_PRUNING
		n_tips = 0;
	    }
        }
    }
    for (i = 0; i < *nedge; i++) {
        if (!ready[i]) continue;
        neworder[nextI] = i + 1;
        nextI++;
    }
}