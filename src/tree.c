#include <glib.h>
#include <mpi.h>
#include <string.h>
#include "commonblock.h"

Node * TREEROOT = NULL;
int tree_node_contains_fckey(Node * node, fckey_t * key) {
    fckey_t tmp;
    fckey_xor(&tmp, &node->key, key);
    fckey_clear(&tmp, 3 * node->order);
    return fckey_is_zero(&tmp);
}

int tree_node_contains_node(Node * node, Node * needle) {
    /* returns 1 if needle and node are identical
     * or needle is a child of node */
    return node->order >= needle->order 
       && tree_node_contains_fckey(node, & needle->key);
}
Node * tree_node_find_image(Node * node, Node * needle) {
    /* returns the first parent of node that is identical to needle */
    Node * image = node;
    while(image) {
        if(fckey_cmp(&image->key, &needle->key) == 0 && 
           image->order == needle->order) {
            return image;
        }
        image = (Node *) image->parent;
    }
    return NULL;
}

Node * tree_detach_subtree(Node * node) {
    if(node->parent) {
        int i;
        for(i = 0; i < node->parent->nchildren; i++) {
            if(node->parent->child[i] == node) break;
        }
        for(i++;i < node->parent->nchildren; i++) {
            node->parent->child[i-1] = node->parent->child[i];
            node->parent->child[i] = NULL;
        }
        node->parent->nchildren --;
    }
    return node;
}

void tree_link(Node * root, Node ** firstleaf, InnerNode ** firstinner) {
    TreeIter * iter = tree_iter_new(root);
    InnerNode * myfirstinner = NULL, * inner = NULL;
    Node * myfirstleaf = NULL, * leaf = NULL;
    Node * cur = tree_iter_next(iter);
    while(cur) {
        if(cur->type == NODE_TYPE_INNER) {
            if(myfirstinner == NULL) {
                myfirstinner = (InnerNode*)cur;
            }
            if(inner != NULL) {
                inner->link = cur;
            }
            inner = (InnerNode*)cur;
        } else {
            if(myfirstleaf == NULL) {
                myfirstleaf = cur;
            }
            if(leaf != NULL) {
                leaf->link = cur;
            }
            leaf = cur;
        }
        cur = tree_iter_next(iter);
    }
    tree_iter_free(iter);
    if(firstinner) *firstinner = myfirstinner;
    if(firstleaf) *firstleaf = myfirstleaf;
}

void tree_destroy(Node * root) {
    /* free the entire thing */
    InnerNode * firstinner = NULL;
    Node * firstleaf = NULL;
    tree_link(root, &firstleaf, &firstinner);
    g_slice_free_chain(InnerNode, firstinner, link);
    g_slice_free_chain(Node, firstleaf, link);
}

static InnerNode * node_to_inner(Node * node) {
    InnerNode * rt = (InnerNode *) node;
    if(node->type == NODE_TYPE_LEAF) {
        rt = g_slice_new(InnerNode);
        memcpy(rt, node, sizeof(Node));
        rt->nchildren = 0;
        rt->type = NODE_TYPE_INNER;
        /* now replace the reference to the node in the parent */
        InnerNode * parent = node->parent;
        if(parent) {
            for(int i = 0; i < parent->nchildren; i++) {
                if(parent->child[i] == node) {
                    parent->child[i] = (Node *)rt;
                }
            }
        }
        g_slice_free(Node, node);
    } else {
        g_error("already inner");
    }
    return rt;
}

static Node * create_leaf(Node * parent, par_t * first) {
    Node * rt = g_slice_new0(Node);
    if(parent->type == NODE_TYPE_LEAF) {
        parent = (Node*) node_to_inner(parent);
    }
    rt->type = NODE_TYPE_LEAF;
    rt->first = first;
    rt->npar = 0;
    rt->order = parent->order - 1;
    rt->key = first->fckey;
    fckey_clear(&rt->key, 3 * rt->order);
    rt->parent = (InnerNode*)parent;
    rt->nchildren= 0;
    if(rt->parent->nchildren >= 8) {
        g_error("too many children");
    }
    rt->parent->child[rt->parent->nchildren] = rt;
    rt->parent->nchildren ++;
    return rt;
}

static void tree_build_subtree(Node * subroot, par_t * first, intptr_t npar) {
    par_t * i = first;
    intptr_t SKIPPED = 0;
    Node * node = (Node*) subroot;
    /* now scan over PAR, creating nodes */
    while(i < first + npar) {
        /* the current particle is no longer in current node,
         * time to close parent nodes */
        while(node != (Node*) subroot->parent && 
            ! tree_node_contains_fckey(node, &i->fckey)) {
            node->npar = i - node->first;
            /* here we shall try to merge the children, if 
             * doing so is intended */
            node = (Node *) node->parent;
        }
        if(node == (Node*) subroot->parent) {
            /* particle not in any nodes, need to skip it */
            SKIPPED ++;
            i ++;
            continue;
        }
        /* ASSERTION: 
         * the current particle i is 
         * 1) in the current node, 
         * 2) not in any of the current children of node 
         *    [guarrenteed by morton key sorting]
         * */
        if(node->type == NODE_TYPE_INNER) {
            /* create a child of this inner node
             * particle will be added later */
            node = create_leaf(node, i);
        } else /* must be a leaf node */
        if(node->npar > CB.NodeSplitThresh && node->order > 0) {
            /* split the node */
            i = node->first;
            node = create_leaf(node, i);
        } else {
            /* add particle to the leaf */
            /* code will fall through */
        }
        int extrastep = 0;
        /* TODO: figure a bigger step size, to directly fill
         * the leaf  */
        node->npar += extrastep + 1;
        i += (1 + extrastep);
    }
    while(node != (Node *) subroot->parent) {
        node->npar = i - node->first;
        node = (Node *) node->parent;
    }
    MPI_Allreduce(MPI_IN_PLACE, &SKIPPED, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    ROOTONLY {
        if(SKIPPED > 0) {
            g_warning("particles out of box: %ld\n", SKIPPED);
        }
    }
}

void tree_build() {
    /* first make the root */
    TREEROOT = (Node *) g_slice_new0(InnerNode);
    TREEROOT->order = FCKEY_BITS;
    TREEROOT->type = NODE_TYPE_INNER;
    TREEROOT->first = &PAR(0);
    tree_build_subtree(TREEROOT, &PAR(0), NPAR);
}

Node * tree_locate_fckey(Node * root, fckey_t * key) {
    TreeIter * iter = tree_iter_new(root);
    Node * node = tree_iter_next(iter);
    while(node) {
        if(tree_node_contains_fckey(node, key)) {
            if(node->type == NODE_TYPE_LEAF) {
                break;
            } else {
                node = tree_iter_next(iter);
            }
        }  else {
            node = tree_iter_next_sibling(iter);
        }
    }
    tree_iter_free(iter);
    return node;
}
TreeIter * tree_iter_new(Node * root) {
    TreeIter * rt = g_slice_new0(TreeIter);
    rt->root = root;
    rt->top = -2;
    return rt;
}
void tree_iter_free(TreeIter * iter) {
    g_slice_free(TreeIter, iter);
}
static Node * tree_iter_real_next(TreeIter * iter, int skip_children) {
    #define TOP ((InnerNode*) (iter->stack[iter->top]))
    #define CTOP iter->c[iter->top]
    if(iter->top == -2) {
        iter->top = -1;
        return iter->root;
    }
    if(iter->top == -1) {
        if (skip_children) return NULL;
        iter->stack[0] = iter->root;
        iter->c[0] = 0;
        iter->top = 0;
        return TOP->child[CTOP];
    }
    Node * current = TOP->child[CTOP];
    if(skip_children || 
       current->type == NODE_TYPE_LEAF || 
       current->nchildren == 0) {
        /* skip children or there is no children, visit sibling */
        if(CTOP + 1 < TOP->nchildren) {
            /* thee is a sibling */
            CTOP ++;
            return TOP->child[CTOP];
        }
        iter->top --;
        /* find first parent sibling */
        while(iter->top >= 0 && 
            CTOP == TOP->nchildren) {
            iter->top --;
        }
        /* exhausted */
        if(iter->top < 0) return NULL;
        return TOP->child[CTOP];
    } else {
        /* do the children */
        CTOP ++;
        iter->top ++;
        TOP = (InnerNode *) current;
        CTOP = 0;
        return TOP->child[CTOP];
    }
}
Node * tree_iter_next(TreeIter * iter) {
    return tree_iter_real_next(iter, FALSE);
}
Node * tree_iter_next_sibling(TreeIter * iter) {
    return tree_iter_real_next(iter, TRUE);
}
