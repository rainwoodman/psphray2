#include <glib.h>
#include <mpi.h>
#include <string.h>
#include "commonblock.h"

static Node * tree_store_append(TreeStore * store, int type) {
    Node * rt;
    if(type == NODE_TYPE_LEAF) {
        rt = (Node*) stack_push(&store->leaf, Node);
    } else {
        rt = (Node*) stack_push(&store->inner, InnerNode);
    }
    rt->type = type;
    return rt;
}
static Node * tree_store_pop(TreeStore * store, Node * node) {
     if(node->type == NODE_TYPE_LEAF) {
        if(node == (Node*) stack_peek(&store->leaf, Node)) {
            return (Node*) stack_pop(&store->leaf, Node);
        } else {
            g_critical("trying to pop a non top leaf node");
            return NULL;
        }
    } else {
        if(node == (Node*) stack_peek(&store->inner, InnerNode)) {
            return (Node*) stack_pop(&store->inner, InnerNode);
        } else {
            g_critical("trying to pop a non top inner node");
            return NULL;
        }
    }

}
void tree_store_init(TreeStore * store) {
    stack_init(&store->inner, InnerNode);
    stack_init(&store->leaf, Node);
}
void tree_store_destroy(TreeStore * store) {
    stack_destroy(&store->leaf);
    stack_destroy(&store->inner);
}

int tree_node_locate_fckey(Node * node, fckey_t * key) {
    /* returns -1 if not in the node, otherwise returns
     * the octrant(0-7) that contains the key */
    fckey_t tmp;
    fckey_xor(&tmp, &node->key, key);
    fckey_rightshift(&tmp, 3 * (node->order - 1));
    if(tmp.a[1] == 0 && tmp.a[0] < 8) {
        return tmp.a[0];
    } else {
        return -1;
    }
}

int tree_node_contains_fckey(Node * node, fckey_t * key) {
    return tree_node_locate_fckey(node, key) != -1;
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

static InnerNode * node_to_inner(TreeStore * store, Node * node) {
    if(node->type == NODE_TYPE_LEAF) {
        InnerNode * rt = (InnerNode*) tree_store_append(store, NODE_TYPE_INNER);
        memcpy(rt, node, sizeof(Node));
        /* memcpy will overwrite type */
        rt->type = NODE_TYPE_INNER;
        /* now replace the reference to the node in the parent */
        InnerNode * parent = node->parent;
        if(parent) {
            int pos = tree_node_locate_fckey((Node*)parent, &rt->key);
            parent->child[pos] = (Node *)rt;
        }
        tree_store_pop(store, node);
        return rt;
    } else {
        g_error("already inner");
    }
}

static Node * create_leaf(TreeStore * store, Node * parent, par_t * first) {
    if(parent->type == NODE_TYPE_LEAF) {
        parent = (Node*) node_to_inner(store, parent);
    }
    Node * rt = tree_store_append(store, NODE_TYPE_LEAF);
    rt->first = first;
    rt->npar = 0;
    rt->order = parent->order - 1;
    rt->key = first->fckey;
    fckey_clear(&rt->key, 3 * rt->order);
    rt->parent = (InnerNode*)parent;
    int pos = tree_node_locate_fckey((Node*)rt->parent, &rt->key);
    rt->parent->child[pos] = rt;
    return rt;
}

static void tree_build_subtree(TreeStore * store, Node * subroot, par_t * first, intptr_t npar) {
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
            node = create_leaf(store, node, i);
        } else /* must be a leaf node */
        if(node->npar > CB.NodeSplitThresh && node->order > 0) {
            /* split the node */
            i = node->first;
            node = create_leaf(store, node, i);
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

Node * tree_build(TreeStore * store, PSystem * psys) {
    /* first make the root */
    Node * root = (Node *) tree_store_append(store, NODE_TYPE_INNER);
    root->order = FCKEY_BITS;
    root->first = &psys->data[0];
    tree_build_subtree(store, root, &psys->data[0], psys->length);
    return root;
}

Node * tree_locate_fckey(Node * root, fckey_t * key) {
    TreeIter iter;
    tree_iter_init(&iter, root);
    Node * node = tree_iter_next(&iter);
    while(node) {
        if(tree_node_contains_fckey(node, key)) {
            if(node->type == NODE_TYPE_LEAF) {
                break;
            } else {
                node = tree_iter_next(&iter);
            }
        }  else {
            node = tree_iter_next_sibling(&iter);
        }
    }
    return node;
}

void tree_iter_init(TreeIter * iter, Node * root) {
    iter->root = root;
    iter->current = NULL;
}

static Node * tree_iter_real_next(TreeIter * iter, int skip_children) {
    if(iter->current == NULL) {
        iter->current = iter->root;
        return iter->current;
    } else {
        if(! skip_children 
        && iter->current->type == NODE_TYPE_INNER
        ) {
            int i = 0;
            while(i < 8) {
                if(((InnerNode *) iter->current)->child[i]) {
                    iter->current = ((InnerNode *) iter->current)->child[i];
                    //g_message("visiting child");
                    return iter->current;
                }
                i++;
            }
        } /* fall through */
        /* visit sibling */
        InnerNode * parent = iter->current->parent;
        while(iter->current != iter->root) {
            int i = tree_node_locate_fckey((Node *)parent, 
                & iter->current->key) + 1;
            while(i < 8) {
                if(parent->child[i] != NULL) {
                    iter->current = parent->child[i];
                    return iter->current;
                }
                i++;
            }
            iter->current = (Node*) parent;
            parent = parent->parent;
        }
        iter->current = NULL;
        return NULL;
    }
}
Node * tree_iter_next(TreeIter * iter) {
    return tree_iter_real_next(iter, FALSE);
}
Node * tree_iter_next_sibling(TreeIter * iter) {
    return tree_iter_real_next(iter, TRUE);
}
