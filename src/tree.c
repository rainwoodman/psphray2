#include <stdint.h>
#include <glib.h>
#include "fckey.h"
#include "par.h"
#include "tree.h"

#include "stack.h"
struct _TreeStore {
    Stack inner;
    Stack leaf;
    PSystem * psys;
    int splitthresh;
};

static void tree_build_subtree(TreeStore * store, Node * subroot, intptr_t ifirst, intptr_t npar, intptr_t * skipped);

/* fetch a node from the tree store,
 * for i > 0, get innernode
 * for i < 0, get leafnode
 * for i ==0, crash */
inline Node * tree_store_get_node(TreeStore * store, intptr_t i) {
    if(i > 0) {
        return (Node*) stack_get(&store->inner, i - 1, InnerNode);
    }
    if(i < 0) {
        return stack_get(&store->leaf, -1 - i , Node);
    }
    return NULL;
}

Node * tree_store_get_leaf_nodes(TreeStore * store, intptr_t *nleaf) {
    if(nleaf) *nleaf = store->leaf.len;
    return (Node*) store->leaf.data;
}

/* the reverse operation of tree_store_get_node
 * returns the index of a given node */
inline intptr_t tree_store_get_index(TreeStore * store, Node * node) {
    if(node->type != NODE_TYPE_INNER) {
        return -1 - (node - (Node*) store->leaf.data);
    } else {
        return 1 + ((InnerNode*)node - (InnerNode*) store->inner.data);
    }
}



/* create a new node of given type, and return the pointer */
static inline Node * tree_store_append(TreeStore * store, int type) {
    Node * rt;
    if(type != NODE_TYPE_INNER) {
        rt = (Node*) stack_push(&store->leaf, Node);
    } else {
        rt = (Node*) stack_push(&store->inner, InnerNode);
    }
    rt->type = type;
    return rt;
}

/* remove a node of given type. the node has to be the last
 * node on the stack, otherwise the stack won't shrink
 * references to the node shall be cleaned before this is called.
 * */
static inline Node * tree_store_pop(TreeStore * store, Node * node) {
     if(node->type != NODE_TYPE_INNER) {
        if(node == (Node*) stack_peek(&store->leaf, Node)) {
            return (Node*) stack_pop(&store->leaf, Node);
        } else {
//            g_critical("trying to pop a non top leaf node");
            node->status = NODE_FREE;
            return NULL;
        }
    } else {
        if(node == (Node*) stack_peek(&store->inner, InnerNode)) {
            return (Node*) stack_pop(&store->inner, InnerNode);
        } else {
//           g_critical("trying to pop a non top inner node");
            node->status = NODE_FREE;
            return NULL;
        }
    }
}

TreeStore * tree_store_alloc() {
    return g_slice_new0(TreeStore);
}
void tree_store_free(TreeStore * store) {
    g_slice_free(TreeStore, store);
}
/* initialize a tree store and build the tree for given psys*/
void tree_store_init(TreeStore * store, PSystem * psys, int splitthresh) {
    stack_init(&store->inner, InnerNode);
    stack_init(&store->leaf, Node);
    store->splitthresh = splitthresh;
    store->psys = psys;
}

/* build a tree and returns number of particles
 * skipped because they are out of box */
intptr_t tree_build(TreeStore * store) {
    intptr_t skipped = 0;
    Node * root = (Node *) tree_store_append(store, NODE_TYPE_INNER);
    root->order = FCKEY_BITS;
    root->ifirst = 0;
    tree_build_subtree(store, root, 0, par_get_length(store->psys), &skipped);
    return skipped;

}
void tree_store_destroy(TreeStore * store) {
    stack_destroy(&store->leaf);
    stack_destroy(&store->inner);
    store->psys = NULL;
}

/* returns the child node position that contains the key
 * or -1 if it is not contained.
 * */
int tree_node_childpos_fckey(Node * node, fckey_t * key) {
    /* returns -1 if not in the node, otherwise returns
     * the octrant(0-7) that contains the key */
    fckey_t tmp;
    fckey_xor(&tmp, &node->key, key);
    if(node->order == 0) {
        g_error("node->order shall be positive");
    }
    fckey_rightshift(&tmp, 3 * (node->order - 1));
    if(tmp.a[1] == 0 && tmp.a[0] < 8) {
        return tmp.a[0];
    } else {
        return -1;
    }
}

/* test if key is contained in node 
 * returns 1 if key is contained in node
 * 0 if not.
 * */
int tree_node_contains_fckey(Node * node, fckey_t * key) {
    if(node->order > 0) {
        return tree_node_childpos_fckey(node, key) != -1;
    } else {
        return fckey_cmp(&node->key, key) == 0;
    }
}

/* test if needle is a subnode of node 
 * returns 1 if needle and node are identical
 * or needle is a child of node */
int tree_node_contains_node(Node * node, Node * needle) {
    return node->order >= needle->order 
       && tree_node_contains_fckey(node, & needle->key);
}

/* returns the first parent of node that is identical to needle
 * identical means both key and order match */
Node * tree_locate_up_image(TreeStore * store, Node * node, Node * needle) {
    Node * image = node;
    while(image) {
        if(fckey_cmp(&image->key, &needle->key) == 0 && 
           image->order == needle->order) {
            return image;
        }
        image = (Node *) tree_store_get_node(store, image->iparent);
    }
    return NULL;
}


/* converts a node to a inner node.
 * the index of the node may change, but this routine takes care
 * of maintaining the reference in the tree from the parent
 * */
static inline InnerNode * node_to_inner(TreeStore * store, Node * node) {
    if(node->type != NODE_TYPE_INNER) {
        InnerNode * rt = (InnerNode*) tree_store_append(store, NODE_TYPE_INNER);
        *((Node*) rt) = *node;
        /* the copy will overwrite type */
        rt->type = NODE_TYPE_INNER;
        /* now replace the reference to the node in the parent */
        InnerNode * parent = (InnerNode*) tree_store_get_node(store, rt->iparent);
        if(parent) {
            int pos = tree_node_childpos_fckey((Node*)parent, &rt->key);
            parent->ichild[pos] = tree_store_get_index(store, (Node*)rt);
        }
        tree_store_pop(store, node);
        return rt;
    } else {
        g_error("already inner");
    }
}

/* add a leaf to the given parent. if parent is not an inner node
 * convert it with node_to_inner.
 * also update the parent's child list to include this newly
 * created leaf node.
 *
 * ifirst and npar will be junk.
 * we also assume the parent node is a InnerNode.
 * */
static inline Node * create_leaf_at(TreeStore * store, InnerNode * parent, int childpos) {
    Node * rt = tree_store_append(store, NODE_TYPE_LEAF);
    rt->order = parent->order - 1;
    rt->key = parent->key;
    fckey_or_with_leftshift(&rt->key, childpos, 3 * rt->order);
    rt->iparent = tree_store_get_index(store, (Node*) parent);
    (parent)->ichild[childpos] = tree_store_get_index(store, rt);
    return rt;
}
static inline Node * create_leaf(TreeStore * store, Node * parent, par_t * i) {
    if(parent->type != NODE_TYPE_INNER) {
        parent = (Node*) node_to_inner(store, parent);
    }
    fckey_t key = i->fckey;
    int pos = tree_node_childpos_fckey(parent, &key);
    Node * rt = create_leaf_at(store, (InnerNode*) parent, pos);
    return rt;
}

Node * tree_split_empty_leaf(TreeStore * store, Node * parent) {
    if(parent->npar != 0) {
        g_error("only split empty leaf nodes!");
    }
    if(parent->type != NODE_TYPE_INNER) {
        parent = (Node*) node_to_inner(store, parent);
    }
    int pos;
    for(pos = 0; pos < 8; pos++) {
        Node * child = create_leaf_at(store, (InnerNode*) parent, pos);
        child->ifirst = parent->ifirst;
        child->npar = 0;
    }
    return parent;
}

/* build a tree as the subtree of subroot, with particles
 * first to first + npar.
 * returns the number of particles that are out of the box
 * in skipped.
 * */
static void tree_build_subtree(TreeStore * store, Node * subroot, intptr_t ifirst, intptr_t npar, intptr_t * skipped) {
    Node * exterior = tree_store_get_node(store, subroot->iparent);
    ParIter iter;
    par_t * i = par_iter_init_range(&iter, store->psys, ifirst, npar);
    intptr_t SKIPPED = 0;
    Node * node = subroot;
    /* now scan over PAR, creating nodes */
    while(i) {
        /* the current particle is no longer in current node,
         * time to close parent nodes */
        while(node !=  exterior && 
            ! tree_node_contains_fckey(node, 
                &i->fckey)) {
            node->npar = par_iter_last_index(&iter) - node->ifirst + 1;
            /* here we shall try to merge the children, if 
             * doing so is intended. not merging any for now */
            node = tree_store_get_node(store, node->iparent);
        }
        if(node == exterior) {
            /* particle is not in any nodes, skip it */
            SKIPPED ++;
            i = par_iter_next(&iter);
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
            node->ifirst = par_iter_last_index(&iter);
            node->npar = 0;
            goto add_particles;
        } 
        /* must be a leaf node */
        if(node->npar > store->splitthresh && node->order > 0) {
            /* rewind, split the node */
            i = par_iter_set(&iter, node->ifirst);
            node = create_leaf(store, node, i);
            node->ifirst = par_iter_last_index(&iter);
            node->npar = 0;
            /* fall through add the particles */
        }
        /* must be a leaf and no split is required,
         * code will fall through */
        add_particles:
        /* add particle to the leaf */

        /* TODO: figure a bigger step size: no need to scan
         * every particle  */
        node->npar +=  1;
        i = par_iter_next(&iter);
    }
    /* close the parent nodes of the last scanned particle */
    while(node != exterior) {
        node->npar = par_iter_last_index(&iter) - node->ifirst + 1;
        node = tree_store_get_node(store, node->iparent);
    }
    if(SKIPPED && skipped) {
        *skipped = SKIPPED;
    }
}

Node * tree_locate_down_fckey(TreeStore * store, Node * start, fckey_t * key) {
    TreeIter iter;
    Node * node = tree_iter_init(&iter, store, start);
    while(node) {
        if(tree_node_contains_fckey(node, key)) {
            if(node->type != NODE_TYPE_INNER) {
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

/* create leaf nodes that terminate the dangling children.*/
void tree_terminate(TreeStore * store) {
    for(intptr_t i = 0; i < store->inner.len; i++) {
        InnerNode * parent = (InnerNode*) tree_store_get_node(store, i + 1);
        for(int j = 0; j < 8; j++) {
            if(parent->ichild[j] != 0) continue;
            Node * node = create_leaf_at(store, parent, j);
            node->ifirst = parent->ifirst;
            node->npar = 0;
        }
    }
}

/* call, if eg, as we go the current node has been splited */
void tree_iter_update_current(TreeIter * iter, Node * node) {
    iter->current = node;
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
                if(((InnerNode *) iter->current)->ichild[i]) {
                    iter->current = tree_store_get_node(iter->store, 
                           ((InnerNode *) iter->current)->ichild[i]);
                    //g_message("visiting child");
                    return iter->current;
                }
                i++;
            }
        } /* fall through */
        /* visit sibling */
        InnerNode * parent = (InnerNode*) tree_store_get_node(iter->store, 
                        iter->current->iparent);
        while(iter->current != iter->root) {
            int i = tree_node_childpos_fckey((Node *)parent, 
                & iter->current->key) + 1;
            while(i < 8) {
                if(parent->ichild[i]) {
                    iter->current = tree_store_get_node(iter->store, 
                            parent->ichild[i]);
                    return iter->current;
                }
                i++;
            }
            iter->current = (Node*) parent;
            parent = (InnerNode *)tree_store_get_node(iter->store, parent->iparent);
        }
        iter->current = NULL;
        return NULL;
    }
}

Node * tree_iter_init(TreeIter * iter, TreeStore * store, Node * root) {
    iter->store = store;
    iter->root = root;
    iter->current = NULL;
    if(iter->root == NULL) {
        iter->root = tree_store_get_node(iter->store, 1);
    }
    return tree_iter_real_next(iter, 0);
}
Node * tree_iter_next(TreeIter * iter) {
    return tree_iter_real_next(iter, FALSE);
}
Node * tree_iter_next_sibling(TreeIter * iter) {
    return tree_iter_real_next(iter, TRUE);
}
