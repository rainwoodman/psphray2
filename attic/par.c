#include <glib.h>
#include <string.h>

#include "par.h"

static struct {
    char name[8];
    size_t elesize;
} PTYPE[256];

unsigned int PRIMARY_MASK = 0;
/* register a ptype. elesize is the size of extra storage, which 
 * can be accessed via Par.data */
void register_ptype(int i, char * name, size_t elesize, int is_primary) {
    strncpy(PTYPE[i].name, name, 8);
    PTYPE[i].name[7] = 0;
    PTYPE[i].elesize = elesize;
    if(is_primary) {
        PRIMARY_MASK |= (1L << i);
    }
}

extern inline int ipos_get_prefix(ipos_t ipos[3], int depth);
extern inline int par_is_primary(Par * par);

PStore * pstore_new(size_t split_limit) {
    PStore * store = g_new0(PStore, 1);
    store->root = g_slice_new0(Node);
    store->merge_limit = -1; /* unused */
    store->split_limit = split_limit;
    store->depth_limit = IPOS_NBITS - 1;
    return store;
}

void par_free(Par * par) {
    g_slice_free1(PTYPE[par->type].elesize + sizeof(Par), par);
}
void par_free_chain(Par * head) {
    Par * q = NULL;
    for(; head; head = q) {
        q = head->next;
        g_slice_free1(PTYPE[head->type].elesize + sizeof(Par), head);
    }
}
/**
 * insert a particle at given position to the store, of type ptype
 * returns a pointer to the particle so that the information can be filled
 */
static void pstore_insert_r(PStore * pstore, Node * node, Par * par, int depth);
Par * pstore_insert(PStore * pstore, ipos_t ipos[3], int ptype) {
    Par * par = g_slice_alloc0(sizeof(Par) + PTYPE[ptype].elesize);
    par->ipos[0] = ipos[0];
    par->ipos[1] = ipos[1];
    par->ipos[2] = ipos[2];
    pstore_insert_r(pstore, pstore->root, par, 0);
    return par;
}
static int ipos_compare(ipos_t a[3], ipos_t b[3]) {
    for(int depth = 0; depth < IPOS_NBITS; depth ++ ) {
        int prefix_a = ipos_get_prefix(a, depth);
        int prefix_b = ipos_get_prefix(b, depth);
        if (prefix_a < prefix_b) {
            return -1;
        } 
        if (prefix_a > prefix_b) {
            return 1;
        }
    }
    return 0;
}
static int child_is_first_nonempty(Node * node, int prefix) {
    int i;
    for(i = prefix - 1; i >= 0; i--) {
        if(node->link[i]->size > 0) return 0;
    }
    return 1;
}
static int child_is_last_nonempty(Node * node, int prefix) {
    int i;
    for(i = prefix + 1; i < 8; i++) {
        if(node->link[i]->size > 0) return 0;
    }
    return 1;
}
static void pstore_insert_r(PStore * pstore, Node * node, Par * par, int depth) {
    node->size ++;
    node->primary_size += par_is_primary(par);
    if(node->link[0]) {
        /* has children, move in */
        int prefix = ipos_get_prefix(par->ipos, depth);
        //g_message("diving in ");
        pstore_insert_r(pstore, node->link[prefix], par, depth + 1);
        /* if last and par has changed, and node->last and par were same,
         * update them */
        if(prefix <= node->first_nonempty_child) {
            node->first_nonempty_child = prefix;
            node->first = node->link[prefix]->first;
        }
        if(prefix >= node->last_nonempty_child) {
            node->last_nonempty_child = prefix;
            node->last = node->link[prefix]->last;
        }
//        g_assert(g_slist_position(node->first, node->last) == node->size - 1);
    } else {
        if(par_is_primary(par) 
        && node->primary_size + 1 > pstore->split_limit && depth < pstore->depth_limit) {
            /* primary par may cause a split */
            int prefix;
            for(prefix = 0; prefix < 8; prefix++) {
                node->link[prefix] = g_slice_new0(Node);
                node->link[prefix]->up = node;
                node->link[prefix]->prefix = prefix;
            }
            //g_message("splitting node %p, for %d", node, node->primary_size);
            node->size = 0;
            node->primary_size = 0;
            node->first_nonempty_child = 8;
            node->last_nonempty_child = -1;
            Par * oldpar = node->first;
            node->last->next= NULL;
            node->first = NULL;
            node->last = NULL;
            Par * p = oldpar, * q;
            for(; p; p = q) {
                q = p->next;
             //   g_message("split reinserting %p into %p", p, node);
                pstore_insert_r(pstore, node, p, depth);
            }
            pstore_insert_r(pstore, node, par, depth);
         //   g_assert(g_slist_position(node->first, node->last) == node->size - 1);
        } else {
            //g_message("node size = %d %d", node->size, node->primary_size);
            int update_first = 0;
            int update_last = 0;
            /* append the par */
            if(node->last == NULL) {
                /* empty list */
                node->first = par;
                node->last = par;
                update_first = 1;
                update_last = 1;
            } else {
                Par * p = node->first;
                Par * pp = NULL;
#ifdef ENABLE_SORT_NODE
                while(p && ipos_compare(par->ipos, p->ipos) >= 0) {
                    pp = p;
                    /* next par must be with greater than par*/
                    g_assert(pp != node->last->next);
                    p = p->next;
                }
#endif
                if(pp == NULL) {
                    /* prepend */
                    par->next = p;
                    node->first = par;
                    update_first = 1;
                } else {
                    /* insert after pp */
                    pp->next = par;
                    par->next = p;
                }
                if(node->last == pp) {
                    node->last = par;
                    update_last = 1;
                }
            }
         //   g_assert(g_slist_position(node->first, node->last) == node->size - 1);
            if(update_first) {
                Par * previous_par = pstore_node_previous_par(node);
                if(previous_par) {
                    previous_par->next = node->first;
                }
            }
            if(update_last) {
                Par * next_par = pstore_node_next_par(node);
                node->last->next = next_par;
            }
//            g_message("direct insert %p into %p", par, node);
 //           g_message("%p %p", node->first, node->last);
        }
    }
}


/**
 * remove a particle by pointer ptr,
 * crash if par is not found
 * */
static void pstore_remove_r(PStore * pstore, Node * node, Par * par, int depth);
void pstore_remove(PStore * pstore, Par * par) {
    pstore_remove_r(pstore, pstore->root, par, 0);
    par->next = NULL;
}

static void pstore_remove_r(PStore * pstore, Node * node, Par * par, int depth) {
    if(node->link[0]) {
        /* inner */
        int prefix = ipos_get_prefix(par->ipos, depth);
        g_assert(prefix >= node->first_nonempty_child);
        g_assert(prefix <= node->last_nonempty_child);
        pstore_remove_r(pstore, node->link[prefix], par, depth + 1);
        if(prefix == node->first_nonempty_child) {
            int i = prefix;
            while(i < 8 && node->link[i]->size == 0) {
                i ++;
            }
            node->first_nonempty_child = i;
            if(i < 8) node->first = node->link[i]->first;
        }
        if(prefix == node->last_nonempty_child) {
            int i = prefix;
            while(i >= 0 && node->link[i]->size == 0) {
                i --;
            }
            node->last_nonempty_child = i;
            if(i >= 0) node->last = node->link[i]->last;
        }
        node->size --;
        node->primary_size -= par_is_primary(par);
        if(node->last_nonempty_child == -1 ||
           node->first_nonempty_child == 8 ||
           node->size == 0
        ) {
            g_assert(node->last_nonempty_child == -1);
            g_assert(node->first_nonempty_child == 8);
            g_assert(node->size == 0);
            /* this node is empty and shall no longer be internal */
            int i;
            for(i = 0; i < 8; i++) {
                g_assert(node->link[i]->size == 0);
                g_slice_free(Node, node->link[i]);
                node->link[i] = NULL;
            }
            node->first = NULL;
            node->last = NULL;
        }
    } else {
        /* external */
        Par * p = NULL;
        if(node->first == par) {
            node->first = par->next;
            if(node->size == 1) node->first = NULL;
        } else {
            for(p = node->first; p && p->next != par; p = p->next) continue;
            p->next = par->next;
            /* if this fails, par is not in this node, and the tree must
             * ve been corrupted */
            g_assert(p != NULL);
        }
        if(par == node->last) {
            node->last = p;
        }
        node->size --;
        node->primary_size -= par_is_primary(par);
        Par * previous_par = pstore_node_previous_par(node);
        Par * next_par = pstore_node_next_par(node);
        if(previous_par) {
            if(node->size) {
                previous_par->next = node->first;
            } else {
                previous_par->next = next_par;
            }
        } 
        if(node->last) {
            node->last->next = next_par;
        }
    }
}

Par * pstore_node_previous_par(Node * node) {
    if(node->up == NULL) return NULL;
    int prefix;
    /* NOTE: we cannot use node->first == node->up->first to 
     * quickly decide if the node is the first non-empty child.
     * when this function is called,
     * node may have been emptied in case there has been a split,
     * but the parent nodes are not modified */
    for(prefix = node->prefix - 1; prefix >= 0; prefix --) {
        Par * child_last_par = node->up->link[prefix]->last;
        if(child_last_par) return child_last_par;
    }
    return pstore_node_previous_par(node->up);
}

Par * pstore_node_next_par(Node * node) {
    if(node->up == NULL) return NULL;
    int prefix;
    /* NOTE: we cannot use node->first == node->up->first to 
     * quickly decide if the node is the first non-empty child.
     * when this function is called,
     * node may have been emptied in case there has been a split,
     * but the parent nodes are not modified */
    for(prefix = node->prefix + 1; prefix < 8; prefix ++) {
        Par * child_first_par = node->up->link[prefix]->first;
        if(child_first_par) return child_first_par;
    }
    return pstore_node_next_par(node->up);
}

/*
 * return the pointer to a particle at the index-th position
 * */
static Par * pstore_get_nearby_r(PStore * pstore, Node * node, ptrdiff_t index);
Par * pstore_get_nearby(PStore * pstore, ptrdiff_t index) {
    return pstore_get_nearby_r(pstore, pstore->root, index);
}
static Par * pstore_get_nearby_r(PStore * pstore, Node * node, ptrdiff_t index) {
    if(node->link[0]) {
        int prefix;
        ptrdiff_t cumsum = 0;
        for(prefix = 0; prefix < 8; prefix++) {
            ptrdiff_t cumsum2 = cumsum + node->link[prefix]->size;
            if(index >= cumsum && index < cumsum2) {
                return pstore_get_nearby_r(pstore, node->link[prefix], index - cumsum);
            }
            cumsum = cumsum2;
        }
        /* if this fails, index is out of range */
        g_assert_not_reached();
    }  else {
        g_assert(index >= 0 && index < node->size);
        Par * par = node->first;
//#ifdef ENABLE_EXACT_GET_NEARBY
        for(; index > 0; index --, par = par->next)
            continue;
//#endif
        return par;
    }
}

/**
 * finish up a pstore, merge the nodes that are too small.
 * also link all particles;
 *
 * unused function.
 */
static void pstore_merge_r(PStore * pstore, Node * node) {
    /* this merging process is probably not that useful.
     * in density calculation we can always check the size of 
     * the children and if they are too small, do it over the parent.
     * */
    short int prefix;
    short int mergethis = 0;
    /* first see if there are children to merge at all*/
    if(!node->link[0]) return;
    /* then see if need to merge (any children is too small */
    for(prefix = 0; prefix < 8; prefix++) {
        if(node->link[prefix]->primary_size  < pstore->merge_limit 
        && node->link[prefix]->size > 0) {
           /* merge only if there are particles in the node,
            * pure empty cells spans empty space.
            * */
            mergethis = 1;
            break;
        }
    }
    if(!mergethis) {
        for(prefix = 0; prefix < 8; prefix++) {
            pstore_merge_r(pstore, node->link[prefix]);
        }
    } else {
        for(prefix = 0; prefix < 8; prefix++) {
            g_slice_free(Node, node->link[prefix]);
            node->link[prefix] = NULL;
        }
    }
}

/**
 * pack particles into a packedpar struct
 *
 * the next pointers are invalid but we do not bother to clear them.
 * 
 */
PackedPar * pstore_pack(Par * first, size_t size) {
    size_t bytes = 0;
    Par * par;
    ptrdiff_t i = 0;
    ptrdiff_t N[256] = {0};
    for(par = first, i = 0; i < size; i++, par = par->next) {
        N[par->type] ++;
    }
    for(i = 0; i < 256; i++) {
        bytes += PTYPE[i].elesize * N[i];
    }
    bytes += size * sizeof(Par);
    PackedPar * pack = g_malloc(sizeof(PackedPar) + bytes);
    pack->bytes = bytes;
    pack->size = size;
    char * ptr = pack->data;
    for(par = first, i = 0; i < size; i++, par = par->next) {
        size_t width = sizeof(Par) + PTYPE[par->type].elesize;
        memcpy(ptr, par, width);
        ptr += width;
    }
    return pack;
}

/**
 * unpack packed par to a GSList kind of structure.
 *
 * this is done in place. do not free any elements. free the entire pack
 * with one g_free, or the free of whatever corresponding allocator
 * */

Par * pstore_unpack(PackedPar * pack) {
    ptrdiff_t i;
    char * ptr = pack->data;
    Par * par = (Par * ) ptr;
    for(i = 0; i < pack->size; i++) {
        size_t width = sizeof(Par) + PTYPE[par->type].elesize;
        ptr += width;
        par->next = (Par *) ptr;
        par = par->next;
    }
    /* remember to set the last next pointer to NULL, it used to point
     * to outside of the data buffer */
    par->next = NULL;
    /* return the first par */
    return (Par *) pack->data;
}

static void pstore_check_r(PStore * pstore, Node * node, int depth, ipos_t x, ipos_t y, ipos_t z);
void pstore_check(PStore * pstore) {
    pstore_check_r(pstore, pstore->root, 0, 0, 0, 0);
}
static void pstore_check_r(PStore * pstore, Node * node, int depth, ipos_t x, ipos_t y, ipos_t z) {
    ipos_t width = 1 << (IPOS_NBITS - depth);
    ipos_t m = 1 << (IPOS_NBITS - depth - 1);
    if(node->link[0]) {
        int prefix;
        for(prefix = 0; prefix < 8; prefix++) {
            if(prefix < node->first_nonempty_child)
                g_assert(node->link[prefix]->size == 0);
            if(prefix > node->last_nonempty_child)
                g_assert(node->link[prefix]->size == 0);
            if(prefix == node->first_nonempty_child)
                g_assert(node->link[prefix]->size > 0);
            if(prefix == node->last_nonempty_child)
                g_assert(node->link[prefix]->size > 0);
            g_assert(node->link[prefix]);
        }
        pstore_check_r(pstore, node->link[0], depth + 1, x, y, z);
        pstore_check_r(pstore, node->link[1], depth + 1, x + m, y, z);
        pstore_check_r(pstore, node->link[2], depth + 1, x, y + m, z);
        pstore_check_r(pstore, node->link[3], depth + 1, x + m, y + m, z);
        pstore_check_r(pstore, node->link[4], depth + 1, x, y, z + m);
        pstore_check_r(pstore, node->link[5], depth + 1, x + m, y, z + m);
        pstore_check_r(pstore, node->link[6], depth + 1, x, y + m, z + m);
        pstore_check_r(pstore, node->link[7], depth + 1, x + m, y + m, z + m);
    } else {
        int prefix;
        for(prefix = 0; prefix < 8; prefix++) {
            g_assert(node->link[prefix] == NULL);
        }
        g_assert(node->size <= pstore->split_limit);
        Par * par = node->first;
        int i;
        if(node->size > 0) {
            g_assert(node->first != NULL);
            g_assert(node->last != NULL);
            Par * prev = pstore_node_previous_par(node);
            if(prev)
                g_assert(prev->next == node->first);
            Par * next = pstore_node_next_par(node);
            g_assert(node->last->next == next);
        }
        size_t pc = 0;
        for(i = 0; i < node->size; i++) {
            g_assert(par != NULL);
            pc += par_is_primary(par);
            if(i == node->size - 1) {
                g_assert(par == node->last);
            }
            g_assert(x <= par->ipos[0]);
            g_assert(y <= par->ipos[1]);
            g_assert(z <= par->ipos[2]);
            g_assert(x + width > par->ipos[0]);
            g_assert(y + width > par->ipos[1]);
            g_assert(z + width > par->ipos[2]);
            par = par->next;
        }
        g_assert(pc == node->primary_size);
    }
    
}
