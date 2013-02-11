#include <glib.h>
#include <string.h>

#include "par.h"

static struct {
    char name[8];
    size_t elesize;
} PTYPE[256];

/* register a ptype. elesize is the size of extra storage, which 
 * can be accessed via Par.data */
int register_ptype(int id, char * name, size_t elesize) {
    int i;
    for(i = 0; i < 256; i++) {
        if(PTYPE[i].name[0] == 0) {
            strncpy(PTYPE[i].name, name, 8);
            PTYPE[i].name[7] = 0;
            PTYPE[i].elesize = elesize;
            return i;
        }
    }
}

extern inline int ipos_get_prefix(ipos_t ipos[3], int depth);
extern inline int par_is_primary(PStore * pstore, Par * par);

PStore * pstore_new(int primarymask) {
    PStore * store = g_new0(PStore, 1);
    store->root = g_slice_new0(Node);
    store->primarymask = primarymask;
    store->merge_limit = 128;
    store->split_limit = 128;
    store->depth_limit = 19;
//IPOS_NBITS;
    return store;
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

static void pstore_insert_r(PStore * pstore, Node * node, Par * par, int depth) {
    node->size ++;
    node->primary_size += par_is_primary(pstore, par);
    if(node->link[0]) {
        /* has children, move in */
        int prefix = ipos_get_prefix(par->ipos, depth);
        //g_message("diving in ");
        int update_par = node->link[prefix]->first == node->first;
        int update_last = node->link[prefix]->last == node->last;
        pstore_insert_r(pstore, node->link[prefix], par, depth + 1);
        /* if last and par has changed, and node->last and par were same,
         * update them */
        if(update_par) node->first = node->link[prefix]->first;
        if(update_last) node->last = node->link[prefix]->last;
    } else {
        if(par_is_primary(pstore, par) 
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
        } else {
            /* append the par */
            if(node->last == NULL) {
                node->first = par;
            } else {
                node->last->next = par;
            }
            node->last = par;
            Par * previous_par = pstore_node_previous_par(node);
            if(previous_par) {
                previous_par->next = par;
            }
            Par * next_par = pstore_node_next_par(node);
            node->last->next = next_par;
            //g_message("direct insert %p into %p", par, node);
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
}

static void pstore_remove_r(PStore * pstore, Node * node, Par * par, int depth) {
    if(node->link[0]) {
        /* inner */
        int prefix = ipos_get_prefix(par->ipos, depth);
        int update_par = node->link[prefix]->first == node->first;
        int update_last = node->link[prefix]->last == node->last;
        pstore_remove_r(pstore, node->link[prefix], par, depth + 1);
        if(update_par) node->first = node->link[prefix]->first;
        if(update_last) node->last = node->link[prefix]->last;
        node->size --;
        node->primary_size -= par_is_primary(pstore, par);
    } else {
        /* external */
        Par * p = NULL;
        if(node->first == par) {
            node->first = par->next;
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
        node->primary_size -= par_is_primary(pstore, par);
        Par * previous_par = pstore_node_previous_par(node);
        Par * next_par = pstore_node_next_par(node);
        if(previous_par) {
            if(node->first) {
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

static Par * pstore_get_r(PStore * pstore, Node * node, ptrdiff_t index);
Par * pstore_get(PStore * pstore, ptrdiff_t index) {
    return pstore_get_r(pstore, pstore->root, index);
}
static Par * pstore_get_r(PStore * pstore, Node * node, ptrdiff_t index) {
    if(node->link[0]) {
        int prefix;
        ptrdiff_t offset = 0;
        for(prefix = 0; prefix < 8; prefix++) {
            if(index - offset < node->link[prefix]->size) {
                return pstore_get_r(pstore, node->link[prefix], index - offset);
            }
            offset += node->link[prefix]->size;
        }
        /* if this fails, index is out of range */
        g_assert_not_reached();
    }  else {
        Par * par;
        for(par = node->first; index > 0; index --, par = par->next)
            continue;
        return par;
    }
}

/**
 * finish up a pstore, merge the nodes that are too small.
 * also link all particles 
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


