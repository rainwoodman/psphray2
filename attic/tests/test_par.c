#include <glib.h>
#include "par.h"

typedef struct {
    float t;
} GasProp;

void test_creation();
void test_pack();
void test_remove();
void pstore_check_consistency(PStore * pstore);

int main(int argc, char * argv[]) {
    register_ptype(0, "gas", sizeof(GasProp), FALSE);
    register_ptype(1, "dark", 0, TRUE);
    test_creation();
    test_pack();
    test_remove();
    return 0;
}

PStore * pstore;

void test_creation() {
    pstore = pstore_new(4);
    g_random_set_seed(0);
    int i;
    int NPAR = 100000;
    g_message("creating a pstore and populating");
    for(i = 0; i < NPAR; i++) {
        ipos_t ipos[3];
        ipos[0] = g_random_int_range(0, IPOS_LIMIT);
        ipos[1] = g_random_int_range(0, IPOS_LIMIT);
        ipos[2] = g_random_int_range(0, IPOS_LIMIT);
        pstore_insert_par(pstore, ipos, 0);
        pstore_insert_par(pstore, ipos, 1);
    }
    pstore_check_consistency(pstore);
    g_message("added %td particles", pstore->root->size);
}

void test_pack() {
    size_t NPAR = pstore->root->size;
    int i;
    g_message("obtainng the pack for root");
    PackedPar * pack = pstore_pack_create_from_node(pstore->root);
    g_message("sorting");
    pstore_pack_sort(pack);
    PackedPar * sel = pstore_pack_create_from_selection(pack, 0, pack->size);
    for(i = 0; i < pstore->root->size; i++) {
        g_assert(ipos_compare(
            pstore_pack_get(pack, i)->ipos,
            pstore_pack_get(sel, i)->ipos) == 0);

        Par * par = pstore_pack_get(pack, i);
        if(i != pstore->root->size - 1) {
            g_assert(ipos_compare(par->ipos, 
                    pstore_pack_get(pack, i + 1)->ipos) <= 0);
        }
    }
    pstore_pack_free(sel);
    pstore_pack_free(pack);
}
void test_remove() {
    Par * p = pstore->root->first;
    g_message("removing first par at " PAR_FMT, PAR_PRINT(p[0]));
    g_message("total length = %d", g_slist_length((GSList *) p));
    pstore_remove_par(pstore, p);
    g_message("removed length = %d", g_slist_length((GSList *)p));
    p = pstore->root->first;
    g_message("total length = %d", g_slist_length((GSList *)p));
    g_message(PAR_FMT, PAR_PRINT(p[0])); 

    g_message("removing some more points");
    size_t NPAR = pstore->root->size * 0.5;
    int i;
    for(i = 0; i < NPAR; i++) {
        int ind = g_random_int_range(0, pstore->root->size);
        Par * p = pstore_get_nearby_par(pstore, ind);
        pstore_remove_par(pstore, p);
        pstore_free_par(p);
    }
    pstore_check_consistency(pstore);
}

static void pstore_check_r(PStore * pstore, Node * node, int depth, ipos_t x, ipos_t y, ipos_t z, Node ** last_ext);
void pstore_check_consistency(PStore * pstore) {
    Node * last_ext = NULL;
    pstore_check_r(pstore, pstore->root, 0, 0, 0, 0, &last_ext);
}

static void pstore_check_r(PStore * pstore, Node * node, int depth, ipos_t x, ipos_t y, ipos_t z, Node ** last_ext) {
    ipos_t width = 1 << (IPOS_NBITS - depth);
    ipos_t m = 1 << (IPOS_NBITS - depth - 1);
    if(node->last) {
        g_assert(node->last->next == pstore_node_next_par(node));
    }
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
        pstore_check_r(pstore, node->link[0], depth + 1, x, y, z, last_ext);
        pstore_check_r(pstore, node->link[1], depth + 1, x + m, y, z, last_ext);
        pstore_check_r(pstore, node->link[2], depth + 1, x, y + m, z, last_ext);
        pstore_check_r(pstore, node->link[3], depth + 1, x + m, y + m, z, last_ext);
        pstore_check_r(pstore, node->link[4], depth + 1, x, y, z + m, last_ext);
        pstore_check_r(pstore, node->link[5], depth + 1, x + m, y, z + m, last_ext);
        pstore_check_r(pstore, node->link[6], depth + 1, x, y + m, z + m, last_ext);
        pstore_check_r(pstore, node->link[7], depth + 1, x + m, y + m, z + m, last_ext);
    } else {
        int prefix;
        for(prefix = 0; prefix < 6; prefix++) {
            g_assert(node->link[prefix] == NULL);
        }
        g_assert(node->primary_size <= pstore->split_limit);
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
        if(*last_ext) {
            g_assert((*last_ext)->next == node);
            g_assert(node->prev == (*last_ext));
        }
        *last_ext = node;
    }
}

