#include <glib.h>
#include <mpi.h>
#include <stdlib.h>
#include "commonblock.h"

PSystem PAR_IN = { "INPUT", 0 };

static intptr_t searchsorted (void * key, void * array, 
      size_t len, size_t elsize, GCompareFunc compare);

void par_destroy(PSystem * psys) {
    g_free(psys->base);
    psys->data = NULL;
    psys->size = 0;
    psys->length = 0;
}

void par_reserve(PSystem * psys, size_t size, size_t before) {
    psys->base = g_malloc((size + before)* sizeof(par_t));
    psys->size = size;
    psys->data = psys->base + before;
}

    static int _par_compare_fckey(par_t * par1, par_t * par2) {
        return fckey_cmp(&par1->fckey, &par2->fckey);
    }

void par_sort_by_fckey(PSystem * psys) {
    qsort(psys->data, psys->length, 
            sizeof(par_t), (GCompareFunc) _par_compare_fckey);
}

    static int _par_compare_fckey2(par_t * par, fckey_t * target) {
        return fckey_cmp(&par->fckey, target);
    }

intptr_t par_search_by_fckey(PSystem * psys, fckey_t * key) {
    return searchsorted(key, psys->data, psys->length, 
        sizeof(par_t), (GCompareFunc) _par_compare_fckey2);
}

par_t * par_append(PSystem * psys, intptr_t add) {
    /* returns a pointer where you can put the new data in */
    /* warning any index reference to the psys will be invalid,
     * if add is negative, truncates the psys and returns 
     * a pointer to the truncated data */
    if(psys->length + add <= psys->size) {
        psys->length += add;
    } else {
        g_error("not enough space on rear needed %ld has %ld\n",
             add, psys->size - psys->length);
    }
    if(add > 0) {
        return psys->data + psys->length - add;
    } else {
        return psys->data + psys->length;
    }
}

par_t * par_prepend(PSystem * psys, intptr_t add) {
    /* if add is negative, truncates the psys and returns 
     * a pointer to the truncated data */
    if(add <= psys->data - psys->base) {
        psys->data -= add;
        psys->length += add;
        psys->size += add;
    } else {
        g_error("not enough space on front need %ld has %ld\n",
          add, psys->data - psys->base);
    }
    g_message(" new par base = %p %p", psys->data, psys->data + add);
    if(add > 0) return psys->data;
    else return psys->data + add;
}

void par_update_igm(par_t * i) {
    i->IGMmass = i->mass * (1.0 - get_cloud_fraction(i->rho * (CB.a * CB.a * CB.a)));
}

static intptr_t searchsorted (void * target, void * array, 
      size_t len, size_t elsize, GCompareFunc compare) {

    intptr_t left;
    intptr_t middle;
    intptr_t right;
    intptr_t cmpval;

    char * data = array;
    g_return_val_if_fail(array != NULL, -1);
    g_return_val_if_fail(compare != NULL, -1);

    left = 0;
    right = len - 1;

    cmpval = compare(target, data + (elsize * right));
    if(cmpval > 0) {
        return len;
    }

    cmpval = compare(target, data);
    if(cmpval < 0) {
        return 0;
    }

    while (left <= right) {
        middle = (left + right) / 2;
        cmpval = compare(target, data + (elsize * middle));
        if (cmpval > 0) {
            left = middle + 1;
        } else if (cmpval < 0) {
            right = middle - 1;
        } else {
            return middle;
        }
    }

    return left;
}
