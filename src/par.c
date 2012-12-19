#include <glib.h>
#include <stdint.h>
#include <stdlib.h>
#include "fckey.h"
#include "par.h"

struct _PSystem {
    char name[8];
    par_t * base;
    par_t * data;
    size_t length;
    size_t size;
};



static intptr_t searchsorted (void * key, void * array, 
      size_t len, size_t elsize, GCompareFunc compare);

void par_init(PSystem * psys, char * name) {
    strncpy(psys->name, name, 7);
}

PSystem * par_alloc() {
    return g_slice_new0(PSystem);
}
void par_free(PSystem * psys) {
    g_slice_free(PSystem, psys);
}

void par_destroy(PSystem * psys) {
    if(psys->base) g_free(psys->base);
    psys->base = NULL;
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
par_t * par_index(PSystem * psys, intptr_t i) {
    if(i >= 0) return &psys->data[i];
    else      return &psys->data[psys->length + i];
}
size_t par_get_length(PSystem * psys) {
    return psys->length;
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
    if(add > 0) return psys->data;
    else return psys->data + add;
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

    while (left < right) {
        middle = left + ((right-left) >> 1);
        cmpval = compare(target, data + (elsize * middle));
        if (cmpval > 0) {
            left = middle + 1;
            /* left element is always less than target*/
        } else {
            right = middle;
        }
    }
    return left;
}

par_t * par_iter_init(ParIter * iter, PSystem * psys) {
    iter->psys = psys;
    iter->i = 0;
    iter->end = psys->length;
    return par_iter_set(iter, 0);
}
par_t * par_iter_init_range(ParIter * iter, PSystem * psys, intptr_t first, size_t npar) {
    iter->psys = psys;
    iter->end = first + npar;
    return par_iter_set(iter, first);
}

par_t * par_iter_next(ParIter * iter) {
    if(iter->i < iter->end) {
        iter->last_i = iter->i;
        iter->i = iter->i + 1;
        return &iter->psys->data[iter->last_i];
    }
    return NULL;
}
/* reset the iter to a given position */
par_t * par_iter_set(ParIter * iter, intptr_t index) {
    iter->last_i = index;
    iter->i = index + 1;
    return &iter->psys->data[iter->last_i];
}

intptr_t par_iter_last_index(ParIter * iter) {
    return iter->last_i;
}
