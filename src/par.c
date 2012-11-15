#include <glib.h>
#include <mpi.h>
#include <stdlib.h>
#include "commonblock.h"

par_t * _PAR = NULL;
size_t NPAR = 0;
size_t _PAR_size = 0;
static par_t * _PARbase = NULL;

par_t * _PARin = NULL;
size_t NPARin = 0;

static intptr_t searchsorted (void * key, void * array, 
      size_t len, size_t elsize, GCompareFunc compare);

void par_allocate_input(size_t size) {
    _PARin = g_malloc(size* sizeof(par_t));
}

void par_free_input() {
    g_free(_PARin);
    _PARin = NULL;
}

void par_allocate(size_t size, size_t before) {
    _PARbase = g_malloc((size + before)* sizeof(par_t));
    _PAR_size = size;
    _PAR = _PARbase + before;
}

    static int _par_compare_fckey(par_t * par1, par_t * par2) {
        return fckey_cmp(&par1->fckey, &par2->fckey);
    }

void par_sort_by_fckey(int which) {
    if(which == PAR_BUFFER_IN) {
        qsort(_PARin, NPARin, sizeof(par_t), (GCompareFunc) _par_compare_fckey);
    } else {
        qsort(_PAR, NPAR, sizeof(par_t), (GCompareFunc) _par_compare_fckey);
    }
}

    static int _par_compare_fckey2(par_t * par, fckey_t * target) {
        return fckey_cmp(&par->fckey, target);
    }

intptr_t par_search_by_fckey(fckey_t * key, int which) {
    if(which == PAR_BUFFER_IN) {
        return searchsorted(key, _PARin, NPARin, sizeof(par_t), 
        (GCompareFunc) _par_compare_fckey2);
    } else {
        return searchsorted(key, _PAR, NPAR, sizeof(par_t), 
        (GCompareFunc) _par_compare_fckey2);

    }
}

size_t par_grow(size_t addsize) {
    if(NPAR + addsize <= _PAR_size) {
        NPAR += addsize;
    } else {
        g_error("not enough space in PARbase\n");
    }
    return NPAR - addsize;
}

void par_shift(size_t addsize) {
    if(addsize <= _PAR - _PARbase) {
        _PAR -= addsize;
        _PAR_size += addsize;
    } else {
        g_error("not enough space in PARbase\n");
    }
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

    return middle;
}
