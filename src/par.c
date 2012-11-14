#include <glib.h>
#include <mpi.h>

#include "commonblock.h"

GArray * _PAR = NULL;

void par_preallocate(size_t size) {
    if(_PAR != NULL) {
        g_array_set_size(_PAR, size);
        g_array_set_size(_PAR, 0);
    } else {
        _PAR = g_array_sized_new(FALSE, FALSE, sizeof(par_t), size);
    }
}

void par_increase_size(size_t addsize) {
    g_array_set_size(_PAR, _PAR->len + addsize);
}

