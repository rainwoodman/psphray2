#include <glib.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include "mpiu.h"
int ThisTask;
int NTask;
int PrevTask, NextTask;
MPI_Datatype MPI_PTRDIFF;

void abort() {
    void exit(int);
    MPI_Abort(MPI_COMM_WORLD, 1);
    /* shut up the compiler */
    exit(1);
}

void print_handler(const gchar * string) {
    fputs(string, stdout);
    fflush(stdout);
}

static void log_handler
    (const gchar *log_domain,
    GLogLevelFlags log_level,
    const gchar *message,
    gpointer unused_data) {

    g_log_default_handler(log_domain, log_level, message, unused_data);
    if((log_level & G_LOG_FLAG_FATAL) && NTask > 1) {
        g_on_error_stack_trace ("");
        abort();
    }
}


void mpiu_module_init() {
    g_log_set_default_handler(log_handler, NULL);
    g_set_print_handler(print_handler);
    if(sizeof(ptrdiff_t) == sizeof(long long)) {
        MPI_PTRDIFF = MPI_LONG_LONG;
    } else if(sizeof(ptrdiff_t) == sizeof(long)) {
        MPI_PTRDIFF = MPI_LONG;
    } else if(sizeof(ptrdiff_t) == sizeof(int)) {
        MPI_PTRDIFF = MPI_INT;
    } else {
        g_assert_not_reached();
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    PrevTask = (ThisTask - 1 + NTask) % NTask,
    NextTask = (ThisTask + 1) % NTask;
}

void mpiu_bcast_string(char ** string) {
    int len;
    ROOTONLY {
        if(string[0]) {
            len = strlen(string[0]);
        } else {
            len = -1;
        }
    }

    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(len == -1) {
        string[0] = NULL;
    } else {
        ROOTONLY { } else {
            string[0] = g_new(char, len + 1);
        }
        MPI_Bcast(string[0], len + 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
}
