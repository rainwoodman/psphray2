#include <glib.h>
#include <mpi.h>

#include "commonblock.h"
struct header_t{
    unsigned int N[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    unsigned int Ntotlow[6];
    int flag_cool;
    int Nfiles;
    double boxsize;
    double OmegaM;
    double OmegaL;
    double h;
    int flag_sft;
    int flag_met;
    unsigned int Ntothigh[6];
    int flag_entropy;
    int flag_double;
    int flag_ic_info;
    int flag_lpt_scalingfactor;
    int unused[12];
};

static MPI_Comm ReaderComm = 0;

static int ReaderRank; /* 0 if needs to read from this task */
static int ReaderSize; /* total number of tasks this reader will scatter to */
static int ReaderColor; /* the rank of the reader group relative to all CB.NReadder*/
static int Nfile;
static int NReader;

#define LEADERONLY if(ReaderRank == 0)
static int get_fmt_field_count(char * fmt) {
    GError * error = NULL;
    GRegex * regex = g_regex_new("%[^d]*d", (GRegexCompileFlags)0, (GRegexMatchFlags)0, &error);
    if(!regex) {
        g_error("regex is wrong: %s\n", error->message);
    }
    GMatchInfo * matchinfo = NULL;
    g_regex_match(regex, CB.inputbase, (GRegexMatchFlags) 0, &matchinfo);
    int count = 1;
    while(g_match_info_next(matchinfo, NULL)) count ++;
    g_match_info_free(matchinfo);
    g_regex_unref(regex);
    return count;
}

static char * format_input_filename(int major, int fid) {
    switch(get_fmt_field_count(CB.inputbase)) {
        case 1:
            return g_strdup_printf(CB.inputbase, major);
        case 2:
            return g_strdup_printf(CB.inputbase, major, fid);
        case 3:
            return g_strdup_printf(CB.inputbase, major, major, fid);
        default:
            g_error("format of inputbase is wrong, need 3, 2 or 1 %%d");
            return NULL;
    }
}

static void snapshot_read_header(int fid, SnapHeader * h) {
    char * filename = format_input_filename(CB.SnapNumMajor, fid);
    GError * error = NULL;
    GMappedFile * file = g_mapped_file_new(filename, FALSE, &error);
    if(!file) {
        g_error("cannot open %s:%s", filename, error->message);
    }
    g_free(filename);
    char * buffer = g_mapped_file_get_contents(file);
    struct header_t * fh = (struct header_t * )( buffer + 4);
    for(int i = 0; i < 6; i++) {
        h[0].N[i]= fh[0].N[i];
        h[0].masstab[i]= fh[0].mass[i];
    }
    h[0].Ngas = fh[0].N[0];
    h[0].NgasTotal = fh[0].Ntotlow[0] 
        + (((int64_t) fh[0].Ntothigh[0]) << 32);
    h[0].a = fh[0].time;
    h[0].BoxSize = fh[0].boxsize;
    h[0].h = fh[0].h;
    h[0].Nfile = fh[0].Nfiles;
    h[0].double_precision = fh[0].flag_double;
    g_mapped_file_unref(file);
}

static size_t snapshot_prepare() {
    /* returns number of pars will be read into this task */
    SnapHeader h;
    ROOTONLY {
        snapshot_read_header(0, &h);
        if(CB.C.h != h.h) {
            g_warning("snapshot file h mismatch");
        }
    }
    MPI_Bcast(&h, sizeof(SnapHeader), MPI_BYTE, 0, MPI_COMM_WORLD);

    Nfile = h.Nfile;
    NReader = CB.NReader;

    CB.BoxSize = h.BoxSize;
    if(NReader > h.Nfile) {
        NReader = h.Nfile;
        ROOTONLY {
            g_warning("NReader bigger than files, using Nfile");
        }
    }

    if(ReaderComm) MPI_Comm_free(&ReaderComm);
    MPI_Comm_split(MPI_COMM_WORLD, ThisTask % NReader, 
           ThisTask / NReader, &ReaderComm);
    ReaderColor = ThisTask % NReader;
    ReaderRank = ThisTask / NReader;
    MPI_Comm_size(ReaderComm, &ReaderSize);
    MPI_Barrier(MPI_COMM_WORLD);

    int64_t Ngas = 0;
    LEADERONLY {
        int fid;
        for(fid = ReaderColor * Nfile / NReader;
            fid < (ReaderColor + 1) * Nfile / NReader;
            fid ++) {
            SnapHeader h;
            snapshot_read_header(fid, &h);
            Ngas += h.Ngas;
            g_debug("fid %d has Ngas %ld", fid, Ngas);
        }
    }
    MPI_Bcast(&Ngas, sizeof(Ngas), MPI_BYTE, 0, ReaderComm);
    int64_t start = ReaderRank * Ngas / ReaderSize;
    int64_t end = (ReaderRank + 1) * Ngas / ReaderSize;

    return end - start;
}

static int snapshot_scatter_block(char * data, size_t N, par_t * recvbuf, 
       int eloffset, int elsize) {
    /* return number of particles shall be received to
     * this task, 
     * if data is not NULL, dispatch the data to recvbuf, given at
     * eloffset and elsize
     * 
     * after all blocks are scattered, we increase NPAR
     * */
    int recvcount = 0;
    recvcount = ((ReaderRank + 1) * N / ReaderSize 
              - ReaderRank * N / ReaderSize);
    if(recvbuf == NULL) return recvcount;

    int * sendcounts = NULL;
    int * displs = NULL;
    MPI_Datatype sendtype;
    MPI_Datatype recvtype, tmprecvtype;

    MPI_Type_contiguous(elsize, MPI_BYTE, &sendtype);
    MPI_Type_commit(&sendtype);
    int one = 1;
    MPI_Aint eloffsetA = eloffset;
    MPI_Type_create_struct(1, &one, &eloffsetA, &sendtype, &tmprecvtype);
    MPI_Type_create_resized(tmprecvtype, 0, sizeof(par_t), &recvtype);
    MPI_Type_commit(&recvtype);

    LEADERONLY {
        sendcounts = g_new0(int, ReaderSize);
        displs = g_new0(int, ReaderSize);
        int i;
        for(i = 0; i < ReaderSize; i++) {
            int start = i * N / ReaderSize;
            int end = (i + 1) * N / ReaderSize;
            displs[i] = start;
            sendcounts[i] = (end - start);
        }
    }

    MPI_Scatterv(data, sendcounts, displs, sendtype, 
              recvbuf, recvcount, recvtype, 0, ReaderComm);

    
    MPI_Type_free(&tmprecvtype);
    MPI_Type_free(&recvtype);
    MPI_Type_free(&sendtype);

    LEADERONLY {
        g_free(sendcounts);
        g_free(displs);
    }
    return recvcount;
}

static void cast_to_float_t(char * in, char * out, size_t N, int elsize) {
    intptr_t i;
    /* this will ensure we do not cast on non LEADERs */
    if(out == NULL) return;
    if(elsize == 8) {
        for(i = 0; i < N; i++) {
            ((float_t *)out)[i] = ((double*)in)[i];
        }
    } else {
        for(i = 0; i < N; i++) {
            ((float_t *)out)[i] = ((float*)in)[i];
        }
    }
}
static void cast_to_long_t(char * in, char * out, size_t N, int elsize) {
    intptr_t i;
    /* this will ensure we do not cast on non LEADERs */
    if(out == NULL) return;
    if(elsize == 8) {
        for(i = 0; i < N; i++) {
            ((long_t*)out)[i] = ((int64_t*)in)[i];
        }
    } else {
        for(i = 0; i < N; i++) {
            ((long_t*)out)[i] = ((int32_t*)in)[i];
        }
    }

}
static void convert_to_fckey_t(char * in, char * out, size_t N, int elsize) {
    if(out == NULL) return;
    int64_t ipos[3];
    for(intptr_t i = 0; i < N; i++) {
        if(elsize == 8) {
            for(int d=0; d < 3; d++) {
                ipos[d] = ((double*)in)[i* 3 + d] / CB.BoxSize * FCKEY_MAX;
            }
        } else {
            for(int d=0; d < 3; d++) {
                ipos[d] = ((float*)in)[i* 3 + d] / CB.BoxSize * FCKEY_MAX;
            }
        }
        fckey_from_ipos(&((fckey_t *)out)[i], ipos);
    }
}

static void snapshot_read_one(int fid, SnapHeader * h) {
    LEADERONLY {
        snapshot_read_header(fid, h);
    }
    MPI_Bcast(h, sizeof(SnapHeader), MPI_BYTE, 0, MPI_COMM_WORLD);

    size_t Ntot = 0;
    for(int i = 0; i < 6; i++) Ntot += h[0].N[i];

    int snap_float_elsize = sizeof(float);

    if(h[0].double_precision) {
        snap_float_elsize = sizeof(double);
    }
    
    char * buffer = NULL;
    GMappedFile * file = NULL;
    char * p = NULL;
    char * cast = NULL;

    LEADERONLY {
        char * filename = format_input_filename(CB.SnapNumMajor, fid);
        g_message("%d (%d/%d on %d) reading %s", ThisTask, ReaderRank, ReaderSize, ReaderColor, filename);
        GError * error = NULL;
        file = g_mapped_file_new(filename, FALSE, &error);
        if(!file) {
            g_error("failed to read file %s:%s", filename, error->message);
        }
        buffer = g_mapped_file_get_contents(file);
        g_free(filename);
        cast = g_malloc(MAX(sizeof(float_t) * 3, sizeof(long_t)) * h[0].Ngas);
    }

    size_t END = NPARin;

    NPARin += snapshot_scatter_block(NULL, h[0].Ngas, NULL, 0, 0);

    TAKETURNS {
        g_message("%02d NPARin = %ld ", ThisTask, NPARin);
    }
    /* now buffer is the file on LEADER */
    p = buffer + 256 + 4 + 4; /* skip header */

    /*pos*/
    convert_to_fckey_t(p + 4, cast, h[0].Ngas, snap_float_elsize);
    snapshot_scatter_block(cast, h[0].Ngas, & PARin(END), 
        offsetof(par_t, fckey), elsizeof(par_t, fckey));
    p += snap_float_elsize * 3 * Ntot + 4 + 4;
    MPI_Barrier(ReaderComm);
    /* vel */
    p += snap_float_elsize * 3 * Ntot + 4 + 4;

    /* id */
    cast_to_long_t(p + 4, cast, h[0].Ngas, CB.IDByteSize);
    snapshot_scatter_block(cast, h[0].Ngas, & PARin(END), 
          offsetof(par_t, id), elsizeof(par_t, id));
    p += CB.IDByteSize * Ntot + 4 + 4;

    /* mass */
    /* mass tab handling */
    if (h[0].masstab[0] == 0) 
        cast_to_float_t(p + 4, cast, h[0].Ngas, snap_float_elsize);
    else
        for(int i = 0; i < h[0].Ngas; i++) {
            ((float_t *) cast)[i] = h[0].masstab[0];
        }
    snapshot_scatter_block(cast, h[0].Ngas, & PARin(END), 
        offsetof(par_t, mass), elsizeof(par_t, mass));
    p += 4 + 4;
    /* mass tab handling */
    for(int i=0; i< 6; i++) {
        if (h[0].masstab[i] == 0) p += snap_float_elsize * h[0].N[i];
    }

    /* ie */
    cast_to_float_t(p + 4, cast, h[0].Ngas, snap_float_elsize);
    snapshot_scatter_block(cast, h[0].Ngas, & PARin(END), 
        offsetof(par_t, T), elsizeof(par_t, T));
    p += snap_float_elsize * h[0].Ngas + 4 + 4;

    LEADERONLY {
        g_free(cast);
        g_mapped_file_unref(file);
        g_message("finished with file %d", fid);
    }
}

void snapshot_read() {
    size_t Ngas = snapshot_prepare();
    NPARin = 0;
    par_allocate_input(Ngas * 1.1);
    int fid;
    for(fid = ReaderColor * Nfile / NReader;
        fid < (ReaderColor + 1) * Nfile / NReader;
        fid ++) {
        SnapHeader h;
        snapshot_read_one(fid, &h);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}


