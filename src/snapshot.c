#include <glib.h>
#include <mpi.h>

#include "commonblock.h"
#include "fckey.h"
#include "par.h"
#include "snapshot.h"
PSystem * PAR_BUFFER_IN;

static MPI_Comm ReaderComm = 0;

static int ReaderRank; /* 0 if needs to read from this task */
static int ReaderSize; /* total number of tasks this reader will scatter to */
static int ReaderColor; /* the rank of the reader group relative to all CB.NReadder*/
static int Nfile;
static int NReader;

#define LEADERONLY if(ReaderRank == 0)

/*
 * count the %d field numbers in a string fmt
 * used to smartly decide the input file name format
 * */
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

/*
 * smartly format a input format file.
 * major is the snapshot number and fid is the file number within a snapshot.
 *
 * %d/%d.%d -> major, major, fid
 * %d.%d -> major, fid,
 * %d -> major
 * */
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

/* read header of file fid and current major snapshot number into h
 * h is a generic format, indippend from header_t */
static void snapshot_read_header(int fid, SnapHeader * h) {
    char * filename = format_input_filename(CB.SnapNumMajor, fid);
    GError * error = NULL;
    GMappedFile * file = g_mapped_file_new(filename, FALSE, &error);
    if(!file) {
        g_error("cannot open %s:%s", filename, error->message);
    } g_free(filename);
    char * buffer = g_mapped_file_get_contents(file);
    struct {
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
    } * fh = buffer + 4;

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

static intptr_t snapshot_prepare(SnapHeader * h) {
    ROOTONLY {
        snapshot_read_header(0, h);
        if(CB.C.h != h->h) {
            g_warning("snapshot file h mismatch");
        }
        if(CB.BoxSize != h->BoxSize) {
            CB.BoxSize = h->BoxSize;
            g_warning("boxsize mismatch!, using %g from snapshot", h->BoxSize);
        }
    }
    MPI_Bcast(h, sizeof(SnapHeader), MPI_BYTE, 0, MPI_COMM_WORLD);

    Nfile = h->Nfile;
    if(CB.NReader == 0) CB.NReader = NTask;
    if(CB.NReader > NTask) {
        ROOTONLY g_warning("too many readers requested, using %d", NTask);
        CB.NReader = NTask;
    }
    if(CB.NReader > h->Nfile) {
        CB.NReader = h->Nfile;
        ROOTONLY g_warning("NReader bigger than files, using Nfile");
    }
    NReader = CB.NReader;

    if(ReaderComm) MPI_Comm_free(&ReaderComm);
    MPI_Comm_split(MPI_COMM_WORLD, ThisTask % NReader, 
           ThisTask / NReader, &ReaderComm);
    ReaderColor = ThisTask % NReader;
    ReaderRank = ThisTask / NReader;
    MPI_Comm_size(ReaderComm, &ReaderSize);
    MPI_Barrier(MPI_COMM_WORLD);

    /* decide how many particles will be written to this task */
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
    MPI_Bcast(&Ngas, 1, MPI_LONG, 0, ReaderComm);
    int64_t start = ReaderRank * Ngas / ReaderSize;
    int64_t end = (ReaderRank + 1) * Ngas / ReaderSize;

    return end - start;
}

/*
 * scatter a data block (size N) from data to par_t within the local reader group,
 * eloffset is the offset in recvbuf, and elsize is the size of an element.
 *
 * if recvbuf == NULL, then return the number of particles will be scattered to
 * this task, so that we can allocate enough bytes.
 *
 * the subroutine uses MPI_Datatype API to scatter the data.
 * */
static int snapshot_scatter_block(char * data, size_t N, par_t * recvbuf, 
       int eloffset, int elsize) {
    int64_t recvcount = 0;
    recvcount = ((ReaderRank + 1) * N / ReaderSize 
              - ReaderRank * N / ReaderSize);
    if(recvcount > G_MAXINT32) {
        g_error("more than 2G particles cannot be sent with MPI");
    }
    if(recvbuf == NULL) {
        return recvcount;
    }
    /* send type is easy */
    MPI_Datatype sendtype;
    MPI_Type_contiguous(elsize, MPI_BYTE, &sendtype);
    MPI_Type_commit(&sendtype);

    /* recv type is messy,
     * first the an element at the given offset, 
     * then grow the size of the struct to match the receive type */
    MPI_Datatype recvtype, tmprecvtype;
    int one = 1;
    MPI_Aint eloffsetA = eloffset;
    MPI_Type_create_struct(1, &one, &eloffsetA, &sendtype, &tmprecvtype);
    MPI_Type_create_resized(tmprecvtype, 0, sizeof(par_t), &recvtype);
    MPI_Type_commit(&recvtype);

    int * sendcounts = NULL;
    int * displs = NULL;

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

/* cast a float kind (double/float) to internal float_t */
static void cast_to_float_t(char * in, char * out, size_t N, int elsize) {
    intptr_t i;
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
/* cast a float kind (double/float) to internal long_t */
static void cast_to_long_t(char * in, char * out, size_t N, int elsize) {
    intptr_t i;
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
static void convert_to_fckey_t(char * in, char * out, size_t N, int elsize, double BoxSize) {
    int64_t ipos[3];
    for(intptr_t i = 0; i < N; i++) {
        if(elsize == 8) {
            for(int d=0; d < 3; d++) {
                ipos[d] = ((double*)in)[i* 3 + d] / BoxSize * (FCKEY_MAX + 1);
            }
        } else {
            for(int d=0; d < 3; d++) {
                ipos[d] = ((float*)in)[i* 3 + d] / BoxSize * (FCKEY_MAX + 1);
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

    par_t * recvbuf = par_append(PAR_BUFFER_IN,  
         snapshot_scatter_block(NULL, h[0].Ngas, NULL, 0, 0));

    /* now buffer is the file on LEADER */
    LEADERONLY {
        p = buffer + 256 + 4 + 4; /* skip header */
    }
    /*pos*/
    LEADERONLY {
        convert_to_fckey_t(p + 4, cast, h[0].Ngas, snap_float_elsize, h[0].BoxSize);
        p += snap_float_elsize * 3 * Ntot + 4 + 4;
    }
    snapshot_scatter_block(cast, h[0].Ngas, recvbuf, 
        offsetof(par_t, fckey), elsizeof(par_t, fckey));
    /* vel */
    LEADERONLY {
        p += snap_float_elsize * 3 * Ntot + 4 + 4;
    }
    /* do not read vel */
    /* id */
    LEADERONLY {
        cast_to_long_t(p + 4, cast, h[0].Ngas, CB.IDByteSize);
        p += CB.IDByteSize * Ntot + 4 + 4;
    }
    snapshot_scatter_block(cast, h[0].Ngas, recvbuf, 
          offsetof(par_t, id), elsizeof(par_t, id));
    /* mass */
    /* mass tab handling */
    LEADERONLY {
        if (h[0].masstab[0] == 0) 
            cast_to_float_t(p + 4, cast, h[0].Ngas, snap_float_elsize);
        else
            for(int i = 0; i < h[0].Ngas; i++) {
                ((float_t *) cast)[i] = h[0].masstab[0];
            }
        int massblocklength = 0;
        /* mass tab handling */
        for(int i=0; i< 6; i++) {
            if (h[0].masstab[i] != 0) massblocklength += h[0].N[i];
        }
        if(massblocklength > 0) {
            p += snap_float_elsize * massblocklength + 4 + 4;
        }
    }
    snapshot_scatter_block(cast, h[0].Ngas, recvbuf, 
        offsetof(par_t, mass), elsizeof(par_t, mass));
    /* ie */
    LEADERONLY {
        cast_to_float_t(p + 4, cast, h[0].Ngas, snap_float_elsize);
        p += snap_float_elsize * h[0].Ngas + 4 + 4;
    }
    snapshot_scatter_block(cast, h[0].Ngas, recvbuf, 
        offsetof(par_t, T), elsizeof(par_t, T));
    /* rho */
    LEADERONLY {
        cast_to_float_t(p + 4, cast, h[0].Ngas, snap_float_elsize);
        p += snap_float_elsize * h[0].Ngas + 4 + 4;
    }
    snapshot_scatter_block(cast, h[0].Ngas, recvbuf, 
        offsetof(par_t, rho), elsizeof(par_t, rho));
    /* ye */
    LEADERONLY {
        cast_to_float_t(p + 4, cast, h[0].Ngas, snap_float_elsize);
        p += snap_float_elsize * h[0].Ngas + 4 + 4;
    }
    #if 0
    snapshot_scatter_block(cast, h[0].Ngas, recvbuf, 
        offsetof(par_t, ye), elsizeof(par_t, ye));
    #endif

    LEADERONLY {
        g_free(cast);
        g_mapped_file_unref(file);
    }
}

void snapshot_read(SnapHeader * h) {
    if(PAR_BUFFER_IN == NULL) {
        PAR_BUFFER_IN = par_alloc();
        par_init(PAR_BUFFER_IN, "INPUT");
    }
    size_t Ngas = snapshot_prepare(h);
    par_reserve(PAR_BUFFER_IN, Ngas * 1.1, 0);
    int fid;
    for(fid = ReaderColor * Nfile / NReader;
        fid < (ReaderColor + 1) * Nfile / NReader;
        fid ++) {
        SnapHeader hh;
        snapshot_read_one(fid, &hh);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    ParIter iter;

    for(par_t * i = par_iter_init(&iter, PAR_BUFFER_IN); i; i = par_iter_next(&iter)) {
        i->IGMmass = i->mass * (1.0 - get_cloud_fraction(i->rho * (h->a * h->a * h->a)));
    }
}


