#include <glib.h>
#include "mpiu.h"

#include "commonblock.h"
#include "par.h"
#include "hydro.h"
#include "snapshot.h"


static MPI_Comm ReaderComm = 0;

static int ReaderRank; /* 0 if needs to read from this task */
static int ReaderSize; /* total number of tasks this reader will scatter to */
static int ReaderColor; /* the rank of the reader group relative to all CB.NReadder*/
static int Nfile;
static int NReader;
extern size_t ptype_get_elesize(int ptype);

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
void snapshot_read_header(int fid, SnapHeader * h) {
    char * filename = format_input_filename(CB.SnapNumMajor, fid);
    GError * error = NULL;
    GMappedFile * file = g_mapped_file_new(filename, FALSE, &error);
    if(!file) {
        g_error("cannot open %s:%s", filename, error->message);
    } g_free(filename);
    char * buffer = g_mapped_file_get_contents(file);
    struct _ {
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
    } * fh = (struct _ * ) (buffer + 4);

    for(int i = 0; i < 6; i++) {
        h[0].N[i]= fh[0].N[i];
        h[0].masstab[i]= fh[0].mass[i];
    }
    h[0].a = fh[0].time;
    h[0].BoxSize = fh[0].boxsize;
    h[0].h = fh[0].h;
    h[0].Nfile = fh[0].Nfiles;
    h[0].double_precision = fh[0].flag_double;
    g_mapped_file_unref(file);
}

static void snapshot_prepare(SnapHeader * h, ptrdiff_t Nglobal[6], int * fidstart, int * fidend) {
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
    int ptype;
    for(ptype = 0; ptype < 6; ptype++) {
        Nglobal[ptype] = 0;
    }
    LEADERONLY {
        *fidstart = ReaderColor * Nfile / NReader;
        *fidend = (ReaderColor + 1) * Nfile / NReader;
        int fid;
        for(fid = *fidstart;
            fid < *fidend;
            fid ++) {
            SnapHeader h;
            snapshot_read_header(fid, &h);
            for(ptype = 0; ptype < 6; ptype++) {
                Nglobal[ptype] += h.N[ptype];
            }
            g_debug("fid %d has Ngas %ld", fid, Nglobal[0]);
        }
    }
    MPI_Bcast(Nglobal, 6, MPI_LONG, 0, ReaderComm);
}

static int start_block(char ** p) {
    int rt;
    rt = *(int*) *p;
    *p += 4;
    return rt;
}
static void check_block(char **p, int expected) {
    int rt;
    rt = *(int*) *p;
    *p += 4;
    if(expected != rt) {
        g_error("block inconsistent");
    }
}
static real_t read_real(char ** p, int double_precision) {
    real_t rt;
    if(double_precision) {
        rt = * (double*) *p;
        *p += 8;
    } else {
        rt = * (float*) *p;
        *p += 4;
    }
    return rt;
}
static uint64_t read_id(char ** p, int longid) {
    uint64_t rt;
    if(longid) {
        rt = * (uint64_t*) *p;
        *p += 8;
    } else {
        rt = * (uint32_t*) *p;
        *p += 4;
    }
    return rt;
}

/*
 * Read a snapshot file fid into 6 PackedPar, one per ptype.
 * */
static void snapshot_to_pack(int fid, PackedPar * pack[6]) {
    GError * error = NULL;
    char * buffer = NULL;
    GMappedFile * file = NULL;
    char * p = NULL;
    char * cast = NULL;
    SnapHeader h;
    snapshot_read_header(fid, &h);

    char * filename = format_input_filename(CB.SnapNumMajor, fid);
    g_message("%d (%d/%d on %d) reading %s", ThisTask, ReaderRank, ReaderSize, ReaderColor, filename);
    file = g_mapped_file_new(filename, FALSE, &error);
    if(!file) g_error("failed to read file %s:%s", filename, error->message);

    buffer = g_mapped_file_get_contents(file);
    g_free(filename);

    ptrdiff_t N[6];
    int ptype;

    for(int ptype = 0; ptype < 6; ptype++) {
        N[ptype] = h.N[ptype];
        pack[ptype] = pstore_pack_create_simple(ptype, N[ptype], TRUE);
    }

    p = buffer + 256 + 4 + 4; /* skip header */
    ptrdiff_t i = 0;
    int check;
    if(N[0] + N[1] + N[2] + N[3] + N[4] + N[5] > 0) {
        check = start_block(&p);
        for(ptype = 0; ptype < 6; ptype++) {
            for(i = 0; i < N[ptype]; i++) {
                Par * par = pstore_pack_get(pack[ptype], i);
                par->ipos[0] = read_real(&p, h.double_precision) / CB.BoxSize * IPOS_LIMIT;
                par->ipos[1] = read_real(&p, h.double_precision) / CB.BoxSize * IPOS_LIMIT;
                par->ipos[2] = read_real(&p, h.double_precision) / CB.BoxSize * IPOS_LIMIT;
            }
        }
        check_block(&p, check);
        check = start_block(&p);
        for(ptype = 0; ptype < 6; ptype++) {
            for(i = 0; i < N[ptype]; i++) {
                Par * par = pstore_pack_get(pack[ptype], i);
                par->vel[0] = read_real(&p, h.double_precision);
                par->vel[1] = read_real(&p, h.double_precision);
                par->vel[2] = read_real(&p, h.double_precision);
            }
        }
        check_block(&p, check);
        check = start_block(&p);
        for(ptype = 0; ptype < 6; ptype++) {
            for(i = 0; i < N[ptype]; i++) {
                Par * par = pstore_pack_get(pack[ptype], i);
                par->id = read_id(&p, CB.IDByteSize == 8);
            }
        }
        check_block(&p, check);
    }

    if( N[0] * (h.masstab[0] != 0)
     +  N[1] * (h.masstab[1] != 0)
     +  N[2] * (h.masstab[2] != 0)
     +  N[3] * (h.masstab[3] != 0)
     +  N[4] * (h.masstab[4] != 0)
     +  N[5] * (h.masstab[5] != 0)
    ) {
        check = start_block(&p);
        for(ptype = 0; ptype < 6; ptype ++) {
            for(i = 0; i < N[ptype]; i++) {
                Par * par = pstore_pack_get(pack[ptype], i);
                if(h.masstab[ptype] == 0) {
                    par->mass = read_real(&p, h.double_precision);
                } else {
                    par->mass = h.masstab[ptype];
                }
            }
        }
        check_block(&p, check);
    }

#define READ(converter, element, ptype) \
    for(i = 0; i < N[ptype]; i++) {\
        Par * par = pstore_pack_get(pack[ptype], i); \
        converter(par)->element = read_real(&p, h.double_precision); \
    }
    if(N[0] > 0) {
        check = start_block(&p);
        READ(AS_GAS, ie, 0);
        check_block(&p, check);
        check = start_block(&p);
        READ(AS_GAS, rho, 0);
        check_block(&p, check);
        check = start_block(&p);
        READ(AS_GAS, ye, 0);
        check_block(&p, check);
        check = start_block(&p);
        READ(AS_GAS, xHI, 0);
        check_block(&p, check);
        check = start_block(&p);
        READ(AS_GAS, sml, 0);
        check_block(&p, check);
        check = start_block(&p);
        READ(AS_GAS, sfr, 0);
        check_block(&p, check);
    }
    if(N[4] > 0) {
        check = start_block(&p);
        READ(AS_STAR, sft, 4);
        check_block(&p, check);
    }
    if(N[0] + N[4] > 0) {
        check = start_block(&p);
        READ(AS_GAS, met, 0);
        READ(AS_STAR, met, 4);
        check_block(&p, check);
    }
    if(N[5] > 0) {
        check = start_block(&p);
        READ(AS_BH, bhmass, 5);
        check_block(&p, check);
        check = start_block(&p);
        READ(AS_BH, bhmdot, 5);
        check_block(&p, check);
    }
    g_mapped_file_unref(file);
}


/* this is a collective MPI routine
 * it reads all snapshot files and returns
 * the distributed particles in a PackedPar list
 * on current local processor.
 * */
PackedPar * snapshot_read(SnapHeader * h) {
    struct queue_element {
        PackedPar * pack;
        ptrdiff_t first;
    };

    ptrdiff_t Nlocal[7];
    ptrdiff_t Nglobal[6];
    int ptype;

    /* this must go first, it sets up Nreader, Nfile */
    int fidstart;
    int fidend;

    snapshot_prepare(h, Nglobal, &fidstart, &fidend);

    for(ptype = 0; ptype < 6; ptype++) {
        Nlocal[ptype] = (ReaderRank + 1) * Nglobal[ptype] / ReaderSize 
                  - ReaderRank * Nglobal[ptype] / ReaderSize;
    }
    Nlocal[6] = -1;
    /* This will be the particles stored locally, after the distribution.
     * */
    PackedPar * pack = NULL;

    LEADERONLY {
        int fid = fidstart;
        GQueue queue[6] = {
                G_QUEUE_INIT, 
                G_QUEUE_INIT, 
                G_QUEUE_INIT, 
                G_QUEUE_INIT, 
                G_QUEUE_INIT, 
                G_QUEUE_INIT
        };
        int rr = 0;
        for(rr = 0; rr < ReaderSize; rr++) {
            ptrdiff_t Nsend[7];
            for(ptype = 0; ptype < 6; ptype++) {
                Nsend[ptype] = (rr + 1) * Nglobal[ptype] / ReaderSize 
                          - rr * Nglobal[ptype] / ReaderSize;
            }
            Nsend[6] = -1;
            PackedPar * packsend = pstore_pack_create_a(Nsend);
            ptrdiff_t iter = 0;
            ptrdiff_t j[6] = {0};
            for(ptype = 0; ptype < 6; ptype ++) {
                struct queue_element * qe;
                while(j[ptype] < Nsend[ptype]) {
                    qe = g_queue_pop_head(&queue[ptype]);
                    while(qe == NULL) {
                        PackedPar * packread[6];
                        snapshot_to_pack(fid, packread);
                        fid ++;
                        int i;
                        for(i = 0; i < 6; i ++) {
                            if(packread[i]->size == 0) {
                                g_free(packread[i]);
                                continue;
                            }
                            struct queue_element * nqe = g_slice_new(struct queue_element);
                            nqe->pack = packread[i];
                            nqe->first = 0;
                            g_queue_push_tail(&queue[i], nqe);
                            g_debug("input queue of ptype %d: length = %d;\n",
                                i, 
                                g_queue_get_length(&queue[i]));
                        }
                        qe = g_queue_pop_head(&queue[ptype]);
                    }
                    for(;
                        j[ptype] < Nsend[ptype] && qe->first < qe->pack->size; 
                        j[ptype]++, qe->first++) {
                        Par * par = pstore_pack_get(qe->pack, qe->first); 
                        pstore_pack_push(packsend, &iter, par);
                    }
                    if(qe->first == qe->pack->size) {
                        g_free(qe->pack);
                        g_slice_free(struct queue_element, qe);
                        /* this chunk has been used up, need to read a new file */
                    } else {
                        g_queue_push_head(&queue[ptype], qe);
                    }
                }
            }
            g_assert(iter == packsend->size);
            /* now send packsend to the receiver */
            if(rr == 0) {
                pack = packsend;
            } else {
                g_debug("a pack is rippen and sent to %d", rr);
                MPI_Send(packsend, packsend->nbytes, MPI_BYTE, rr, rr, ReaderComm);
                g_free(packsend);
            }
        }
        for(ptype = 0; ptype < 6; ptype++) {
            g_assert(g_queue_is_empty(&queue[ptype]));
        }
        g_assert(fid == fidend);
    } else {
        pack = pstore_pack_create_a(Nlocal); 
        MPI_Recv(pack, pack->nbytes, MPI_BYTE, 
                    0, ReaderRank, ReaderComm, MPI_STATUS_IGNORE);
        g_debug("%03d received the particle pack", ThisTask);
    }
    return pack;
}

void pack_stat(PackedPar * pack, ptrdiff_t iter) {
        Par * par;
        ptrdiff_t i;
        uint64_t idmin = 0xfffffff, idmax = 0;
        for(i = 0; i < iter; i++) {
            Par * par = pstore_pack_get(pack, i);
            if(par->id > idmax) idmax = par->id;
            if(par->id < idmin) idmin = par->id;
        }
        g_debug("upto %td min %lu max %lu\n", iter, idmin, idmax);
}
