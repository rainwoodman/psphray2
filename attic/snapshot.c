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

static void snapshot_prepare(SnapHeader * h, ptrdiff_t Nglobal[6]) {
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
        int fid;
        for(fid = ReaderColor * Nfile / NReader;
            fid < (ReaderColor + 1) * Nfile / NReader;
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
    real_t rt;
    if(longid) {
        rt = * (uint64_t*) *p;
        *p += 8;
    } else {
        rt = * (uint32_t*) *p;
        *p += 4;
    }
    return rt;
}

static PackedPar * snapshot_to_pack(int fid, ptrdiff_t Nread[6]) {
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

    size_t Ntot = 0;
    ptrdiff_t N[7];
    int ptype;

    for(int i = 0; i < 6; i++) {
        Ntot += h.N[i];
        N[i] = h.N[i];
        Nread[i] = N[i];
    }
    N[6] = -1;

    PackedPar * temp = pstore_pack_create_a(NULL, N);
    ptrdiff_t iter = 0;
    for(ptype = 0; ptype < 6; ptype++) {
        Par * par = pstore_alloc_par(ptype);
        ptrdiff_t i;
        for(i = 0; i < N[ptype]; i++) {
            pstore_pack_push(temp, &iter, par);
        }
        pstore_free_par(par);
    }
    g_assert(iter == Ntot);

    p = buffer + 256 + 4 + 4; /* skip header */
    ptrdiff_t i = 0;
    int check;
    check = start_block(&p);
    for(i = 0; i < Ntot; i++) {
        Par * par = pstore_pack_get(temp, i);
        par->ipos[0] = read_real(&p, h.double_precision) / CB.BoxSize * IPOS_LIMIT;
        par->ipos[1] = read_real(&p, h.double_precision) / CB.BoxSize * IPOS_LIMIT;
        par->ipos[2] = read_real(&p, h.double_precision) / CB.BoxSize * IPOS_LIMIT;
    }
    check_block(&p, check);
    check = start_block(&p);
    for(i = 0; i < Ntot; i++) {
        Par * par = pstore_pack_get(temp, i);
        par->vel[0] = read_real(&p, h.double_precision);
        par->vel[1] = read_real(&p, h.double_precision);
        par->vel[2] = read_real(&p, h.double_precision);
    }
    check_block(&p, check);
    check = start_block(&p);
    for(i = 0; i < Ntot; i++) {
        Par * par = pstore_pack_get(temp, i);
        par->id = read_id(&p, CB.IDByteSize == 8);
    }
    check_block(&p, check);

    ptrdiff_t Nmass = 0;
    /* mass tab handling */
    for(int i=0; i< 6; i++) {
        if (h.masstab[i] != 0) Nmass += N[i];
    }

    if(Nmass > 0) {
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(h.masstab[par->type] == 0) {
                par->mass = read_real(&p, h.double_precision);
            } else {
                par->mass = h.masstab[par->type];
            }
        }
        check_block(&p, check);
    }
    if(N[0] > 0) {
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 0) continue;
            AS_GAS(par)->ie = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 0) continue;
            AS_GAS(par)->rho = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 0) continue;
            AS_GAS(par)->ye = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 0) continue;
            AS_GAS(par)->xHI = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 0) continue;
            AS_GAS(par)->sml = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 0) continue;
            AS_GAS(par)->sfr = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
    }
    if(N[0] + N[4] > 0) {
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 0) continue;
            AS_GAS(par)->met = read_real(&p, h.double_precision);
        }
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 4) continue;
            AS_STAR(par)->met = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
    }
    if(N[4] > 0) {
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 4) continue;
            AS_STAR(par)->sft = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
    }
    if(N[5] > 0) {
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 5) continue;
            AS_BH(par)->bhmass = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
        check = start_block(&p);
        for(i = 0; i < Ntot; i++) {
            Par * par = pstore_pack_get(temp, i);
            if(par->type != 5) continue;
            AS_BH(par)->bhmdot = read_real(&p, h.double_precision);
        }
        check_block(&p, check);
    }
    g_mapped_file_unref(file);
    return temp;
}
PackedPar * snapshot_read(SnapHeader * h) {


    ptrdiff_t Nlocal[7];
    ptrdiff_t Nsend[7];
    ptrdiff_t Nglobal[6];
    int ptype;

    /* this must go first, it sets up Nreader, Nfile */
    snapshot_prepare(h, Nglobal);

    for(ptype = 0; ptype < 6; ptype++) {
        Nlocal[ptype] = (ReaderRank + 1) * Nglobal[ptype] / ReaderSize 
                  - ReaderRank * Nglobal[ptype] / ReaderSize;
    }
    Nlocal[6] = -1;

    int fid;
    int fidstart = ReaderColor * Nfile / NReader;
    int fidend  = (ReaderColor + 1) * Nfile / NReader;

    PackedPar ** packread = g_newa(PackedPar*,  fidend - fidstart) - fidstart;
    ptrdiff_t (*Nremain) [6] = (ptrdiff_t (*)[6]) 
                g_newa(ptrdiff_t, (fidend - fidstart) * 6) - fidstart;
    int * stale = g_newa(int, fidend - fidstart) - fidstart;

    Par * head[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

    PackedPar * pack = pstore_pack_create_a(NULL, Nlocal); /* return value */

    LEADERONLY for(fid = fidstart; fid < fidend; fid ++) { 
        int ptype;
        ptrdiff_t cursor[6] = {0};
        int receiverrank = 0;
        int i;
        /* first read snapshot into a pack */
        packread[fid] = snapshot_to_pack(fid, Nremain[fid]);

        stale[fid] = 0;
        for(ptype = 0; ptype < 6; ptype ++) {
            if(Nremain[fid][ptype] == 0) stale[fid]++;
        }
        /* single linked list prepend and reverse = O(N),
         * append = O(N**2) */
        for(ptype = 0; ptype < 6; ptype ++) {
            head[ptype] = (Par * ) g_slist_reverse((GSList * ) head[ptype]);
        }
        for(i = 0; i < packread[fid]->size; i++) {
            /* convert pack to 6 link lists, one for each type, 
             * add to the active set */
            Par * p = pstore_pack_get(packread[fid], i);
            p->next = head[p->type];
            head[p->type] = p;
            cursor[p->type] ++;
        }
        for(ptype = 0; ptype < 6; ptype ++) {
            head[ptype] = (Par * ) g_slist_reverse((GSList * ) head[ptype]);
        }
        /* see if we are read to spawn a pack for a worker */
        int rippen;
        do {
            rippen = 0;
            for(ptype = 0; ptype < 6; ptype++) {
                Nsend[ptype] = (receiverrank + 1) * Nglobal[ptype] / ReaderSize 
                          - receiverrank * Nglobal[ptype] / ReaderSize;
            }
            Nsend[6] = -1;
            for(ptype = 0; ptype < 6; ptype ++) {
                if(cursor[ptype] >= Nsend[ptype]) {
                    rippen ++;
                }
            }
            /* not yet !*/
            if(rippen != 6) {
                break;
            }
            PackedPar * packsend = pstore_pack_create_a(NULL, Nsend);
            ptrdiff_t iter = 0; 
            for(ptype = 0; ptype < 6; ptype ++) {
                Par * par;
                ptrdiff_t j; 
                /* push this many each type to the spawning pack */
                for(par = head[ptype], j = 0; j < Nsend[ptype]; j++, par = par->next) {
                    pstore_pack_push(packsend, &iter, par);
                }
                /* and remove them from the active set */
                head[ptype] = par;
                cursor[ptype] -= Nsend[ptype];
                /* see if any of the packread is no longer in use and free it */
                ptrdiff_t deduct = Nsend[ptype];
                int fid2;
                for(fid2 = fidstart; fid2 <= fid; fid2++) {
                    ptrdiff_t size = Nremain[fid2][ptype];
                    if(size >= deduct) {
                        size -= deduct;
                        deduct = 0;
                    } else {
                        deduct -= size;
                        size = 0;
                    }
                    if(size == 0 && Nremain[fid2][ptype] != 0 
                        && stale[fid2] < 6) {
                        /* toggle the staleness if Nremain drops to zero*/
                        stale[fid2] ++;
                    }
                    Nremain[fid2][ptype] = size;
                    /* if any particle type is nolonger in use, then
                     * staleness goes up. when staleness goes to 6 then it's done */
                    if(stale[fid2] == 6) {
                        g_free(packread[fid2]);  
                        packread[fid2] = NULL;
                    }
                }
            }
            /* now send packsend to the receiver */
            if(receiverrank == 0) {
                pack = packsend;
            } else {
                MPI_Send(packsend, pstore_pack_total_bytes(packsend), MPI_BYTE, receiverrank, receiverrank, ReaderComm);
                g_free(packsend);
            }
            receiverrank ++;
        } while(TRUE);
    }
    for(fid = fidstart; fid < fidend; fid++) {
        g_assert(packread[fid] == NULL);
    }
    LEADERONLY {} else {
        MPI_Recv(pack, pstore_pack_total_bytes(pack), MPI_BYTE, 
                    0, ReaderRank, ReaderComm, MPI_STATUS_IGNORE);
    }
    return pack;
}
