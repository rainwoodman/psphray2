#ifndef __PAR_H__
#define __PAR_H__
#include <stdint.h>

typedef int32_t ipos_t;
typedef float real_t;

#define IPOS_NBITS 20
#define IPOS_LIMIT (1L << IPOS_NBITS)
#define IPOS_BAD (-1l)

typedef struct _Par {
    union {
        char type;
        gpointer space;
    };
    struct _Par * next;
    /* do not modify above. for GSList compatibility. */
    ipos_t ipos[3];
    real_t vel[3];
    uint64_t id;
    real_t mass;
    char prop[];
} Par;

struct _ParNode {
    struct _ParNode * link[8];
    struct _ParNode * up;
    Par * first;
    Par * last;
    size_t primary_size; /* total number of primary particles, for split */
    size_t size; /* total number of particles, for mem/lookup */
    void * data; /* extra data associated with a node */
    char first_nonempty_child;
    char last_nonempty_child;
    int prefix; /* duplicated, this is also the position in up->link. avoid a linear search */
    int fill;
};
#define PARNODE_FMT "up %p par %p size %d"
#define PARNODE_PRINT(x) x.up, x.par, x.size
#define PAR_FMT "ipos (%g %g %g) type %d"
#define PAR_PRINT(x) \
        x.ipos[0] * 1.0 / IPOS_LIMIT, \
        x.ipos[1] * 1.0 / IPOS_LIMIT, \
        x.ipos[2] * 1.0 / IPOS_LIMIT, \
        x.type

typedef struct _ParNode Node;

typedef struct {
    int split_limit; /* when to split a node, this does not decide the final shape of the octtree */
    int merge_limit; /* when to merge all child nodes. this is the lower limit of # particles (unused) */
    int depth_limit; /* will never split the tree to deeper than this, shall usually be IPOS_NBITS */
    struct _ParNode * root;
} PStore;

/**
 * returns the prefix of the ipos at given depth.
 *
 * depth starts from 0, which is the highest bit
 */
inline int ipos_get_prefix(ipos_t ipos[3], int depth) {
    int bit = (IPOS_NBITS - 1 - depth);
    return 
      (((ipos[0] >> bit) & 1) << 0) |
      (((ipos[1] >> bit) & 1) << 1) |
      (((ipos[2] >> bit) & 1) << 2) 
        ;
}
inline int ipos_compare(ipos_t a[3], ipos_t b[3]) {
    for(int depth = 0; depth < IPOS_NBITS; depth ++ ) {
        int prefix_a = ipos_get_prefix(a, depth);
        int prefix_b = ipos_get_prefix(b, depth);
        if (prefix_a < prefix_b) {
            return -1;
        } 
        if (prefix_a > prefix_b) {
            return 1;
        }
    }
    return 0;
}
inline int par_is_primary(Par * par) {
    extern unsigned int PAR_PRIMARY_MASK;
    return (PAR_PRIMARY_MASK >> par->type) & 1;
}

void register_ptype(int ptype, char * name, size_t elesize, int is_primary);
PStore * pstore_new(size_t split_limit);

void pstore_free_node(Node * node);
void pstore_free_node_r(Node * node);
Par * pstore_alloc_par(int ptype);
void pstore_free_par(Par * par);
void pstore_free_par_chain(Par * head);
Par * pstore_insert(PStore * pstore, ipos_t ipos[3], int ptype);
void pstore_remove(PStore * pstore, Par * par);
Par * pstore_node_previous_par(Node * node);
Par * pstore_node_next_par(Node * node);
Par * pstore_get_nearby(PStore * pstore, ptrdiff_t index);

/**
 * packed particles,
 * use g_free to free
 * ready for MPI send/recv.
 * */
typedef struct {
    size_t totalbytes; /* total size in bytes, including the tailing index and this header */
    size_t size;   /* total number of particles */
    char data[];  
    /* data contains two parts, first a sequence of variable unitsize particles
     * then an index that points to the starting of each particles   */
} PackedPar;

PackedPar * pstore_pack_create_a(ptrdiff_t size[]);
PackedPar * pstore_pack_create_simple(int ptype, ptrdiff_t size);
PackedPar * pstore_pack(Par * first, size_t size);

void pstore_pack_push(PackedPar * pack, ptrdiff_t * cursor, Par * par);
Par * pstore_pack_get(PackedPar * pack, ptrdiff_t cursor);
void pstore_pack_sort(PackedPar * pack);

#endif
