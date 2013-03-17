#ifndef __PAR_H__
#define __PAR_H__
#include <stdint.h>

typedef int32_t ipos_t;

typedef float real_t;

#define IPOS_NBITS 20
#define IPOS_LIMIT (1L << IPOS_NBITS)
#define IPOS_BAD (-1l)

/* this printf macro shall be disabled for non-debugging
 * builds.
 * it is non-portable (alloca), though seems to work fine
 * with icc. */
#define ipos_str(ipos) ({\
    char * result = alloca(IPOS_NBITS+1); \
    int i; \
    for(i = 0; i < IPOS_NBITS; i++) \
        result[i] = "01234567"[ipos_get_prefix(ipos, i)]; \
    result[i] = 0; \
    (char*) result; \
    })
#if 0
inline char * ipos_str(ipos_t ipos[3]) {
    char * result = malloc(IPOS_NBITS+1); 
    int i; 
    for(i = 0; i < IPOS_NBITS; i++) \
        result[i] = "01234567"[ipos_get_prefix(ipos, i)]; 
    result[i] = 0; 
    return result;
}
#endif
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
 * packed particles,
 * use g_free to free
 * ready for MPI send/recv.
 * */
typedef struct {
    size_t nbytes; /* total size in bytes, including the tailing index and this header */
    size_t size;   /* total number of particles */
    char data[];  
    /* data contains two parts, first a sequence of variable unitsize particles
     * then an index that points to the starting of each particles   */
} PackedPar;

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
inline void ipos_set_prefix(ipos_t ipos[3], int depth, int prefix) {
    int bit = (IPOS_NBITS - 1 - depth);
    ipos[0] = (ipos[0] & ~(1 << bit)) | (((prefix >> 0) & 1) << bit);
    ipos[1] = (ipos[1] & ~(1 << bit)) | (((prefix >> 1) & 1) << bit);
    ipos[2] = (ipos[2] & ~(1 << bit)) | (((prefix >> 2) & 1) << bit);
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
inline void ipos_bisect(ipos_t a[3], ipos_t b[3], ipos_t c[3]) {
    g_assert(ipos_compare(a, b) <= 0);
    int carry = 0;
    c[0] = 0; c[1] = 0; c[2] = 0;
    for(int depth = 0; depth < IPOS_NBITS; depth ++ ) {
        int prefix_a = ipos_get_prefix(a, depth);
        int prefix_b = ipos_get_prefix(b, depth);
        int temp = prefix_a + prefix_b + (carry << 3);
        int prefix_c = temp >> 1;
        carry = temp & 1;
        ipos_set_prefix(c, depth, prefix_c);
    }
}

inline void ipos_immediate_next(ipos_t a[3], ipos_t next[3]) {
    int carry = 1;
    for(int depth = IPOS_NBITS - 1; depth >= 0; depth -- ) {
        int prefix_a = ipos_get_prefix(a, depth);
        int temp = prefix_a + carry;
        carry = temp >> 3;
        ipos_set_prefix(next, depth, temp & 0x7);
    }

}

inline int par_is_primary(Par * par) {
    extern unsigned int PAR_PRIMARY_MASK;
    return (PAR_PRIMARY_MASK >> par->type) & 1;
}

void register_ptype(int ptype, char * name, size_t elesize, int is_primary);

/* these functions create /allocate a pstore */
PStore * pstore_new(size_t split_limit);
PStore * pstore_unpack(PackedPar * pack);

Par * pstore_insert_par(PStore * pstore, ipos_t ipos[3], int ptype);
void pstore_remove_par(PStore * pstore, Par * par);
void pstore_free_par(Par * par);
void pstore_remove_node(PStore * pstore, Node * node);
void pstore_free_node(Node * node);
Par * pstore_node_previous_par(Node * node);
Par * pstore_node_next_par(Node * node);
Par * pstore_get_nearby_par(PStore * pstore, ptrdiff_t index);


/* these function allocate/create a PackedPar structure */
PackedPar * pstore_pack_create_a(ptrdiff_t size[]);
PackedPar * pstore_pack_create_simple(int ptype, ptrdiff_t size, int filled);
PackedPar * pstore_pack_create_from_node(Node * node);
void pstore_pack_free(PackedPar * pack);

void pstore_pack_push(PackedPar * pack, ptrdiff_t * cursor, Par * par);
Par * pstore_pack_get(PackedPar * pack, ptrdiff_t cursor);
void pstore_pack_sort(PackedPar * pack);
ptrdiff_t pstore_pack_searchsorted_left(PackedPar * pack, ipos_t a[3]);
#endif
