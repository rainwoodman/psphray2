#include <stdint.h>

typedef int32_t ipos_t;
#define IPOS_NBITS 20
#define IPOS_LIMIT (1L << IPOS_NBITS)
#define IPOS_BAD (-1l)

typedef struct _Par {
    union {
        char type;
        gpointer space;
    };
    struct _Par * next;
    ipos_t ipos[3];
    uint64_t id;
    float mass;
    char data[];
} Par;

struct _ParNode {
    struct _ParNode * link[8];
    struct _ParNode * up;
    Par * first;
    Par * last;
    size_t primary_size; /* total number of primary particles, for split */
    size_t size; /* total number of particles, for mem/lookup */
    void * data; /* extra data associated with a node */
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
    int primarymask;  /* the mask of the primary ptype */
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
inline int par_is_primary(PStore * pstore, Par * par) {
    return ((1 << par->type) & pstore->primarymask) != 0;
}

int register_ptype(int id, char * name, size_t elesize);
PStore * pstore_new(int primarymask);

Par * pstore_insert(PStore * pstore, ipos_t ipos[3], int ptype);
void pstore_remove(PStore * pstore, Par * par);
Par * pstore_node_previous_par(Node * node);
Par * pstore_node_next_par(Node * node);
Par * pstore_get(PStore * pstore, ptrdiff_t index);

/**
 * packed particles,
 * use g_free to free
 * ready for MPI send/recv.
 * */
typedef struct {
    size_t size;   /* size of data in number of particles */
    size_t bytes;  /* size of data in bytes */
    char data[];
} PackedPar;

PackedPar * pstore_pack(Par * first, size_t size);
Par * pstore_unpack(PackedPar * pack);
