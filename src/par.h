typedef struct { 
    union {
        uint64_t x;
        uint64_t y;
        uint64_t z;
    };
    uint64_t a[3];
} fckey_t;

typedef uint64_t long_t;
typedef double float_t;
typedef float_t istate_t[3];
typedef float_t dstate_t[3];

typedef struct {
    union {
        fckey_t fckey; /* position */
        float_t pos[3];
    };
    long_t id; /* for debugging */
    float_t mass; 
    float_t T;
    istate_t istate;
} par_t;

extern GArray * _PAR;
#define PAR(i) (g_array_index(_PAR, par_t, ((i)<0)?((i)+_PAR->len):(i)))
#define PAREND (g_array_index(_PAR, par_t, (_PAR->len)))
void par_increase_size(size_t add);
void par_preallocate(size_t size);
#define elsizeof(type, el) sizeof(((type *) 0)->el)
#define NPAR _PAR->len
