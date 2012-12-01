
typedef uint64_t long_t;

#define float_t myfloat_t
typedef double float_t;

typedef float_t istate_t[3];
typedef float_t dstate_t[3];

typedef struct {
    fckey_t fckey; /* position */
    long_t id; /* for debugging */
    float_t mass;
    float_t rho;
    float_t T;
    float_t IGMmass;
    istate_t istate;
} par_t;

typedef struct {
    char name[8];
    par_t * base;
    par_t * data;
    size_t length;
    size_t size;
} PSystem;

void par_reserve(PSystem * psys, size_t size, size_t before);
par_t * par_append(PSystem * psys, intptr_t add);
par_t * par_prepend(PSystem * psys, intptr_t add);
void par_destroy(PSystem * psys);
void par_sort_by_fckey(PSystem * psys);
intptr_t par_search_by_fckey(PSystem * psys, fckey_t * key);

#define elsizeof(type, el) sizeof(((type *) 0)->el)

