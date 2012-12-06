
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

typedef struct _PSystem PSystem;

typedef struct {
    PSystem * psys;
    intptr_t i;
    intptr_t end;
    intptr_t last_i;
} ParIter;

PSystem * par_alloc();
void par_free(PSystem * psys);
void par_init(PSystem * psys, char * name);
void par_reserve(PSystem * psys, size_t size, size_t before);
par_t * par_append(PSystem * psys, intptr_t add);
par_t * par_prepend(PSystem * psys, intptr_t add);

void par_destroy(PSystem * psys);
void par_sort_by_fckey(PSystem * psys);
intptr_t par_search_by_fckey(PSystem * psys, fckey_t * key);
size_t par_get_length(PSystem * psys);

par_t * par_index(PSystem * psys, intptr_t index);

par_t * par_iter_init(ParIter * iter, PSystem * psys);
par_t * par_iter_init_range(ParIter * iter, PSystem * psys, intptr_t first, size_t npar);
par_t * par_iter_set(ParIter * iter, intptr_t index);
intptr_t par_iter_last_index(ParIter * iter);
inline par_t * par_iter_next(ParIter * iter);
#define elsizeof(type, el) sizeof(((type *) 0)->el)

