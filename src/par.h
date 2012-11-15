
typedef uint64_t long_t;

typedef double float_t;

typedef float_t istate_t[3];
typedef float_t dstate_t[3];

typedef struct {
    fckey_t fckey; /* position */
    long_t id; /* for debugging */
    float_t mass; 
    float_t T;
    istate_t istate;
} par_t;

extern par_t * _PAR;
extern size_t NPAR;
extern size_t _PAR_size;
#define PAR(i) (_PAR[((signed)(i)<0)?((i)+NPAR):(i)])

extern par_t * _PARin;
extern size_t NPARin;
#define PARin(i) (_PARin[((signed)(i)<0)?((i)+NPARin):(i)])
size_t par_grow(size_t add);
void par_shift(size_t add);
void par_allocate(size_t size, size_t before);
void par_allocate_input(size_t size);
void par_free_input();
#define PAR_BUFFER_IN 0
#define PAR_BUFFER_MAIN 1
void par_sort_by_fckey(int which);
intptr_t par_search_by_fckey(fckey_t * key, int which);

#define elsizeof(type, el) sizeof(((type *) 0)->el)
