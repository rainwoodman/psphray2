typedef union {
    uint64_t a[2];
    uint8_t c[16];
} fckey_t;

#define FCKEY_BITS 16
#define FCKEY_MAX ((1L << FCKEY_BITS) - 1)

#define FCKEY_FMT "%lX%0.16lX"
#define FCKEY_PRINT(x) x.a[1], x.a[0]

int fckey_cmp(fckey_t * key1, fckey_t * key2);
void fckey_clear(fckey_t * key, int bits);
void fckey_to_ipos(fckey_t * key, int64_t pos[3]);
void fckey_from_ipos(fckey_t * key, int64_t pos[3]);
void fckey_center(fckey_t * center, fckey_t * key1, fckey_t * key2);
void fckey_set_max(fckey_t * key);

