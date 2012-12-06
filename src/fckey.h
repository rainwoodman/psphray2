typedef union {
    uint64_t a[2];
    uint8_t c[16];
} fckey_t;

extern unsigned int FCKEY_BITS;
#define FCKEY_MAX ((uint64_t) ((((uint64_t)1L) << FCKEY_BITS) - 1L))

#define FCKEY_FMT "%0.*lo%0.*lo%0.*lo"
#define FCKEY_PRINT_PREFIX(x, n) \
_fckey_print_width(x, 2, n), _fckey_print_value(x, 2, FCKEY_BITS - (n)), \
_fckey_print_width(x, 1, n), _fckey_print_value(x, 1, FCKEY_BITS - (n)), \
_fckey_print_width(x, 0, n), _fckey_print_value(x, 0, FCKEY_BITS - (n))
#define FCKEY_PRINT(x) FCKEY_PRINT_PREFIX(x, FCKEY_BITS)

uint64_t _fckey_print_value(fckey_t key, int part, int rightshiftdigits);
int _fckey_print_width(fckey_t key, int part, int digits);

int fckey_cmp(fckey_t * key1, fckey_t * key2);
void fckey_clear(fckey_t * key, int bits);
void fckey_set(fckey_t * key, int bits);
void fckey_to_ipos(fckey_t * key, int64_t pos[3]);
void fckey_from_ipos(fckey_t * key, int64_t pos[3]);
void fckey_center(fckey_t * center, fckey_t * key1, fckey_t * key2);
void fckey_set_max(fckey_t * key);
int fckey_is_zero(fckey_t * key);
void fckey_xor(fckey_t * result, fckey_t * key1, fckey_t * key2);
void fckey_rightshift(fckey_t * key, int offset);
void fckey_or_with_leftshift(fckey_t * key, uint32_t value, int t);
