#include <stdlib.h>
#include <stdint.h>
#include "fckey.h"
unsigned int FCKEY_BITS = 40;

uint64_t _fckey_print_value(fckey_t key, int part, int rightshiftdigits) {
    fckey_rightshift(&key, rightshiftdigits * 3);
    switch(part) {
        case 0:
            return key.a[0] & 077777777777; /* 33 bits */
        case 1:
            return ((key.a[0] >> 33) & 077777777777) + ((key.a[1] & 3) << 31); /* 31 + 2 = 33 bits */
        case 2:
            return (key.a[1] >> 2);
    }
    abort();
}

int _fckey_print_width(fckey_t key, int part, int digits) {
    switch(part) {
        case 2:
            if (digits > 22) 
                return digits - 22;
            else
                return 0;
        case 1:
            if(digits > 11) {
                if(digits < 22)
                    return digits - 11;
                else
                    return 10;
            } else
                return 0;
        case 0:
            if(digits < 11)
                return digits;
            else
                return 11;
    }
    abort();
}


/* t is in bits */
void fckey_or_with_leftshift(fckey_t * key, uint32_t value, int t) {
    if(t < 64) {
    /* low is affected */
       key->a[0] |= (uint64_t)value << t;
    }
    if(t > 32) {
    /* high is affected */
        t -= 64;
        if(t < 0) {
           key->a[1] |= (uint64_t)value >> (-t);
        } else {
           key->a[1] |= (uint64_t)value << t;
        }
    }
}

void fckey_rightshift(fckey_t * key, int offset) {
    if(offset == 0) return;
    if(offset < 64) {
        key->a[0] >>= offset;
        key->a[0] |= (key->a[1] & ((1L << offset) - 1)) << (64 - offset);
        key->a[1] >>= offset;
    } else {
        offset -= 64;
        key->a[0] = key->a[1] >> offset;
        key->a[1] = 0;
    }
}

static uint32_t utab[];
static inline void fckey_fill(fckey_t * key, int64_t x, int offset) {
    int i;
    for (i = 0; i < FCKEY_BITS; i+= 8) {
        fckey_or_with_leftshift(key, utab[(uint8_t) x], (i * 3 + offset));
        x >>= 8;
    }
}
static uint8_t ctab[];
static inline void fckey_extract(const fckey_t * key, int64_t *x, int offset) {
    int i, b;
    fckey_t temp;
    uint8_t raw = 0;
    static fckey_t mask = {
        {0x9249249249249249, 0x4924924924924924},
    };
    temp = *key;
    fckey_rightshift(&temp, offset);
    temp.a[0] = temp.a[0] & mask.a[0];
    temp.a[1] = temp.a[1] & mask.a[1];
    *x = 0;
    for(b = 0, i=0; b < FCKEY_BITS; b+=8, i+=3) {
        raw = temp.c[i];
        raw |= temp.c[i+1];
        raw |= temp.c[i+2];
        *x |= (int64_t )ctab[raw] << b;
    }
}

void fckey_clear(fckey_t * key, int bits) {
    if(bits >= 64) {
        key->a[0] = 0;
        bits -= 64;
        key->a[1] = (key->a[1] >> bits) << bits;
    } else {
        key->a[0] = (key->a[0] >> bits) << bits;
    }
}
void fckey_set(fckey_t * key, int bits) {
    if(bits >= 64) {
        key->a[0] = -1L;
        bits -= 64;
        key->a[1] = (key->a[1] >> bits) << bits;
        key->a[1] |= (1L << bits) - 1;
    } else {
        key->a[0] = (key->a[0] >> bits) << bits;
        key->a[0] |= (1L << bits) - 1;
    }
}

static int uint64_compare(uint64_t * p1, uint64_t * p2) {
    return (*p1 > *p2) - (*p2 > *p1);
}
static uint64_t uint64_center(uint64_t p1, uint64_t p2) {
    return (p1 >> 1) + (p2 >> 1) + ((p1 & 1) && (p2 & 1));
}
void fckey_center(fckey_t * center, fckey_t * key1, fckey_t * key2) {
    /* this is quite rough if the high bit is */
    if(key1->a[1] == key2->a[1]) {
        center->a[1] = key1->a[1];
        center->a[0] = uint64_center(key1->a[0], key2->a[0]);
    } else {
        center->a[1] = uint64_center(key1->a[1], key2->a[1]);
        center->a[0] = uint64_center(key1->a[0], key2->a[0]);
        if((key1->a[1] & 1) != (key2->a[1] & 1)) {
            center->a[0] += (1uL<<63) - 1;
        }
    }
}
int fckey_cmp(fckey_t * key1, fckey_t * key2) {
    int r1 = uint64_compare(&key1->a[1], &key2->a[1]);
    if(r1 == 0) {
        return uint64_compare(&key1->a[0], &key2->a[0]);
    } else {
        return r1;
    }
}
void fckey_to_ipos(fckey_t * key, int64_t pos[3]) {
    int d;
    for(d=0; d < 3; d++) 
        fckey_extract(key, &pos[d], d);
}
void fckey_from_ipos(fckey_t * key, int64_t pos[3]) {
    int d;
    key->a[0] = 0;
    key->a[1] = 0;
    for(d=0; d < 3; d++)
        fckey_fill(key, (uint64_t) pos[d] & FCKEY_MAX, d);
}

void fckey_set_zero(fckey_t * key) {
    key->a[0] = 0;
    key->a[1] = 0;
}
void fckey_minus_one(fckey_t * key) {
    if(key->a[0] == 0) {
        key->a[0] = -1;
        key->a[1] --;
    } else {
        key->a[0] --;
    }
}
void fckey_set_max(fckey_t * key) {
    static int64_t pos[3];
    pos[0] = FCKEY_MAX;
    pos[1] = FCKEY_MAX;
    pos[2] = FCKEY_MAX;
    fckey_from_ipos(key, pos);
}
void fckey_xor(fckey_t * result, fckey_t * key1, fckey_t * key2) {
    result->a[0] = key1->a[0] ^ key2->a[0];
    result->a[1] = key1->a[1] ^ key2->a[1];
}
int fckey_is_zero(fckey_t * key) {
    return (key->a[0] == 0) && (key->a[1] == 0);
}
#if 0
#define FMT8 "%0.16lX "
#define FMT16 FMT8 FMT8
void test(int bits) {
    fckey_t key = {0, 0};
    fckey_fill(&key, 0xdeadbeefdeadbeef, bits);
    int64_t x;
    fckey_extract(&key, &x, bits);
    printf( "%d " FMT16 " %0.16lX \n", bits, key.a[1], key.a[0], x);
}
void test2(int bits) {
    fckey_t key1 = {0, 0}, key2 = {0, 0}, center;
    fckey_or_with_leftshift(&key1, 0xdeadbeefdeadbeef, bits);
    fckey_or_with_leftshift(&key2, 0xdeadbeefdeadbeef - 1, bits);
    fckey_center(&center, &key1, &key2);
    printf( "%d " FMT16  FMT16 "\n", bits, 
            key1.a[1], key1.a[0], center.a[1], center.a[0]);
}
void test3() {
    fckey_t key1;
    fckey_set_max(&key1);
    printf( FCKEY_FMT "\n", FCKEY_PRINT(key1));
    printf( "%lo %lo\n", key1.a[0], key1.a[1]);
    key1.a[0] = 0345670123456701234567L;
    key1.a[1] = 01234567012345670123456712L << 2;
    printf( FCKEY_FMT "\n", FCKEY_PRINT(key1));
}
void main() {
    int i = 0;
    test3(i);
}
#endif

static uint32_t utab[256] = {
    0x000000, 0x000001, 0x000008, 0x000009,
    0x000040, 0x000041, 0x000048, 0x000049,
    0x000200, 0x000201, 0x000208, 0x000209,
    0x000240, 0x000241, 0x000248, 0x000249,
    0x001000, 0x001001, 0x001008, 0x001009,
    0x001040, 0x001041, 0x001048, 0x001049,
    0x001200, 0x001201, 0x001208, 0x001209,
    0x001240, 0x001241, 0x001248, 0x001249,
    0x008000, 0x008001, 0x008008, 0x008009,
    0x008040, 0x008041, 0x008048, 0x008049,
    0x008200, 0x008201, 0x008208, 0x008209,
    0x008240, 0x008241, 0x008248, 0x008249,
    0x009000, 0x009001, 0x009008, 0x009009,
    0x009040, 0x009041, 0x009048, 0x009049,
    0x009200, 0x009201, 0x009208, 0x009209,
    0x009240, 0x009241, 0x009248, 0x009249,
    0x040000, 0x040001, 0x040008, 0x040009,
    0x040040, 0x040041, 0x040048, 0x040049,
    0x040200, 0x040201, 0x040208, 0x040209,
    0x040240, 0x040241, 0x040248, 0x040249,
    0x041000, 0x041001, 0x041008, 0x041009,
    0x041040, 0x041041, 0x041048, 0x041049,
    0x041200, 0x041201, 0x041208, 0x041209,
    0x041240, 0x041241, 0x041248, 0x041249,
    0x048000, 0x048001, 0x048008, 0x048009,
    0x048040, 0x048041, 0x048048, 0x048049,
    0x048200, 0x048201, 0x048208, 0x048209,
    0x048240, 0x048241, 0x048248, 0x048249,
    0x049000, 0x049001, 0x049008, 0x049009,
    0x049040, 0x049041, 0x049048, 0x049049,
    0x049200, 0x049201, 0x049208, 0x049209,
    0x049240, 0x049241, 0x049248, 0x049249,
    0x200000, 0x200001, 0x200008, 0x200009,
    0x200040, 0x200041, 0x200048, 0x200049,
    0x200200, 0x200201, 0x200208, 0x200209,
    0x200240, 0x200241, 0x200248, 0x200249,
    0x201000, 0x201001, 0x201008, 0x201009,
    0x201040, 0x201041, 0x201048, 0x201049,
    0x201200, 0x201201, 0x201208, 0x201209,
    0x201240, 0x201241, 0x201248, 0x201249,
    0x208000, 0x208001, 0x208008, 0x208009,
    0x208040, 0x208041, 0x208048, 0x208049,
    0x208200, 0x208201, 0x208208, 0x208209,
    0x208240, 0x208241, 0x208248, 0x208249,
    0x209000, 0x209001, 0x209008, 0x209009,
    0x209040, 0x209041, 0x209048, 0x209049,
    0x209200, 0x209201, 0x209208, 0x209209,
    0x209240, 0x209241, 0x209248, 0x209249,
    0x240000, 0x240001, 0x240008, 0x240009,
    0x240040, 0x240041, 0x240048, 0x240049,
    0x240200, 0x240201, 0x240208, 0x240209,
    0x240240, 0x240241, 0x240248, 0x240249,
    0x241000, 0x241001, 0x241008, 0x241009,
    0x241040, 0x241041, 0x241048, 0x241049,
    0x241200, 0x241201, 0x241208, 0x241209,
    0x241240, 0x241241, 0x241248, 0x241249,
    0x248000, 0x248001, 0x248008, 0x248009,
    0x248040, 0x248041, 0x248048, 0x248049,
    0x248200, 0x248201, 0x248208, 0x248209,
    0x248240, 0x248241, 0x248248, 0x248249,
    0x249000, 0x249001, 0x249008, 0x249009,
    0x249040, 0x249041, 0x249048, 0x249049,
    0x249200, 0x249201, 0x249208, 0x249209,
    0x249240, 0x249241, 0x249248, 0x249249,
};
static uint8_t ctab[256] = {
0x00,  0x01,  0x08,  0x09, 
0x40,  0x41,  0x48,  0x49, 
0x02,  0x03,  0x0a,  0x0b, 
0x42,  0x43,  0x4a,  0x4b, 
0x10,  0x11,  0x18,  0x19, 
0x50,  0x51,  0x58,  0x59, 
0x12,  0x13,  0x1a,  0x1b, 
0x52,  0x53,  0x5a,  0x5b, 
0x80,  0x81,  0x88,  0x89, 
0xc0,  0xc1,  0xc8,  0xc9, 
0x82,  0x83,  0x8a,  0x8b, 
0xc2,  0xc3,  0xca,  0xcb, 
0x90,  0x91,  0x98,  0x99, 
0xd0,  0xd1,  0xd8,  0xd9, 
0x92,  0x93,  0x9a,  0x9b, 
0xd2,  0xd3,  0xda,  0xdb, 
0x04,  0x05,  0x0c,  0x0d, 
0x44,  0x45,  0x4c,  0x4d, 
0x06,  0x07,  0x0e,  0x0f, 
0x46,  0x47,  0x4e,  0x4f, 
0x14,  0x15,  0x1c,  0x1d, 
0x54,  0x55,  0x5c,  0x5d, 
0x16,  0x17,  0x1e,  0x1f, 
0x56,  0x57,  0x5e,  0x5f, 
0x84,  0x85,  0x8c,  0x8d, 
0xc4,  0xc5,  0xcc,  0xcd, 
0x86,  0x87,  0x8e,  0x8f, 
0xc6,  0xc7,  0xce,  0xcf, 
0x94,  0x95,  0x9c,  0x9d, 
0xd4,  0xd5,  0xdc,  0xdd, 
0x96,  0x97,  0x9e,  0x9f, 
0xd6,  0xd7,  0xde,  0xdf, 
0x20,  0x21,  0x28,  0x29, 
0x60,  0x61,  0x68,  0x69, 
0x22,  0x23,  0x2a,  0x2b, 
0x62,  0x63,  0x6a,  0x6b, 
0x30,  0x31,  0x38,  0x39, 
0x70,  0x71,  0x78,  0x79, 
0x32,  0x33,  0x3a,  0x3b, 
0x72,  0x73,  0x7a,  0x7b, 
0xa0,  0xa1,  0xa8,  0xa9, 
0xe0,  0xe1,  0xe8,  0xe9, 
0xa2,  0xa3,  0xaa,  0xab, 
0xe2,  0xe3,  0xea,  0xeb, 
0xb0,  0xb1,  0xb8,  0xb9, 
0xf0,  0xf1,  0xf8,  0xf9, 
0xb2,  0xb3,  0xba,  0xbb, 
0xf2,  0xf3,  0xfa,  0xfb, 
0x24,  0x25,  0x2c,  0x2d, 
0x64,  0x65,  0x6c,  0x6d, 
0x26,  0x27,  0x2e,  0x2f, 
0x66,  0x67,  0x6e,  0x6f, 
0x34,  0x35,  0x3c,  0x3d, 
0x74,  0x75,  0x7c,  0x7d, 
0x36,  0x37,  0x3e,  0x3f, 
0x76,  0x77,  0x7e,  0x7f, 
0xa4,  0xa5,  0xac,  0xad, 
0xe4,  0xe5,  0xec,  0xed, 
0xa6,  0xa7,  0xae,  0xaf, 
0xe6,  0xe7,  0xee,  0xef, 
0xb4,  0xb5,  0xbc,  0xbd, 
0xf4,  0xf5,  0xfc,  0xfd, 
0xb6,  0xb7,  0xbe,  0xbf, 
0xf6,  0xf7,  0xfe,  0xff, 
};
