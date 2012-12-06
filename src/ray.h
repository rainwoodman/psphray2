typedef struct {
    uint64_t id;
    double pos[3];
    double dir[3];
} ray_t;
typedef struct {
    int domainindex;
    int rayindex;
    intptr_t nodeindex;
    double tE;
    double tL;
} intersect_t;

#define RAY_FMT "%ld (%g,%g,%g)->(%g,%g,%g)"
#define RAY_PRINT(x) x.id, x.pos[0], x.pos[1], x.pos[2], x.dir[0], x.dir[1], x.dir[2]
extern ray_t * R;
extern int NRay;

