typedef struct _Domain {
    PSystem psys;
    TreeStore * treestore;
    int index;
    int prev;
    int next;
} Domain;

typedef struct {
    int HostTask;
    int Color; /* index on the HostTask */
} DomainTable;

void domain_init();
void domain_decompose();
void domain_build_tree();
void domain_destroy();

#define PAR(color, i) (D[color].psys.data[((signed)(i)<0)?((i)+D[color].psys.length):(i)])
#define NPAR(color)  D[color].psys.length

extern Domain * D;
extern int NColor;
