typedef struct {
    PSystem psys;
    TreeStore treestore;
    Node * tree;
    int HostTask;
} Domain;

void domain_decompose();
extern Domain D;
#define PAR(i) (D.psys.data[((signed)(i)<0)?((i)+D.psys.length):(i)])
#define NPAR  D.psys.length
