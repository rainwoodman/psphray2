typedef struct {
    PSystem psys;
    TreeStore treestore;
    Node * tree;
    int HostTask;
} Domain;

void domain_decompose(Domain * domain);
void domain_build_tree(Domain * domain);
void domain_destroy(Domain * domain);
#define PAR(d, i) (d->psys.data[((signed)(i)<0)?((i)+d->psys.length):(i)])
#define NPAR(d)  d->psys.length
