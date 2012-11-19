#define NODE_TYPE_INNER 'N'
#define NODE_TYPE_LEAF 'L'
struct _InnerNode;
struct _Node;

#define NODE_COMMON \
    char type; \
    char order; \
    char complete; \
    char nchildren; \
    int32_t npar;  \
    par_t * first; \
    fckey_t key; \
    struct _InnerNode * parent; \
    struct _Node * link; 
/* link is used in to free the tree */

typedef struct _Node {
    NODE_COMMON;
} Node;

typedef struct _InnerNode {
    NODE_COMMON;
    Node * child[8];
} InnerNode;

#define NODE_FMT \
    "{ " \
    "%c%s " \
    FCKEY_FMT "/"\
    "%d " \
    "fst: " "%d " \
    "np: " "%d " \
    "nc: " "%d " \
    "}"

#define NODE_PRINT(x) \
    (x).type, \
    (x).complete?"C":"I", \
    FCKEY_PRINT((x).key), \
    (x).order, \
    (int)((x).first - &PAR(0)), \
    (x).npar, \
    (x).nchildren

typedef struct {
  Node * stack[128];
  int8_t c[128];
  Node * root;
  int top;
} TreeIter;

void tree_build();
#define tree_free() (tree_destroy(TREEROOT), TREEROOT = NULL)
void tree_destroy(Node * root);
void tree_link(Node * root, Node ** firstleaf, InnerNode ** firstinner);
Node * tree_locate_fckey(Node * root, fckey_t * key);
int tree_node_contains_fckey(Node * node, fckey_t * key);
int tree_node_contains_node(Node * node, Node * needle);
Node * tree_node_find_image(Node * node, Node * needle);

TreeIter * tree_iter_new(Node * root);
Node * tree_iter_next(TreeIter * iter);
Node * tree_iter_next_sibling(TreeIter * iter);
void tree_iter_free(TreeIter * iter);

extern Node * TREEROOT;

