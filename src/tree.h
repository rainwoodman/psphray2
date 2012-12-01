#define NODE_TYPE_INNER 'N'
#define NODE_TYPE_LEAF 'L'
struct _InnerNode;
struct _Node;

#define NODE_COMMON \
    char type; \
    char order; \
    char complete; \
    char unused; \
    int32_t npar;  \
    int32_t ifirst; \
    fckey_t key; \
    union { \
        struct _InnerNode * parent; \
        intptr_t iparent; \
    };

typedef struct _Node {
    NODE_COMMON;
} Node;

typedef struct _InnerNode {
    NODE_COMMON;
    union {
        Node * child[8];
        intptr_t ichild[8];
    };
} InnerNode;

typedef struct {
    Stack inner;
    Stack leaf;
    PSystem * psys;
} TreeStore;

#define NODE_FMT \
    "{ " \
    "%c%s " \
    FCKEY_FMT "/"\
    "%d " \
    "fst: " "%d " \
    "np: " "%d " \
    "nc: " "%d " \
    "}"
static inline int tree_node_nchildren(Node * node) {
    if(node->type == NODE_TYPE_INNER) {
        int i = 0, p = 0;
        for(i = 0; i < 8; i++) {
            p += (((InnerNode*) node)->child[i] != NULL);
        }
        return p;
    } else {
        return 0;
    }
}
#define NODE_PRINT(x) \
    (x).type, \
    (x).complete?"C":"I", \
    FCKEY_PRINT_PREFIX(x.key, FCKEY_BITS - (x).order), \
    (x).order, \
    (x).ifirst, \
    (x).npar, \
    tree_node_nchildren((Node*)&(x))

typedef struct {
  Node * root;
  Node * current;
} TreeIter;

Node * tree_build(TreeStore * store);

void tree_store_init(TreeStore * store, PSystem * psys);
void tree_store_destroy(TreeStore * store);

Node * tree_locate_fckey(Node * root, fckey_t * key);
int tree_node_contains_fckey(Node * node, fckey_t * key);
int tree_node_contains_node(Node * node, Node * needle);
Node * tree_node_find_image(Node * node, Node * needle);

void tree_iter_init(TreeIter * iter, Node * root);
Node * tree_iter_next(TreeIter * iter);
Node * tree_iter_next_sibling(TreeIter * iter);

