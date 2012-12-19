#define NODE_TYPE_INNER '>'
#define NODE_TYPE_LEAF '+'

#define NODE_FREE 'X'
#define NODE_LOCAL 'L'
#define NODE_INCOMPLETE 'I'
#define NODE_REMOTE 'R'

/* GHOST is a special leaf node, it stores which domain
 * is intersecting this node from ifirst to npar */

struct _InnerNode;
struct _Node;

#define NODE_COMMON \
    char type; \
    char order; \
    char status; \
    char unused; \
    int32_t iparent; \
    int32_t npar;  \
    int32_t ifirst; \
    fckey_t key; \

typedef struct _Node {
    NODE_COMMON;
} Node;

typedef struct _InnerNode {
    NODE_COMMON;
    int32_t ichild[8];
} InnerNode;

/* hide implementation details */
typedef struct _TreeStore TreeStore;

#define NODE_FMT \
    "{ " \
    "%c%c " \
    FCKEY_FMT "/" \
    "%d " \
    "fst: " "%d " \
    "np: " "%d " \
    "nc: " "%d " \
    "}"
static inline int tree_node_nchildren(Node * node) {
    if(node->type == NODE_TYPE_INNER) {
        int i = 0, p = 0;
        for(i = 0; i < 8; i++) {
            p += (((InnerNode*) node)->ichild[i] != 0);
        }
        return p;
    } else {
        return 0;
    }
}
#define NODE_PRINT(x) \
    (x).type, \
    (x).status?(x).status:'?', \
    FCKEY_PRINT_PREFIX(x.key, FCKEY_BITS - (x).order), \
    (x).order, \
    (x).ifirst, \
    (x).npar, \
    tree_node_nchildren((Node*)&(x))

typedef struct {
  TreeStore * store;
  Node * root;
  Node * current;
} TreeIter;

TreeStore * tree_store_alloc();
void tree_store_free(TreeStore * store);
void tree_store_init(TreeStore * store, PSystem * psys, int splitthresh);
void tree_store_destroy(TreeStore * store);
Node * tree_store_get_node(TreeStore * store, intptr_t i);
intptr_t tree_store_get_index(TreeStore * store, Node * node);

#define tree_store_root(store) tree_store_get_node(store, 1)

intptr_t tree_build(TreeStore * store);
Node * tree_split_empty_leaf(TreeStore * store, Node * parent);
void tree_terminate(TreeStore * store);

Node * tree_locate_down_fckey(TreeStore * store, Node * start, fckey_t * key);
Node * tree_locate_up_image(TreeStore * store, Node * start, Node * needle);

int tree_node_childpos_fckey(Node * node, fckey_t * key);
int tree_node_contains_fckey(Node * node, fckey_t * key);
int tree_node_contains_node(Node * node, Node * needle);

Node * tree_iter_init(TreeIter * iter, TreeStore * store, Node * root);
Node * tree_iter_next(TreeIter * iter);
void tree_iter_update_current(TreeIter * iter, Node * node);
Node * tree_iter_next_sibling(TreeIter * iter);

