#include "commonblock.h"
#include "liangbarsky.h"
#include "tree.h"
ray_t * R;
int NRay;

static inline void node_get_aabb(Node * node, double pos[3], double size[3]) {
    int64_t ipos[3];
    fckey_to_ipos(&node->key, ipos);
    for(int d=0; d < 3; d++) {
        pos[d] = CB.BoxSize * ipos[d] / (FCKEY_MAX + 1);
        size[d] = CB.BoxSize / (1L << (FCKEY_BITS - node->order));
    }
}

int ray_intersect_node(ray_t * ray, Node * node, double * tE, double * tL) {
    double pos[3], size[3];
    node_get_aabb(node, pos, size);
    return LiangBarsky(pos, size, ray->pos, ray->dir, tE, tL);
}

intersect_t * ray_intersect_tree(ray_t * ray, TreeStore * store, Node * start) {
    TreeIter iter = {store, start, NULL};
    double tE, tL;
    Stack stack;
    stack_init(&stack, intersect_t);
    Node * node = tree_iter_next(&iter);
    while(node) {
        if(!ray_intersect_node(ray, node, &tE, &tL)) {
            node = tree_iter_next_sibling(&iter);
            continue;
        } 
        if(node->complete) {
            intersect_t * item = stack_push(&stack, intersect_t);
            item->domainindex = 0;
            item->rayindex = 0;
            item->nodeindex = 0;
            item->tE = tE;
            item->tL = tL;
        }
        node = tree_iter_next(&iter);
    }
}
