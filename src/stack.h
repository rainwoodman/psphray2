/*
 * used by tree.c
 * a growable stack.
 * */
typedef struct {
    void * data;
    int elsize;
    intptr_t len;
    intptr_t size;
} Stack;

#define stack_push(stack, type) ((type*) stack_push1(stack))
#define stack_pop(stack, type) ((type*) stack_pop1(stack))
#define stack_peek(stack, type) ((type*) stack_peek1(stack))
#define stack_get(stack, i, type) ((type*) stack_get1(stack, i))
#define stack_destroy(stack) ((stack)->data?g_free((stack)->data):0, (stack)->data = NULL)
#define stack_init(stack, type) stack_init1(stack, sizeof(type))

static inline void * stack_steal(Stack * stack, size_t * len) {
    stack->data = NULL;
    if(len) *len = stack->len;
    stack->size = 0;
    stack->len = 0;
}

static inline void stack_init1(Stack * stack, size_t elsize) {
    stack->data = NULL;
    stack->elsize = elsize;
    stack->len = 0;
    stack->size = 0;
}

static inline void * stack_push1(Stack * stack) {
    if(stack->data == NULL) {
        stack->len = 0;
        stack->size = 0;
    } 
    if(stack->len == stack->size) {
        if(stack->size == 0) stack->size += 1024;
        else if(stack->size >= 1048576) {
            stack->size += 1048576;
        } else {
            stack->size *= 2;
        }
        stack->data = g_realloc(stack->data, stack->size * stack->elsize);
    }
    void * rt = (void*)(((char*)stack->data) + (stack->len ++) * stack->elsize);
    /* shut up the compiler. we do not want string.h */
    void * memset(void*, int, size_t);
    memset(rt, 0, stack->elsize);
    return rt;
}

static inline void * stack_peek1(Stack * stack) {
    g_assert(stack->len > 0);
    return (void*)(((char*)stack->data) + (stack->len - 1) * stack->elsize);
}
static inline void * stack_get1(Stack * stack, intptr_t i) {
    return (void*)(((char*)stack->data) + i * stack->elsize);
}

static inline void * stack_pop1(Stack * stack) {
    g_assert(stack->len > 0);
    return (void*)(((char*)stack->data) + (--stack->len) * stack->elsize);
}

