typedef struct {
    void * data;
    int elsize;
    intptr_t len;
    intptr_t size;
} Stack;

#define stack_push(stack, type) ((type*) stack_push1(stack))
#define stack_pop(stack, type) ((type*) stack_pop1(stack))
#define stack_peek(stack, type) ((type*) stack_peek1(stack))
#define stack_destroy(stack) (g_free((stack)->data), (stack)->data = NULL)
#define stack_init(stack, type) stack_init1(stack, sizeof(type))

static void stack_init1(Stack * stack, size_t elsize) {
    stack->data = NULL;
    stack->elsize = elsize;
    stack->len = 0;
    stack->size = 0;
}

static void * stack_push1(Stack * stack) {
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
    memset(rt, 0, stack->elsize);
    return rt;
}

static void * stack_peek1(Stack * stack) {
    g_assert(stack->len > 0);
    return (void*)(((char*)stack->data) + (stack->len - 1) * stack->elsize);
}

static void * stack_pop1(Stack * stack) {
    g_assert(stack->len > 0);
    return (void*)(((char*)stack->data) + (--stack->len) * stack->elsize);
}

