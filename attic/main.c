#include <glib.h>
#include <stdint.h>
	
#include "par.h"

int main(int argc, char * argv[]) {
    PStore * pstore = pstore_new(0xff);
    int i;
    g_random_set_seed(0);
    for(i = 0; i < 1000000; i++) {
        ipos_t ipos[3];
        ipos[0] = g_random_int_range(0, IPOS_LIMIT);
        ipos[1] = g_random_int_range(0, IPOS_LIMIT);
        ipos[2] = g_random_int_range(0, IPOS_LIMIT);
    //    g_message("inserting %d", i);
        pstore_insert(pstore, ipos, 0);
    }
    Par * p = pstore_get(pstore, 0);
    g_message(PAR_FMT, PAR_PRINT(p[0])); 
    pstore_remove(pstore, p);
    p = pstore_get(pstore, 0);
    g_message(PAR_FMT, PAR_PRINT(p[0])); 
	return 0;
}
