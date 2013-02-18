
/* name has to match g_key_file_get_xxxx, type is a ctype. 
 * two functions will be defile 
 * dxxxxx(block, key, targetaddress, default) :key is optional 
 * _xxxxx(block ,key, targetaddress) : key is required
 *    */
#define SCHEMA(name, type) \
static int __ ## name (GKeyFile * keyfile, char * group, char * key, type * value, type def, int required) { \
    GError * error = NULL; \
    if(!g_key_file_has_group(keyfile, group) \
        || !g_key_file_has_key(keyfile, group, key, &error)) { \
        if(error != NULL) { \
            g_error("paramfile %s/%s:%s", group, key, error->message); \
        } \
        if(required) { \
            g_error("paramfile %s/%s is required", group, key); \
        } \
        g_key_file_set_ ## name(keyfile, group, key, def); \
        value[0] = def; \
        return FALSE; \
    } \
    value[0] = g_key_file_get_ ## name(keyfile, group, key, &error); \
    if(error) { \
        g_error("error %s\n", error->message); \
    } \
    return TRUE; \
} \
static int d ## name (GKeyFile * keyfile, char * group, char * key, type * value, type def) { \
    return __ ## name(keyfile, group, key, value, def, 0); \
} \
static int _ ## name (GKeyFile * keyfile, char * group, char * key, type * value) { \
    return __ ## name(keyfile, group, key, value, (type) 0, 1); \
}

SCHEMA(double, double);
SCHEMA(integer, int);
SCHEMA(string, char *) ;
