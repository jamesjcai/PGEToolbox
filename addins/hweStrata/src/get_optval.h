#define MFLEN 80

/* possible values for vtype */
#define OV_FLAG          0
#define OV_INT           1
#define OV_LONG          2
#define OV_FLOAT         3
#define OV_DOUBLE        4
#define OV_UNSIGNED_INT  5
#define OV_UNSIGNED_LONG 6
#define OV_STRING        7 

int get_optval(char * argv, char *opt, int vtype, void *val);

char * cmdString(int argc, char ** argv);
