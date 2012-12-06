typedef struct {
    double a;
    double h;
    int64_t NgasTotal;
    int32_t Ngas;
    int N[6];
    int Nfile;
    double BoxSize;
    int double_precision;
    double masstab[6];
} SnapHeader;

void snapshot_read();

extern PSystem * PAR_BUFFER_IN;


