#ifndef __SNAPSHOT_H__
#define __SNAPSHOT_H__
typedef struct {
    double a;
    double h;
    int N[6];
    int Nfile;
    unsigned int flag_double:1;
    unsigned int flag_cool:1;
    unsigned int flag_sfr:1;
    double BoxSize;
    double masstab[6];
} SnapHeader;

PackedPar * snapshot_read(SnapHeader * h);
#endif
