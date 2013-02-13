#ifndef __SNAPSHOT_H__
#define __SNAPSHOT_H__
typedef struct {
    double a;
    double h;
    int N[6];
    int Nfile;
    double BoxSize;
    int double_precision;
    double masstab[6];
} SnapHeader;

PackedPar * snapshot_read(SnapHeader * h);
#endif
