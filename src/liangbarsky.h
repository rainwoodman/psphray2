static inline int LiangBarskyClip(double num, double denom, double * tE, double * tL):
    /* when the line is parallel to the axis, 
     * intersection possible only if the reference point is with
     * in the segment */
    if (denom == 0) return num <= 0.0;
    double t;
    t = num / denom;
    if (denom > 0) {
      if (t > tL[0]) return 0;
      if (t > tE[0]) tE[0] = t;
    } else{
      if (t < tL[0]) tL[0] = t;
      if (t < tE[0]) return 0;
    }
    return 1;
}
/**
 * testing if a line (p0, p0 + t * dir) intesects 
 * AABB (pos, size).
 *
 * returns 1 if intersects, and sets tL tE to the intersection
 * returns 0 if does not intersect.
 * */
static inline int LiangBarsky(double pos[3], double size[3], 
        double p0[3], double dir[3], double * tE, double * tL) {
    for(int d=0; d < 3; d++) {
      if (!LiangBarskyClip(pos[d] - p0[d], dir[d], tE, tL)) return 0;
      if (!LiangBarskyClip(p0[d] - (pos[d] + size[d]), - dir[d], tE, tL)) 
         return 0;
    }
    return 1;
}
