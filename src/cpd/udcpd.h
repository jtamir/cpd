#ifndef UDCPD_H
#define UDCPD_H 1

#ifdef __cplusplus
extern "C" {
#endif



#define KDIMS 2
#define KTDIMS 3
#define Y_DIM 0u
#define Z_DIM 1u
#define T_DIM 2u

#define SCALE_MINDISTANCE_K 0.9
#define SCALE_MINDISTANCE_T 0.9
#define THRESH_MINDISTANCE_K 0.25
#define THRESH_MINDISTANCE_T 0.1


/* min distance shape options */
enum shape_opt {CROSS, L1_BALL, L2_BALL, CONES, PLANE_AND_CONES};

/* min distance relaxation options */
#define DIST_K 0
#define DIST_KT 1

extern void genUDCPD(const long *dims, _Complex float *masks, const int *feasiblePoints, 
            const double FOVRatio, const double C, const long shapeOpt, const float initial_mindist);


#ifdef __cplusplus
}
#endif

#endif // UDCPD_H

