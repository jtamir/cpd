#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "num/multind.h"

#include "misc/mri.h"
#include "misc/debug.h"
#include "misc/misc.h"

#include "misc.h"
#include "udcpd.h"
#include "vdcpd.h"

#ifdef __cplusplus
extern "C" {
#endif


void genVDCPD(const long dims[], complex float *pattern_cfl, 
              const int *feasiblePointsD, const float FOVRatio, 
              const float C, const int shapeOpt, const float initial_mindist, 
              const float vd_exp, const int maxR){
    assert(FOVRatio > 0);
    assert(maxR > 0);
    assert(vd_exp > 0);
    int min_region_dim = 2; // Minimimum dimension of a region
    int skip_reg = (int) ceilf((double) min_region_dim / ((double) MIN(dims[0], dims[1]) / (double) maxR));
    debug_printf(DP_INFO, "---------- Variable Density CPD -------- \n");
    debug_printf(DP_INFO, "Rmax:\t%d\n", maxR);
    debug_printf(DP_INFO, "FOVRatio:\t%f\n", FOVRatio);
    debug_printf(DP_INFO, "Density falloff:\t1/%f\n", vd_exp);
    debug_printf(DP_INFO, "C:\t%f\n", C);
    debug_printf(DP_INFO, "Shape option:\t%d\n", shapeOpt);
    debug_printf(DP_INFO, "Region skip:\t%d\n", skip_reg);
    debug_printf(DP_INFO, "---------------------------------------- \n");
    

    long f_size = dims[0]*dims[1];
    double alpha_sq = pow( (double) dims[Y_DIM] / (double) dims[Z_DIM], 2);
    int *is_sampled = (int *) xmalloc(dims[Y_DIM]*dims[Z_DIM]*sizeof(int));
    double * kr = (double *) xmalloc(dims[Y_DIM]*dims[Z_DIM]*sizeof(double));
    int *region     = (int *) xmalloc(dims[Y_DIM]*dims[Z_DIM]*sizeof(int));
    debug_printf(DP_INFO, "Region assignments...\n");
    const int corner_cut = 1;
    int iky, ikz;
    for( iky = 0; iky < dims[Y_DIM] ; iky++ )
    for( ikz = 0; ikz < dims[Z_DIM] ; ikz++ ){
        long ik = iky + ikz * dims[Y_DIM];
        kr[ik] = sqrt( pow( (double) iky - (double) dims[Y_DIM] / 2.0,2) + alpha_sq*pow((double) ikz - (double) dims[Z_DIM] / 2.0,2));
        if( !corner_cut || kr[ik] <= (double) dims[Y_DIM]/2 ){
            region[ik] = MIN(1 + (int) (skip_reg * floorf( (double) maxR / (double) skip_reg * pow(kr[ik]/(dims[Y_DIM]/2),vd_exp))), maxR);
        }else{
            // Do not sample
            region[ik] = -1;
        }
    }
    int *feasiblePointsReg = xmalloc(f_size*sizeof(int));
    memset(pattern_cfl, 0, sizeof(complex float)*dims[Y_DIM]*dims[Z_DIM]*dims[T_DIM]);
    debug_printf(DP_INFO, "Region-wise CPD\n");
    int reg;
    for( reg = 1 ; reg <= maxR ; reg++ ){
        debug_printf(DP_INFO, "%d/%d\n", reg, maxR);
        complex float *masks_reg = xmalloc(f_size*reg*sizeof(complex float));
        memcpy(feasiblePointsReg, feasiblePointsD, f_size*sizeof(int));
        /*
        */
        long ik, t;
        for( ik = 0 ; ik < f_size ; ik++ ){
            feasiblePointsReg[ik] = (region[ik] == reg);
        }
        long dims_tmp[DIMS];
        md_copy_dims(DIMS, dims_tmp, dims);
        dims_tmp[T_DIM] = reg;
        
        genUDCPD(dims_tmp, masks_reg, feasiblePointsReg, FOVRatio, C, shapeOpt, initial_mindist);
        debug_printf(DP_INFO, "Copy to pat\n");
        for( t = 0 ; t < dims[T_DIM] ; t++ )
        for( ik = 0 ; ik < dims[Y_DIM]*dims[Z_DIM] ; ik++ ){
            if( masks_reg[ik + dims[Y_DIM]*dims[Z_DIM]*(t % reg)] )
                pattern_cfl[ik + dims[Y_DIM]*dims[Z_DIM]*t] = 1;
        }
        /*
        */
        free(masks_reg);
    }
    debug_printf(DP_INFO, "Cleanup...\n");
    free(feasiblePointsReg);
    free(kr); 
    free(is_sampled);
    free(region);
}

#ifdef __cplusplus
}
#endif
    
