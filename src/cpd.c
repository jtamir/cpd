
#include <getopt.h>
#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include <math.h>

#include "num/multind.h"
#include "num/flpmath.h"

#include "misc/misc.h"
#include "misc/mri.h"
#include "misc/debug.h"
#include "misc/mmio.h"
#include "misc/opts.h"

#include "cpd/udcpd.h"
#include "cpd/vdcpd.h"

static const char* usage_str = "-t segments <input> <assignment>";
static const char* help_str = 	"UD CPD segmentation\n";

int main_cpd(int argc, char* argv[])
{
	float FOVRatio = 1;
	long nt = 4;
	float C = 1;
	float initial_mindist = 1000; // must be large enough
	long shapeOpt = CONES;
	bool vardens = false;
	int maxR = -1;
	float vd_exp = 1;

	const struct opt_s opts[] = {
		OPT_FLOAT('v', &vd_exp, "exp", "Variable density exponent"),
		OPT_INT('M', &maxR, "regs", "regions at kmax"),
		OPT_INT('d', &debug_level, "level", "debug level"),
		OPT_LONG('t', &nt, "T", "Number of segments"),
		OPT_FLOAT('F', &FOVRatio, "FOV", "FOV Ratio"),
		OPT_FLOAT('D', &initial_mindist, "dist", "Initial min distance"),
		OPT_FLOAT('C', &C, "C", "Scaling spatial vs temporal"),
		OPT_LONG('s', &shapeOpt, "S", "Shape opt"),
	};
	cmdline(&argc, argv, 2, 2, usage_str, help_str, ARRAY_SIZE(opts), opts);

	if( maxR > 0 ){
		vardens = true;
	}

	long dimsF[DIMS];
	complex float* feasiblePoints = load_cfl(argv[1], DIMS, dimsF);
	long dims[DIMS];
	md_copy_dims(DIMS, dims, dimsF);
	dims[T_DIM] = nt;

	long f_size = md_calc_size(DIMS, dimsF);
	int *feasiblePointsI = xmalloc(f_size*sizeof(int));
	for( long i = 0 ; i < f_size ; i++ ){
		feasiblePointsI[i] = (int) (feasiblePoints[i] != 0.0);
	}

	debug_print_dims(DP_INFO, DIMS, dims);
	complex float *pattern_cfl = create_cfl(argv[2], DIMS, dims);
	if( vardens ){
		genVDCPD(dims, pattern_cfl, feasiblePointsI, FOVRatio, C, shapeOpt, initial_mindist, vd_exp, maxR); 
	}else{
		genUDCPD(dims, pattern_cfl, feasiblePointsI, FOVRatio, C, shapeOpt, initial_mindist);
	}

	// cleanup
	unmap_cfl(DIMS, dimsF, feasiblePoints);
	unmap_cfl(DIMS, dims, pattern_cfl);
	free(feasiblePointsI);

	return 0;
}




