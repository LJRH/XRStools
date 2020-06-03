#include<map>
#include<vector>
// #include<openmpi/mpi.h>

void lutprod(
	     int n_1, float *lut1,
	     int n_2, float *lut2,
	     int na2,
	     int nb2,
	     int dim1,
	     int dim2,
	     float *reponse_pixel,
	     int     &n_result,
	     float * &result,
	     float * rois
	     );

void lutprod4reponse(
		     int n_1, float *lut1,
		     int n_2, float *lut2,
		     int na1,
		     int na2,
		     int nb2,
		     int dim1,
		     int dim2,
		     float *reponse_pixel,
		     int dim_sol ,
		     float *solution,
		     int     &n_result,
		     float * &result,
		     int simmetrizza,
		     float * rois
		     );
  
