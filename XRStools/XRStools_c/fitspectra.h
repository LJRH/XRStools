void spectra_roi_by_line_c( int Nsp,
			  double *spettro_byline,
			  double* error_sum,
			  double* frequencies_sum,
			  int Nenes,
			  float* enes_data,
			  float* denominator ,
			  int ny,
			    int nx,
			  float* mms,
			  float* MASK,

			  float mine,
			  float de,
			  float discard_threshold ,
			  float threshold_fraction,
			  float hline,
			  float slopeline,
			  float DHoverDI,
			  float deltaEref,
			  int useref,
			  int weight_by_response,
			  float Xintercept ,
			  float CRX ,
			  float fNmiddle,
			  float Xslope ,
			  int Nresp,
			  double* response_line_intensity
			  ) ;

