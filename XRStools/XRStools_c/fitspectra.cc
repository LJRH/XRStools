#include<math.h>  
#include<stdio.h>
#include<stdlib.h>
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
			  )  {

  {
    static int pass=0;
      if(!pass) {
	printf(" DENTRO spectra_roi_by_line_c \n" );
	pass++;
      }
  }
  int mask_npix=0;
  for(int i=0; i< ny*nx; i++ ) {
    if(MASK[i]) mask_npix++;
  }
  for(int iE=0; iE<Nenes; iE++) {

    {
      static int pass=0;
	if(!pass) {
	  printf(" DENTRO spectra_roi_by_line_c   Nenes %d \n" , Nenes);
	  pass++;
	}
    }
    
    float *mm = mms + iE*ny*nx;
    float ene = enes_data[iE];
    float deno = denominator[iE];
    
    int mm_npix=0;
    for(int i=0; i< ny*nx; i++ ) {
      if(mm[i]) mm_npix++;
    }   
    if(discard_threshold>0) {
      if(mm_npix>threshold_fraction*mask_npix) {
	continue;
      }
    }
    // one energy, one 2D image from the stack, one value for the denominator

    for(int iy=0; iy< ny ; iy++) {
      for(int ix=0; ix< nx ; ix++) {

	{
	  static int pass=0;
	    if(!pass) {
	      printf(" DENTRO spectra_roi_by_line_c   ny, nx %d %d\n" , ny, nx);
	      pass++;
	    }
	}
	
	if( MASK[ iy*nx+ix]==0 )  continue;


	{
	  static int pass=0;
	    if(!pass) {
	      printf(" DENTRO spectra_roi_by_line_c  MASK OK  ny, nx %d %d\n" , ny, nx);
	      pass++;
	    }
	}


	
	// Assigning an energy to the pixel
	float H0 = iy- ( hline + ix * slopeline) ; // distance in pixels from the line
	float ii = H0 /DHoverDI                  ; // The above distance converted in energy steps
	float  E = ene - ii *deltaEref                  ; // Corrected energy
	
	//  calculating the longitudinal position along the line ( to weight with the projected response_line_intensity)
	float freq;
	if( useref &&  weight_by_response) {
	  int posx_inref =   (int) round(   (ix-( (Xintercept-CRX+ fNmiddle*Xslope ) + ii * Xslope  )   )) ;
	  if (posx_inref>=0 and posx_inref< Nresp ) {
	    freq = response_line_intensity[posx_inref]; //  a weigthing factor ;
	  } else {
	    freq = 0;
	  }
	}else {
	  freq = 1.0;
	}

	{
	  static int pass=0;
	    if(!pass) {
	      printf(" DENTRO spectra_roi_by_line_c freq %e\n" , freq);
	      pass++;
	    }
	}

	
	if(freq) {
	  float fipos = (E-mine)/de  ;                      //  the position in pixel units of the contribution to the spectra array (which starts from mine)
	  int ipos  = (int)fipos      ;                  // the integer part of fipos
	  float f = fipos - ipos  ;                         // the fractional residu of fipos


	  // printf(" ,mm[iy*nx+ix] %e , freq %e , ipos %d , Nenes %d    \n",mm[iy*nx+ix] , freq, ipos, Nenes      );
	  
	  if( ipos>0 && ipos < Nsp   ) {              // If I am within the range of the spectra
	    frequencies_sum[ipos] += (1-f)* freq  ;    //     I distribute the contribution : 100% if f=0, to ipos , with the weigth given by response
	    spettro_byline [ipos] += (1-f)* mm[iy*nx+ix]; // Same thing : distributing intensity to spectro_..
	    
	    frequencies_sum[ipos+1] +=  f*freq ;        // Same thing as above, just 100% if f=1 because we are distributing to the upper pixel
	    spettro_byline [ipos+1] +=  f* mm[iy*nx+ix] ; 

	    static double sum=0;
	    sum+=  mm[iy*nx+ix]    ; 
	    printf(" SUM %e \n", sum);
	    // Calculating the error by hoping that the final result  be  gaussian
	    error_sum[ipos] += (1-f)*(1-f)* mm[iy*nx+ix] /deno   ;
	    error_sum[ipos+1] += f*f* mm[iy*nx+ix] /deno   ;
	  }
	}
      }
    }
  }
}

