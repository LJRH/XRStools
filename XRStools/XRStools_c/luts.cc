#include<map>
#include<vector>
#include<stdio.h>
#include<string.h>
#include<math.h>
// #include<mpi.h>

#include<emmintrin.h>
#define FLOAT_TO_INT(out,in)  \
     out=_mm_cvtss_si32(_mm_load_ss(&(in)));


#define max(a,b) (((a)>(b))? (a):(b))
#define min(a,b) (((a)<(b))? (a):(b))


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
	     ) {

  std::vector<float> res;
  int i1,j1;
  float y0,y1;
  float Y0,Y1;
  int iY0, iY1;
  float  fiY0, fiY1;
  float f1,f2;
  // printf("qui\n");
  // double SSUM = 0;
  for(int il1 =0; il1<n_1; il1++ ) {

    FLOAT_TO_INT(i1,lut1[il1*5+0   ]);
    FLOAT_TO_INT(j1,lut1[il1*5+1   ]);
    
    f1 = lut1[il1*5+2   ];
    y0 = lut1[il1*5+3   ];
    y1 = lut1[il1*5+4   ];
    
    Y0 = dim1*y0;
    Y1 = dim1*y1;
    
    float tmp;
    FLOAT_TO_INT(  iY0 ,  tmp=ceil(Y0)   );
    FLOAT_TO_INT(  iY1 ,  tmp=floor(Y1)   );



    fiY0 = min(  iY0, Y1 );
    fiY1 = max(  iY1, Y0 );


    
    // facts = (reponse_pixel[ iY0:iY1]).sum(axis=0)
    {
      int i2,j2;
      float x0,x1;
      float X0,X1;
      int iX0, iX1;
      float  fiX0, fiX1;
      for(int il2 =0; il2<n_2; il2++ ) {
	
	FLOAT_TO_INT(i2,lut2[il2*5+0   ]);
	FLOAT_TO_INT(j2,lut2[il2*5+1   ]);


	if(rois[ i1*na2+i2 ]!=1.0) {
	  // printf("SALTO roi %e i1,i2, na1, na2    %d %d %d   \n" ,rois[ i1*na2+i2 ] , i1,i2,  na2 );
	  continue;
	} else {
	  // printf("VADO  roi %e i1,i2, na2    %d %d %d   \n" ,rois[ i1*na2+i2 ] , i1,i2,  na2 );
	}
	
	f2 = lut2[il2*5+2   ];
	x0 = lut2[il2*5+3   ];
	x1 = lut2[il2*5+4   ];
	
	X0 = dim2*x0;
	X1 = dim2*x1;
	float tmp;
	FLOAT_TO_INT(  iX0 ,  tmp=ceil(X0)   );
	FLOAT_TO_INT(  iX1 ,  tmp=floor(X1)   );
	
	float Fatt =0.0;
	fiX0 = min(  iX0, X1 );
	fiX1 = max(  iX1, X0 );
	
	
	float yfact;	  
	for(int  idy=iY0-1; idy<iY1+1; idy++) {

	  yfact=1.0;
	  if(idy==iY0-1) {
	    yfact = (fiY0-Y0);
	    if (yfact< 1.0e-8 ) continue;
	  }
	  if(idy==iY1) {
	    yfact = (Y1-fiY1);
	    if (yfact< 1.0e-8 ) continue;
	  }

	  // if ( idy==25 ) {
	  //   printf("  idy, idx , dim1, dim2  %d %d %d %d\n", idy, Y0, Y1, y0, y1, iY0, iY1);
	  //   exit(0);
	  // }

	  
	  for(int  idx=iX0; idx<iX1; idx++) {

	    if ( idy<0 || idy>= dim1 ||  idx<0 || idx>= dim2  ) {
	      printf("  idy, idx , dim1, dim2  %d %d %d %d\n", idy, idx , dim1, dim2);
	    } 
	    
	    Fatt +=  reponse_pixel[ idy*dim2  +    idx ]*yfact ; 
	  }
	  
	  
	  if ( fiX0-X0>1.0e-8) {
	    int where ;
	    float tmp;
	    FLOAT_TO_INT(where , tmp=floor(X0)) ; 
	    Fatt += (fiX0-X0)  *  reponse_pixel[ idy*dim2  +    where ] *yfact ; 
	  }
	  if(iX1>=iX0) {
	    if(X1-fiX1>1.0e-8) {
	      int where = iX1;
	      Fatt += (X1-fiX1)    *  reponse_pixel[ idy*dim2  +    where ] *yfact; 
	    }
	  }
	}
	
	Fatt /= (X1-X0)*(Y1-Y0)  ;
	res.push_back(    i1*na2+i2    );
	res.push_back(    j1*nb2+j2    );
	res.push_back(    f1*f2* Fatt    );
	// SSUM += f1*f2* Fatt ; 
	// SSUM += Fatt ; 
      }
    }
  }
  // printf("qua %e \n", SSUM);

  n_result =  res.size()/3;
  result = new float [n_result*3];
  memcpy(result, &(res[0]) , n_result*3*sizeof(float));



  
}

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
		     ) {


#define MATRIX(i,j,k,l)  matrix[     ((( (i)  )*na2   +(j)   )*dim1   + (k) )*dim2     +(l)  ]   
  std::vector<float>  matrix( na1*na2* dim1*dim2 ,0.0 ) ;

  std::vector<float> res;
  int i1,j1;
  float y0,y1;
  float Y0,Y1;
  int iY0, iY1;
  float  fiY0, fiY1;
  float f1,f2;
  // int simmetrizza = 1;
  // printf("qui\n");
  // double SSUM = 0;
  for(int il1 =0; il1<n_1; il1++ ) {
    
    FLOAT_TO_INT(i1,lut1[il1*5+0   ]);
    FLOAT_TO_INT(j1,lut1[il1*5+1   ]);
    
    f1 = lut1[il1*5+2   ];
    y0 = lut1[il1*5+3   ];
    y1 = lut1[il1*5+4   ];
    
    Y0 = dim1*y0;
    Y1 = dim1*y1;
    
    float tmp;
    FLOAT_TO_INT(  iY0 ,  tmp=ceil(Y0)   );
    FLOAT_TO_INT(  iY1 ,  tmp=floor(Y1)   );
   
    fiY0 = min(  iY0, Y1 );
    fiY1 = max(  iY1, Y0 );
    
    // facts = (reponse_pixel[ iY0:iY1]).sum(axis=0)
    {
      int i2,j2;
      float x0,x1;
      float X0,X1;
      int iX0, iX1;
      float  fiX0, fiX1;
      for(int il2 =0; il2<n_2; il2++ ) {
	
	FLOAT_TO_INT(i2,lut2[il2*5+0   ]);
	FLOAT_TO_INT(j2,lut2[il2*5+1   ]);
	

	if( fabs(solution[ j1*nb2+j2 ])<1.0e-30) continue;
       
	
	f2 = lut2[il2*5+2   ];
	x0 = lut2[il2*5+3   ];
	x1 = lut2[il2*5+4   ];
	
	X0 = dim2*x0;
	X1 = dim2*x1;
	float tmp;
	FLOAT_TO_INT(  iX0 ,  tmp=ceil(X0)   );
	FLOAT_TO_INT(  iX1 ,  tmp=floor(X1)   );
	
	float Fatt =0.0;
	fiX0 = min(  iX0, X1 );
	fiX1 = max(  iX1, X0 );
	
	
	float yfact;	  
	for(int  idy=iY0-1; idy<iY1+1; idy++) {
	  
	  yfact=1.0/( (X1-X0)*(Y1-Y0));
	  if(idy==iY0-1) {
	    yfact = (fiY0-Y0)/( (X1-X0)*(Y1-Y0));
	    if (yfact< 1.0e-8 ) continue;
	  }
	  if(idy==iY1) {
	    yfact = (Y1-fiY1)/( (X1-X0)*(Y1-Y0));
	    if (yfact< 1.0e-8 ) continue;
	  }
	  
	  // if ( idy==25 ) {
	  //   printf("  idy, idx , dim1, dim2  %d %d %d %d\n", idy, Y0, Y1, y0, y1, iY0, iY1);
	  //   exit(0);
	  // }

	  yfact *= f1*f2;
	  if (simmetrizza) yfact/=8;

	  if (yfact< 1.0e-10 ) continue;
	  
	  for(int  idx=iX0; idx<iX1; idx++) {
	    
	    if ( idy<0 || idy>= dim1 ||  idx<0 || idx>= dim2  ) {
	      printf("CHE E STA ROBA  ??? idy, idx , dim1, dim2  %d %d %d %d %e\n", idy, idx , dim1, dim2, yfact);
	    } else {

	      float term =    solution[ j1*nb2+j2    ]*  yfact ;
	      MATRIX( i1,i2,  idy,idx   ) += term  ;
	      if(simmetrizza) {
		MATRIX( i1,i2,  dim1-1 -idy,idx   ) += term  ;
		MATRIX( i1,i2,  idy,dim1-1 -idx   ) += term  ;
		MATRIX( i1,i2,  dim1-1 -idy,dim1-1 -idx   ) += term  ;
		
		MATRIX( i1,i2,  idx,idy   ) += term  ;
		MATRIX( i1,i2,  dim1-1 -idx,idy   ) += term  ;
		MATRIX( i1,i2,  idx,dim1-1 -idy   ) += term  ;
		MATRIX( i1,i2,  dim1-1 -idx,dim1-1 -idy   ) += term  ;
		
	      }
	    }
	    // Fatt +=  reponse_pixel[ idy*dim2  +    idx ]*yfact ; 
	  }
	  	  
	  if ( fiX0-X0>1.0e-8) {
	    int where ;
	    float tmp;
	    FLOAT_TO_INT(where , tmp=floor(X0)) ; 
	    // Fatt += (fiX0-X0)  *  reponse_pixel[ idy*dim2  +    where ] *yfact ;

	    if ( idy<0 || idy>= dim1 ||  where<0 || where>= dim2  ) {
	      printf("CHE E STA ROBA ???   idy, where , dim1, dim2  %d %d %d %d yfact %e\n", idy, where , dim1, dim2, yfact);
	    } else {
	      
	      
	      float term =   solution[ j1*nb2+j2    ]* (fiX0-X0)  * yfact  ; 
	      
	      MATRIX( i1,i2,  idy,where) +=   term ;
	      if(simmetrizza) {
		
		MATRIX( i1,i2,  dim1-1 -idy,where   ) += term  ;
		MATRIX( i1,i2,  idy,dim1-1 -where   ) += term  ;
		MATRIX( i1,i2,  dim1-1 -idy,dim1-1 -where   ) += term  ;
		
		MATRIX( i1,i2,  where,idy   ) += term  ;
		MATRIX( i1,i2,  dim1-1 -where,idy   ) += term  ;
		MATRIX( i1,i2,  where,dim1-1 -idy   ) += term  ;
		MATRIX( i1,i2,  dim1-1 -where,dim1-1 -idy   ) += term  ;
	      }
	    }

	  }
	  if(iX1>=iX0) {
	    if(X1-fiX1>1.0e-8) {
	      int where = iX1;
	      // Fatt += (X1-fiX1)    *  reponse_pixel[ idy*dim2  +    where ] *yfact;
	      if ( idy<0 || idy>= dim1 ||  where<0 || where>= dim2  ) {
		printf("CHE E STA ROBA   idy, where , dim1, dim2  %d %d %d %d yfact %e\n", idy, where , dim1, dim2, yfact);
	      } else {

		float term =   solution[ j1*nb2+j2    ]*  (X1-fiX1)  * yfact  ; 
	      
		MATRIX( i1,i2,  idy,where) +=  term  ;
		if(simmetrizza) {
		  
		  MATRIX( i1,i2,  dim1-1 -idy,where   ) += term  ;
		  MATRIX( i1,i2,  idy,dim1-1 -where   ) += term  ;
		  MATRIX( i1,i2,  dim1-1 -idy,dim1-1 -where   ) += term  ;
		  
		  MATRIX( i1,i2,  where,idy   ) += term  ;
		  MATRIX( i1,i2,  dim1-1 -where,idy   ) += term  ;
		  MATRIX( i1,i2,  where,dim1-1 -idy   ) += term  ;
		  MATRIX( i1,i2,  dim1-1 -where,dim1-1 -idy   ) += term  ;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  for(int i1 = 0; i1<na1 ; i1++) {
    for(int i2 = 0; i2<na2 ; i2++) {
      
      if(rois[ i1*na2+i2 ]!=1.0){
	continue;
      } else {
	// printf(" %e %d %d %d %d \n", rois[ i1*na2+i2 ], i1,i2,na1, na2);
      }
      
      for(int iy = 0; iy<dim1 ; iy++) {
	for(int ix = 0; ix<dim2 ; ix++) {
	  if( MATRIX(i1,i2,iy,ix)!=0.0) {
	    res.push_back(    i1*na2+i2    );
	    res.push_back(    iy*dim1+ix    );
	    res.push_back(    MATRIX(i1,i2,iy,ix)   );
	  }
	}
      }
    }
  }
  n_result =  res.size()/3;
  result = new float [n_result*3];
  memcpy(result, &(res[0]) , n_result*3*sizeof(float));
#undef MATRIX  
}
  


