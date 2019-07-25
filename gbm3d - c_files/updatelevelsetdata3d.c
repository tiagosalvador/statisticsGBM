#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

#define PI acos(-1.0)

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
 /* 
  *  updatelevelsetdata3d(presence,grains,ID,S,alpha,beta,option,N0)
  *
  *  CAUTION: MODIFIES MATLAB INPUT *IN PLACE*.
  */

  mxArray *indices, *grainlevsetvals, *grainconvvals1, *grainconvvals2, *locs;
  double *pindices, *plocs, *pgrainlevsetvals, *pgrainconvvals1, *pgrainconvvals2, *id, *S;
  double sum, mink, st, m, a, b, alpha, beta;
  int N,i,j,k,ell,dims,n,nograins,gind,gind2,idk,idell, option, N0;
  double temp1[1000], temp2[1000], phi[1000], minphi[1000];
  double aux; /* auxiliar variable to deal with integer overflow when N0 ~ 50000 */
  
  dims = mxGetM(prhs[0]); /* Number of pixels. */
  n = (int) sqrt(dims);   /* Dimension of grid. */
  N = mxGetM(prhs[1]);    /* Number of grains. */  

  id = (double *) mxGetData(prhs[2]);  /* List of grain IDs. */
  S = (double *) mxGetData(prhs[3]); /* Surface tension matrix. */
  
  alpha = mxGetScalar(prhs[4]);
  beta = mxGetScalar(prhs[5]);
  option = (int) mxGetScalar(prhs[6]);
  N0 = (int) mxGetScalar(prhs[7]);

  for (j=0;j<dims;j++){ /* Loop over pixels. */
    indices = mxGetCell(prhs[0],j); /* Grains near this pixel. */
    pindices = (double *) mxGetData(indices);
    locs = mxGetCell(prhs[0],dims+j); /* Location of pixel in grain's data. */
    plocs = (double *) mxGetData(locs);
    nograins = mxGetN(indices); /* Number of grains near this pixel. */
    
    if (nograins > 1000)
            printf("nograins = %d\n",nograins);
    
    for (k=0;k<nograins;k++){ /* Loop over grains. */
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainconvvals1 = mxGetCell(prhs[1],2*N+gind-1);
      pgrainconvvals1 = (double *) mxGetData(grainconvvals1);
      grainconvvals2 = mxGetCell(prhs[1],3*N+gind-1);
      pgrainconvvals2 = (double *) mxGetData(grainconvvals2);
      temp1[k] = pgrainconvvals1[i];
      temp2[k] = pgrainconvvals2[i];
    }

    /* These lines implement the redistribution step in Esedoglu-Otto algorithm: */
    /* Form the "phi" functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      sum = 0.0;
      idk = (int) id[gind-1]; /* id of the k-th grain in the local list. */
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          gind = (int) pindices[ell]; /* Index of grain in list of all grains. */
          idell = (int) id[gind-1]; /* id of the ell-th grain in the local list. */
          aux = (double) min(idk,idell)*(2*N0-min(idk,idell)-1)/2-N0-1+max(idk,idell);
          st = S[(int) aux];
          /*st = S[(idk-1)*N0+idell-1];*/
          switch (option)
          {
              case 1:{
                  m = 1.0;
                  break;
              }
              case 2:{
                  m  = 1/st;
                  break;
              }
          }
          a = sqrt(PI)*sqrt(alpha)/(alpha-beta)*(st-beta/m);
          b = sqrt(PI)*sqrt(beta)/(alpha-beta)*(-st+alpha/m);
          sum = sum + a*temp1[ell] + b*temp2[ell];
        }
      }
      phi[k] = sum;
    }
    
    /* Minimization over the "phi" functions involved in forming level set functions: */
    for (k=0;k<nograins;k++){
      mink = 1e100;
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          mink = min( mink , phi[ell] );
        }
      }
      minphi[k] = mink;
    }    

    /* Form the level set functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainlevsetvals = mxGetCell(prhs[1],N+gind-1);
      pgrainlevsetvals = (double *) mxGetData(grainlevsetvals);

      pgrainlevsetvals[i] = (minphi[k] - phi[k]);
      if (nograins==1) pgrainlevsetvals[i] = temp1[k]; /*why?????*/
    }
    
  } /* (for j). Loop over pixels ends. */
    
}