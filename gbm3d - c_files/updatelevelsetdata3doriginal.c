#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
 /* 
  *  updatelevelsetdata3doriginal(presence,grains,ID,S,N0);
  *
  *  CAUTION: MODIFIES MATLAB INPUT *IN PLACE*.
  */

  mxArray *indices, *grainlevsetvals, *grainconvvals, *locs;
  double *pindices, *plocs, *pgrainlevsetvals, *pgrainconvvals, *id, *S;
  double sum, mink, st;
  int N,i,j,k,ell,dims,nograins,gind,gind2,idk,idell,N0;
  double temp[1000], phi[1000], minphi[1000];
  double aux; /* auxiliar variable to deal with integer overflow when N0 ~ 50000 */
  
  dims = mxGetM(prhs[0]); /* Number of pixels. */
  N = mxGetM(prhs[1]);    /* Number of grains. */  

  id = (double *) mxGetData(prhs[2]);  /* List of grain IDs. */
  S = (double *) mxGetData(prhs[3]); /* Surface tension matrix. */
  
  N0 = (int) mxGetScalar(prhs[4]);
  
  for (j=0;j<dims;j++){ /* Loop over pixels. */
    indices = mxGetCell(prhs[0],j); /* Grains near this pixel. */
    pindices = (double *) mxGetData(indices);
    locs = mxGetCell(prhs[0],2*dims+j); /* Location of pixel in grain's data. */
    plocs = (double *) mxGetData(locs);
    nograins = mxGetN(indices); /* Number of grains near this pixel. */
   
    for (k=0;k<nograins;k++){ /* Loop over grains. */
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainconvvals = mxGetCell(prhs[1],2*N+gind-1);
      pgrainconvvals = (double *) mxGetData(grainconvvals);
      temp[k] = pgrainconvvals[i];
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
          sum = sum + st*temp[ell];
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
      if (nograins==1) pgrainlevsetvals[i] = temp[k];
    }
    
  } /* (for j). Loop over pixels ends. */
    
}
