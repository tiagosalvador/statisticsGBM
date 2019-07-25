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
    MATLAB INTERFACE:
      presence = get_nhd_grains(grains,dim1*dim2*dim3)
    where dim1, dim2 and dim3 are the dimensions of the grid.
    The output "presence" is an dim1*dim2*dim3-by-2 cell array:
      presence{index,1} = Grains in a nhd. of pixel at loc. index.
      presence{index,2} = Used to locate this pixel in each grains's 
                          data structure.
    Example:
      If k=presence{index,1}(i) and j=presence{index,2}(i), then
      grains{k,1}(j) = index. 
  */
  mxArray *gind, *gloc, *cell1, *cell2;
  mwSize ndim, vecdim[2];
  int i,j,k,m,n,N,dims;
  int **ind, **loc;
  int *nopresent; /* Number of grains present in a nhd. of each pixel. */
  double *pgind, *pgloc; /* C pointer versions of gind, gloc, gcval. */
  int index, location;
  double *p1, *p2;
  int *limit;
  
  N = mxGetM(prhs[0]); /* Number of grains. */
  dims = (int) mxGetScalar(prhs[1]); /* Number of pixels. */

  /* Local, C version of data structure: */
  ind = (int **) calloc(dims,sizeof(int *));
  loc = (int **) calloc(dims,sizeof(int *));
  nopresent = (int *) calloc(dims,sizeof(int));
  limit = (int *) calloc(dims,sizeof(int));
  
  /* Initial allocation: 5 neighbors. */
  for (i=0;i<dims;i++){
    limit[i] = 5;
    ind[i] = (int *) calloc(limit[i],sizeof(int));
    loc[i] = (int *) calloc(limit[i],sizeof(int));
  }
  
  /* Loop over grains: */
  for (k=0;k<N;k++){
    gind = mxGetCell(prhs[0],k);
    pgind = (double *) mxGetData(gind);
        
    m = mxGetM(gind); /* Number of pixels in buffered grain */
    
    for (i=0;i<m;i++) { /* Loop over pixels in buffered grain. */
      index = pgind[i]-1;  /* Grid index of current pixel. */
      j = nopresent[index]; /* Number of grains already at pixel. */            
      /* Reallocate memory, in chunks, if needed: */
      if (j>=limit[index]) {
        limit[index] = limit[index] + 5;
        ind[index] = realloc( ind[index], limit[index] * sizeof(int) );
        loc[index] = realloc( loc[index], limit[index] * sizeof(int) );
      } /* end if */
      ind[index][j] = k;              /* k-th grain is near pix. "index". */
      loc[index][j] = i;              /* Pix. "index" in i-th loc. */
      nopresent[index] = j+1;
    } /* end for i. */
  } /* end for k. */
  
  /* Copy data to MATLAB structure: */
  ndim = 2;
  vecdim[0] = dims;
  vecdim[1] = 2;
  plhs[0] = mxCreateCellArray(ndim,vecdim);
  for (index=0;index<dims;index++){
    
    m = nopresent[index];
    cell1 = mxCreateDoubleMatrix(1,m,mxREAL);
    cell2 = mxCreateDoubleMatrix(1,m,mxREAL);
    p1 = mxGetPr(cell1);
    p2 = mxGetPr(cell2);
    
    /* Place data in MATLAB data structure: */
    for (i=0;i<m;i++){
      p1[i] = ind[index][i]+1;
      p2[i] = loc[index][i]+1;
    }
    mxSetCell(plhs[0],index,cell1);
    mxSetCell(plhs[0],index+dims,cell2);
  }      

  /* Free dynamically allocated memory: */    
  for (i=0;i<dims;i++){
    free(ind[i]);
    free(loc[i]);
  }
  free(ind);
  free(loc);
  free(nopresent);
  free(limit); 
}