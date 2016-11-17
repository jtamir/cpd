
#include <stdio.h>
#include <stdarg.h>
#include <string.h> /* memset */
#include "misc.h"
#include <stdlib.h>
#include <math.h>

#include "misc/misc.h"


#ifdef __cplusplus
extern "C" {
#endif



/* convert subscript to column major linear index */
long
sub2ind(const unsigned int D, const long *strides, const long *sub)
{
    long res = 0;
    unsigned int i;
    for( i = 0 ; i < D ; i++ ){
        res += sub[i]*strides[i];
    }
    return res;
}

/* mul2 */
void mul2(long N, double *dst, const int *src1,  const double *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i]*src2[i];
    }
}



/* sum */
int sumi(const long N, const int *src)
{
    int res = 0.0f;
    long i;
    for( i = 0 ; i < N ; i++ ){
        res += src[i];
    }
    return res;
}



/* convert linear index to MD subscript */
void
ind2sub(const unsigned int D, const long *dims, long *sub, const long ind)
{
    
    if( D == 1 ){
        sub[0] = ind;
    }else if( D == 2){
        sub[1] = ind / dims[0];
        sub[0] = ind - sub[1]*dims[0];
    }else if( D == 3){
        sub[2] = ind / (dims[0]*dims[1]);
        sub[1] = (ind - sub[2]*dims[0]*dims[1]) / dims[0];
        sub[0] = (ind - sub[2]*dims[0]*dims[1] - sub[1]*dims[0]);
    }
}


static void rpermute(const long n, long *a) {
    long k;
    for (k = 0; k < n; k++)
        a[k] = k;
    for (k = n-1; k > 0; k--) {
        long j = rand() % (k+1);
        long temp = a[j];
        a[j] = a[k];
        a[k] = temp;
    }
}

/* randomly permutes the elements of the array perm of length n */
void randperm( const long n, long perm[])
{
    long *ind = xmalloc(n*sizeof(long));
    long *tmp = xmalloc(n*sizeof(long));
    long k;
    rpermute(n, ind);
    for( k = 0 ; k < n ; k++ ){
        tmp[k] = perm[ind[k]];
    }
    for( k = 0 ; k < n ; k++ ){
        perm[k] = tmp[k];
    }
    free(ind);
    free(tmp);
}
                    
/* returns if samples is in bounds */
int in_bounds( const long D, long *sample, const long *dims ){
    if( D == 1 ){
        return (sample[0] >= 0 && sample[0] < dims[0]);
    }else{
        return (sample[0] >= 0 && sample[0] < dims[0])
                  && in_bounds(D-1, &sample[1], &dims[1]);
    }
}


/* mod function */
long mod(const long x, const long y){
    long r = x % y;
    return r < 0 ? r + y : r;
}


#ifdef __cplusplus
}
#endif

