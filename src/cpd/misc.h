long sub2ind(const unsigned int D, const long *strides, const long *sub);

void mul2(long N, double *dst, const int *src1,  const double *src2);
int sumi(const long N, const int *src);
void ind2sub(const unsigned int D, const long *dims, long *sub, const long ind);

void randperm( const long n, long perm[]);
int in_bounds( const long D, long *sample, const long *dims );
long mod( const long t, const long nt);

