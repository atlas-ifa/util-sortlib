/* A few sorting routines */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* The matrix we want...
 *             ushort  short  int  float  double  long_double  *string(?)
 *  qsort2
 *  insort2
 *  quantile
 *  wquantile
 */

// Some test code...
/* 
  nasty=0
   ( for i in {1..20} 30 50 ; do sortlib $i 20000 $nasty ; done ; \
     for i in {1,2,3,5}00 ; do sortlib $i 5000 $nasty ; done ; \
     for i in {1,2,3,5}000 ; do sortlib $i 100 $nasty ; done ; \
     for i in {1,2,3,5}0000 ; do sortlib $i 10 $nasty ; done ; ) | tee /tmp/foo

  nasty=0
   ( for i in {1..20} 30 50 ; do sortest $i 20000 $nasty ; done ; \
     for i in {1,2,3,5}00 ; do sortest $i 5000 $nasty ; done ; \
     for i in {1,2,3,5}000 ; do sortest $i 100 $nasty ; done ; \
     for i in {1,2,3,5}0000 ; do sortest $i 10 $nasty ; done ; ) | tee /tmp/foo
*/

/* Timings on Dell XPS, random data: */
/*    qsort2:     t ~  3.9[nsec] * (N log_2(N)) */
/*    heap_sort:  t ~  3.8[nsec] * (N log_2(N)) */
/*    qsort_unix: t ~  5.7[nsec] * (N log_2(N)) */
/*    insort2:    t ~  25[nsec] * N + 0.16[nsec]*N^2 */
/*    medec:      t ~  20[nsec] * N */

/* Timings on Dell XPS, block sorted: */
/*    qsort2:     t ~  1.6[nsec] * (N log_2(N)) */
/*    heap_sort:  t ~  3.1[nsec] * (N log_2(N)) */
/*    qsort_unix: t ~  1.7[nsec] * (N log_2(N)) */
/*    insort2:    t ~  8.4[nsec] * N + 0.14[nsec]*N^2 */
/*    medec:      t ~  5.9[nsec] * N */

/* Timings on Dell XPS, already sorted: */
/*    qsort2:     t ~  0.8[nsec] * (N log_2(N)) */
/*    heap_sort:  t ~  3.2[nsec] * (N log_2(N)) */
/*    qsort_unix: t ~  1.6[nsec] * (N log_2(N)) */
/*    insort2:    t ~  0.8[nsec] * N */
/*    medec:      t ~  5.0[nsec] * N */

typedef double DATA;

void qsort(void *base, size_t nmemb, size_t size,
	   int (*compar)(const void *, const void *));

/* Pseudo-random number algorithm of Marsaglia: rand = xorshift32(&rand) */
uint32_t xorshift32(uint32_t *state);

#define TEST		/* A main to do timing testing */
#ifdef TEST
#include <sys/time.h>
#endif


#define MAXSTACK 256
#define NSTOP 15
// #define NSTOP 8
#define QSORT_NMIN 100	/* If this few, just sort by insertion */

void insort2(int n, DATA *key, int *idx);

/* Pseudo-random number algorithm of Marsaglia: rand = xorshift32(rand) */
uint32_t xorshift32(uint32_t *state)
{
   uint32_t x = *state;
   x ^= x << 13;
   x ^= x >> 17;
   x ^= x << 5;
   return x;
}

/* Quicksort, possibly carry an integer array: 1980 (John Tonry) */
void qsort2(int n, DATA *x, int *idx)
{
   DATA key, kl, kr, km, temp;
   int l, r, m, lstack[MAXSTACK], rstack[MAXSTACK], sp, randomness, itemp;
   uint32_t ran=3141592653;
   register int i, j;
   int mgtl, lgtr, rgtm;

/* Insertion sorting when n is small */
   if(n < QSORT_NMIN) {
      insort2(n, x, idx);
      return;
   }

   for(i=randomness=0; i<n-1; i++) randomness += x[i+1] > x[i];
   if(randomness > n/2) randomness = n - randomness;

   sp = 0;
   lstack[sp] = 0;
   rstack[sp] = n-1;

   while(sp >= 0) {
/* Sort a subrecord off the stack */
      l = lstack[sp];
      r = rstack[sp--];
      m = (l + r) / 2;
/* Set KEY = median of X(L), X(M), X(R) */
      kl = x[l];
      km = x[m];
      kr = x[r];
      mgtl = km > kl;
      rgtm = kr > km;
      lgtr = kl > kr;

      if(randomness < n/10) {
/* Set KEY = random point between X(L), X(R) */
	 ran = xorshift32(&ran);
	 m = l + ran / 4294967296.0 * (r-l);
	 key = x[m];

      } else {
/* Set KEY = median of X(L), X(M), X(R) */
	 if(mgtl ^ rgtm) {
	    if(mgtl ^ lgtr) key = kr;
	    else            key = kl;
	 } else {
	    key = km;
	 }
      }

      i = l;
      j = r;
      if(idx == NULL) {
	 while(1) {
/* Find a big record on the left */
	    while(x[i] < key) i++;

/* Find a small record on the right */
	    while(x[j] > key) j--;

	    if(i >= j) break;
/* Exchange records */
	    temp = x[i];
	    x[i++] = x[j];
	    x[j--] = temp;
	 }
      } else {
	 while(1) {
	    while(x[i] < key) i++;
	    while(x[j] > key) j--;
	    if(i >= j) break;
	    temp = x[i];
	    x[i++] = x[j];
	    x[j--] = temp;
	    itemp = idx[i];
	    idx[i++] = idx[j];
	    idx[j--] = itemp;
	 }
      }

/* Subfile is partitioned into two halves, left .le. right */
/* Push the two halves on the stack */
      if(j-l+1 > NSTOP) {
	 lstack[++sp] = l;
	 rstack[sp] = j;
      }
      if(sp >= MAXSTACK) break;	// Stack overflow, fall back on insertion
      if(r-j > NSTOP) {
	 lstack[++sp] = j+1;
	 rstack[sp] = r;
      }
      if(sp >= MAXSTACK) break;	// Stack overflow, fall back on insertion
   }

/* Insertion sorting to finish up */
   insort2(n, x, idx);
   return;
}

/* insertion sort, possibly carry an integer array 020130 John Tonry */
void insort2(int n, DATA *key, int *idx)
{
   register int i, j, k, itmp;
   register DATA tmp;

   if(n < 2) return;
   if(idx == NULL) {
      for(j=n-2; j>=0; j--) {
	 for(i=j+1, k=j; i<n; i++) {
	    if(key[j] <= key[i]) break;
	    k = i;
	 }
	 if(k != j) {
	    tmp = key[j];
	    memmove(key+j, key+j+1, (k-j)*sizeof(DATA));
	    key[k] = tmp;
	 }
      }
   } else {
      for(j=n-2; j>=0; j--) {
	 for(i=j+1, k=j; i<n; i++) {
	    if(key[j] <= key[i]) break;
	    k = i;
	 }
	 if(k != j) {
	    tmp = key[j];
	    memmove(key+j, key+j+1, (k-j)*sizeof(DATA));
	    key[k] = tmp;
	    itmp = idx[j];
	    memmove(idx+j, idx+j+1, (k-j)*sizeof(int));
	    idx[k] = itmp;
	 }
      }
   }
   return;
}

/* Simple median using insertion sort */
DATA median(int n, DATA *key, int *idx)
{
   DATA tmp;
   insort2(n, key, idx);
   tmp = 0.5 * (key[n/2]+key[(n-1)/2]);		/* Median */
   return(tmp);
}

/* Return weighted median (or other quantile).  If wtot!=0 do not re-sort arrays. */
DATA wmedian(int n, DATA *value, DATA *wgt, int *idx, DATA *wtot, DATA quantile)
#if 0	/* argument explanation */
   int n, 	/* number of data points */
   DATA *value,	/* values of data points: will be sorted in place if wtot==0 */
   DATA *wgt, 	/* weights of data points wgt[idx[i]] corresponds to value[i] */
   int *idx, 	/* index array: will be filled and sorted in place if wtot==0 */
   DATA *wtot, 	/* Total weight: sigma = sqrt(n/wtot).  Set to 0 initially, */
		/*   non-zero for eval of different quantiles without re-sort */
   DATA quantile/* Return this quantile: 0.5 for median */
#endif
{
   int i;
   DATA wprev, wthis, xm0, xm1, wm0, wm1;

/* Degenerate cases */
   if(n <= 0) {
      return(0.0);
   } else if(n == 1) {
      return(value[0]);
   }

/* If wtot initialized to 0, sort value and index arrays */
   if(*wtot == 0.0) {
      for(i=0; i<n; i++) {
	 *wtot += wgt[i];
	 idx[i] = i;
      }
      qsort2(n, value, idx);
   }

/* Get the weighted mean and quartiles */
   wprev = 0.0;
   xm0 = value[0];
   xm1 = value[n-1];
   wm0 = 0.0;
   wm1 = *wtot;
   for(i=0; i<n; i++) {
      wthis = wgt[idx[i]];
      if(wprev < quantile*(*wtot) && wprev+wthis >= quantile*(*wtot)) {
	 if(i > 0)   xm0 = 0.5*(value[i-1] + value[i]);
	 else        xm0 = 0.5*(3*value[i] - value[i+1]);
	 if(i < n-1) xm1 = 0.5*(value[i] + value[i+1]);
	 else        xm1 = 0.5*(-value[i-1] + 3*value[i]);
	 wm0 = wprev / (*wtot);
	 wm1 = (wprev + wthis) / (*wtot);
	 break;
      }
      wprev += wthis;
   }
   return(xm0 + (xm1-xm0) * (quantile-wm0)/(wm1-wm0));
}


/* median_sift() moves central points which straddle goal to start of array */
/* There are nnew points left in the array that must be sorted, when sorted the
 * median will be found at goal-nbelow.  nbelow is the number of points from
 * the original array that lie below these nnew new ones (and are now pushed
 * past the nnew points at the beginning of the array). */
/* Typical usage:
 *
 *    median_sift(n, buf, n/2, 10, &nnew, &nbelow);
 *    median(nnew, buf, NULL);
 *    med = 0.5*(buf[(n-1)/2-nbelow] + buf[n/2-nbelow]);
 *
 * Note that median_sift() can be called multiple times, e.g.
 *
 *    median_sift(n, buf, n/2, 10, &nnew, &nbelow);
 *    median_sift(nnew, buf,  n/2-nbelow, 10, &nnew, &nbelow2);
 *    median(nnew, buf, NULL);
 * 
 * The median will be found at goal-nbelow (and prev for even), e.g.
 *
 *    med = 0.5*(buf[(n-1)/2-nbelow-nbelow2] + buf[n/2-nbelow-nbelow2]);
 */

/* 090315 John Tonry */
static DATA *sift_posts=NULL;
static int *sift_ngelt=NULL;
static int nsift=0;
void median_sift(int n, DATA *x, int goal, int nbin,
		int *nnew, int *nbelow)
{
   int i, j, k, nabove;
   DATA tmp;

   if(n < nbin+1) {
      *nnew = n;
      *nbelow = 0;
      return;
   }

/* Allocate some buffer space */
   if(nbin > nsift) {
      if(sift_posts != NULL) free(sift_posts);
      if(sift_ngelt != NULL) free(sift_ngelt);
      sift_posts = (DATA *)calloc(nbin+1, sizeof(DATA));
      sift_ngelt = (int *)calloc(nbin+1, sizeof(int));
      nsift = nbin;
   }

/* Sample some points */
   for(i=1; i<nbin; i++) sift_posts[i] = x[(i*(n-1))/nbin];
/* and sort them */
   insort2(nbin-1, sift_posts+1, NULL);

/* Here's what we're calculating
 *    sift_posts[0]           = minimum value
 *    sift_posts[1:nbin-1]    = suggested division points amongst the values
 *    sift_posts[nbin]        = maximum value
 *    
 *    sift_ngelt[0]           = number below sift_posts[1]
 *    sift_ngelt[1:nbin-1]    = number >= sift_posts[j] and < sift_posts[j+1]
 *      (i.e. equal to sift_posts[j] goes into sift_ngelt[j])
 */

   sift_posts[0] = sift_posts[1];
   sift_posts[nbin] = sift_posts[nbin-1];
   for(i=0; i<=nbin; i++) sift_ngelt[i] = 0;
/* Count up all the points and find the min and max */
   for(i=0; i<n; i++) {
      if(x[i] < sift_posts[1]) {
	 if(x[i] < sift_posts[0]) sift_posts[0] = x[i];
	 sift_ngelt[0] += 1;
      } else if(x[i] >= sift_posts[nbin-1]) {
	 if(x[i] >= sift_posts[nbin]) sift_posts[nbin] = x[i];
	 sift_ngelt[nbin-1] += 1;
      } else {
	 for(j=1; j<nbin-1; j++) {
	    if(x[i] >= sift_posts[j] && x[i] < sift_posts[j+1]) {
	       sift_ngelt[j] += 1;
	       break;
	    }
	 }
      }
   }

/* Figure out the points which can be ignored for the median */
/* The points between sift_posts[j] and sift_posts[k] contain the goal point */
/* There are nbelow and nabove=i points which don't need to be examined */
   for(j=0, *nbelow=0; j<nbin; j++) {
      if(*nbelow + sift_ngelt[j] >= goal) break;
      *nbelow += sift_ngelt[j];
   }
   for(k=nbin, nabove=0; k>=1; k--) {
      if(nabove + sift_ngelt[k-1] >= n-goal-1) break;
      nabove += sift_ngelt[k-1];
   }
   

/* Swap the points between sift_posts[j] and sift_posts[k] to the beginning */
   for(i=0, *nnew=0; i<n; i++) {
      if(x[i] >= sift_posts[j] && x[i] < sift_posts[k]) {
	 tmp = x[i];
	 x[i] = x[*nnew];
	 x[*nnew] = tmp;
	 *nnew += 1;
      }
   }

   return;
}

/* wmedian_sift() moves points which straddle goal to start of array */
/* There are "nnew" points left in the array with "below" summed
 * weight already accounted for, so the desired goal will be found at
 * "goal-below" summed weight in the truncated array.  The weight
 * array must sum to 1. */

/* Typical usage:
 *
 *    wmedian_sift(n, buf, wgt, 0.5, 10, &nnew, &below);
 *    med = wmedian(nnew, buf, wgt, NULL, &wtot, (0.5-below)/(1-below));
 *
 * Note that median_sift() can be called multiple times, e.g.
 *
 *    median_sift(n, buf, n/2, 10, &nnew, &nbelow);
 *    median_sift(nnew, buf,  n/2-nbelow, 10, &nnew, &nbelow2);
 *    median(nnew, buf, NULL);
 * 
 * The median will be found at goal-nbelow (and prev for even), e.g.
 *
 *    med = 0.5*(buf[(n-1)/2-nbelow-nbelow2] + buf[n/2-nbelow-nbelow2]);
 */

/* 090315 John Tonry */
void wmedian_sift(int n, DATA *x, DATA *wgt, double goal, int nbin,
		int *nnew, double *below)
{
   int i, j, loend, hiend;
   DATA tmp, lo, hi, xhi;

   if(n < nbin+1) {
      *nnew = n;
      *below = 0;
      return;
   }

/* Allocate some buffer space */
   if(nbin > nsift) {
      if(sift_posts != NULL) free(sift_posts);
      if(sift_ngelt != NULL) free(sift_ngelt);
      sift_posts = (DATA *)calloc(nbin, sizeof(DATA));
      sift_ngelt = (int *)calloc(nbin, sizeof(int));
      nsift = nbin;
   }

/* Sample some points */
   for(i=1; i<nbin; i++) sift_posts[i] = x[(i*(n-1))/nbin];
/* and sort them */
   insort2(nbin-1, sift_posts+1, NULL);

/* Here's what we're calculating
 *    sift_posts[0]           = minimum value
 *    sift_posts[1:nbin-1]    = suggested division points amongst the values
 *    sift_posts[nbin]        = maximum value
 *    
 *    sift_ngelt[0]           = sum weight below sift_posts[1]
 *    sift_ngelt[1:nbin-1]    = sum weight in [ sift_posts[j],sift_posts[j+1] )
 *      (i.e. equal to sift_posts[j] goes into sift_ngelt[j])
 */

   sift_posts[0] = sift_posts[1];
   sift_posts[nbin] = sift_posts[nbin-1];
   for(i=0; i<=nbin; i++) sift_ngelt[i] = 0;
/* Count up all the points and find the min and max */
   for(i=0; i<n; i++) {
      if(x[i] < sift_posts[1]) {
	 if(x[i] < sift_posts[0]) sift_posts[0] = x[i];
	 sift_ngelt[0] += wgt[i];
      } else if(x[i] >= sift_posts[nbin-1]) {
	 if(x[i] >= sift_posts[nbin]) sift_posts[nbin] = x[i];
	 sift_ngelt[nbin-1] += wgt[i];
      } else {
	 for(j=1; j<nbin-1; j++) {
	    if(x[i] >= sift_posts[j] && x[i] < sift_posts[j+1]) {
	       sift_ngelt[j] += wgt[i];
	       break;
	    }
	 }
      }
   }

/* Figure out the points which can be ignored */
/* The points between sift_posts[j] and sift_posts[j+1] contain the goal */
   for(j=0, *below=0.0; j<nbin; j++) {
      if(*below + sift_ngelt[j] >= goal) break;
      *below += sift_ngelt[j];
   }

/* Swap the points between sift_posts[j] and sift_posts[k] to the beginning */
/* Ensure that the goal is not achieved at an endpoint */
   loend = hiend = lo = hi = 0;
   for(i=0, *nnew=0; i<n; i++) {
      if(x[i] >= sift_posts[j] && x[i] < sift_posts[j+1]) {
	 tmp = x[i];
	 x[i] = x[*nnew];
	 x[*nnew] = tmp;
	 if(*below + wgt[i] > goal) loend = 1;
	 if(*below + sift_ngelt[j] - wgt[i] < goal) hiend = 1;
	 if(nnew == 0) {
	    lo = hi = x[i];
	 } else {
	    if(x[i] < lo) lo = x[i];
	    if(x[i] > hi) hi = x[i];
	 }
	 *nnew += 1;
      }
   }

/* Augment the new set if loend or hiend is set */
   if(loend) {		// goal is met by the smallest entry
      j = 0;
      xhi = sift_posts[0];
      for(i=*nnew; i<n; i++) {
	 if(x[i] <= lo && x[i] >= xhi) {	// biggest but smaller than lo
	    xhi = x[i];
	    j = i;
	 }
      }
      if(j > 0) {
	 tmp = x[j];
	 x[j] = x[*nnew];
	 x[*nnew] = tmp;
	 *nnew += 1;
	 *below -= wgt[j];
      }
   }

   if(hiend) {		// goal is met by the largest entry
      j = 0;
      xhi = sift_posts[nbin];
      for(i=*nnew; i<n; i++) {
	 if(x[i] >= hi && x[i] <= xhi) {	// smallest but bigger than hi
	    xhi = x[i];
	    j = i;
	 }
      }
      if(j > 0) {
	 tmp = x[j];
	 x[j] = x[*nnew];
	 x[*nnew] = tmp;
	 *nnew += 1;
      }
   }

   return;
}

#define MEDEC_NOSIFT     80	/* don't decimate at all */
#define MEDEC_ONEPASS  1000	/* decimate twice */
#define MEDEC_TWOPASS 20000	/* decimate thrice */

/* Use median_sift() to decimate and find the median */
DATA medec(int n, DATA *buf)
{
   int nnew, nbelow, nnew2, nbelow2, nnew3, nbelow3;
   DATA med;

   if(n < MEDEC_NOSIFT) {
      qsort2(n, buf, NULL);
      med = 0.5*(buf[(n-1)/2]+buf[n/2]);
      return(med);
   }

   median_sift(n, buf, n/2, 10, &nnew, &nbelow);

   if(n < MEDEC_ONEPASS) {
      qsort2(nnew, buf, NULL);
      med = 0.5*(buf[(n-1)/2-nbelow]+buf[n/2-nbelow]);
   } else if(n < MEDEC_TWOPASS) {
      median_sift(nnew, buf, n/2-nbelow, 10, &nnew2, &nbelow2);
      qsort2(nnew2, buf, NULL);
      med = 0.5*(buf[(n-1)/2-nbelow-nbelow2]+buf[n/2-nbelow-nbelow2]);
   } else {
      median_sift(nnew, buf, n/2-nbelow, 10, &nnew2, &nbelow2);
      median_sift(nnew2, buf, n/2-nbelow-nbelow2, 10, &nnew3, &nbelow3);
      qsort2(nnew3, buf, NULL);
      med = 0.5*(buf[(n-1)/2-nbelow-nbelow2-nbelow3]+buf[n/2-nbelow-nbelow2-nbelow3]);
   }
   return(med);
}

/* Use median_sift() to decimate and find a quantile q */
DATA quantile(int n, DATA *buf, double q)
{
   int nnew, nbelow, nbelow2, goal;
   double frac;
   DATA quantile;

   quantile = 0;
   if(n <= 0) return(quantile);
   if(n == 1) return(buf[0]);

   if(q > 1) q = 1;
   if(q < 0) q = 0;

   goal = q*(n-1) - 0.49999999999;
   frac = q*(n-1) - goal;
   if(goal >= n-1) goal--;

   if(n < MEDEC_NOSIFT) {
      qsort2(n, buf, NULL);
      quantile = (1-frac)*buf[goal] + frac*buf[goal+1];
      return(quantile);
   }

   median_sift(n, buf, goal, 10, &nnew, &nbelow);

   if(n < MEDEC_ONEPASS) {
      qsort2(nnew, buf, NULL);

   } else if(n < MEDEC_TWOPASS) {
      median_sift(nnew, buf, goal-nbelow, 10, &nnew, &nbelow2);
      nbelow += nbelow2;
      qsort2(nnew, buf, NULL);

   } else {
      median_sift(nnew, buf, goal-nbelow, 10, &nnew, &nbelow2);
      nbelow += nbelow2;
      median_sift(nnew, buf, goal-nbelow, 10, &nnew, &nbelow2);
      nbelow += nbelow2;
      qsort2(nnew, buf, NULL);
   }

   quantile = (1-frac)*buf[goal-nbelow] + frac*buf[goal-nbelow+1];

   return(quantile);
}



#ifdef TEST

/* Median using Heapsort algorithm (based on Num.Rec, ex sextractor). */
DATA hmedian(DATA *ra, int n)
{
   int		l, j, ir, i;
   DATA	rra;
   if (n<2) return *ra;
   ra--;	/* Bleah, fake up fortran indexing */
   for (l = ((ir=n)>>1)+1;;) {
      if (l>1) {
	 rra = ra[--l];
      } else {
	 rra = ra[ir];
	 ra[ir] = ra[1];
	 if (--ir == 1) {
	    ra[1] = rra;
	    return n&1? ra[n/2+1] : 0.5*(ra[n/2]+ra[n/2+1]);
	 }
      }
      for (j = (i=l)<<1; j <= ir;) {
	 if (j < ir && ra[j] < ra[j+1]) ++j;
	 if (rra < ra[j]) {
	    ra[i] = ra[j];
	    j += (i=j);
	 } else {
	    j = ir + 1;
	 }
      }
      ra[i] = rra;
   }
}

int dblcmp(DATA *v1, DATA *v2)
{
   if(*v1 < *v2) return(-1);
   if(*v1 > *v2) return(1);
   return(0);
}

int main(int argc, char **argv)
{
   int i, j, k, n, niter, nasty;
   double time, time0, dummy, med, medok, medok1, medok2;
   DATA *std, *buf;
   struct timeval tv0, tv1;
   DATA median(int n, DATA *key, int *idx);

   n = 1000;
   niter = 1000;
   nasty = 0;
   if(argc > 1) sscanf(argv[1], "%d", &n);
   if(argc > 2) sscanf(argv[2], "%d", &niter);
   if(argc > 3) sscanf(argv[3], "%d", &nasty);

   std = (DATA *)calloc(n*niter, sizeof(DATA));
   for(i=0; i<n*niter; i++) std[i] = rand() / ((DATA)RAND_MAX);
   if(nasty == 3) { // nasty=3 -> reverse sorted
      for(k=0; k<niter; k++) {
	 for(j=0; j<n; j++) std[j+k*n] = -std[j+k*n];
	 qsort2(n, std+k*n, NULL);
	 for(j=0; j<n; j++) std[j+k*n] = -std[j+k*n];
      }
   } else if(nasty) {
      j = n / 7;
      for(k=0; k<niter; k++) {
	 qsort2(n, std+k*n, NULL); // nasty=1 -> sorted
	 if(nasty > 1) { // nasty=2 -> swap blocks around
	    memmove(std+k*n+0*j, std+(k+1)*n-1-1*j, j*sizeof(DATA));
	    memmove(std+k*n+1*j, std+(k+1)*n-1-2*j, j*sizeof(DATA));
	    memmove(std+k*n+2*j, std+(k+1)*n-1-3*j, j*sizeof(DATA));
	 }
      }
   }
   buf = (DATA *)calloc(n, sizeof(DATA));

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      dummy += k * buf[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("n= %5d", n);

/* Test qsort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      qsort2(n, buf, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   medok = 0.5*(buf[(n-1)/2]+buf[n/2]);
   printf("  qsort2= %8.3f %8.5f", 1e6*(time-time0)/niter, medok);
   medok1 = buf[(n-1)/2];
   medok2 = buf[n/2];

/* Test Linux qsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      qsort(buf, n, sizeof(DATA), (__compar_fn_t) dblcmp);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf[(n-1)/2]+buf[n/2]);
   printf("  qsort= %8.3f %8.5f", 1e6*(time-time0)/niter, med);

/* Test simple insertion median */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      med = median(n, buf, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  insort= %8.3f %8.5f", 1e6*(time-time0)/niter, med);

/* Test heap sort median */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      med = hmedian(buf, n);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  heap= %8.3f %8.5f", 1e6*(time-time0)/niter, med);

/* Test decimated median */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      med = medec(n, buf);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  medec= %8.3f %8.5f", 
	  1e6*(time-time0)/niter, med);

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
//      med = quantile(n, buf, 0.99);
      med = quantile(n, buf, 0.50);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %8.5f", 
	  1e6*(time-time0)/niter, med);


   printf(" [med usec]  %8.5f %8.5f\n", medok1, medok2);

   return(0);
}

#endif
