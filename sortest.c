/* test sorting routines */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include "tsort.h"

int main(int argc, char **argv)
{
   int i, j, k, n, niter, nasty, which;
   double time, time0, dummy, med;
   struct timeval tv0, tv1;
   double *wgt;
//   int *idx;
//   short int *buf;

   n = 1000;
   niter = 1000;
   nasty = 0;
   which = 1111111;
   if(argc > 1) sscanf(argv[1], "%d", &n);
   if(argc > 2) sscanf(argv[2], "%d", &niter);
   if(argc > 3) sscanf(argv[3], "%d", &nasty);
   if(argc > 4) sscanf(argv[4], "%d", &which);


   wgt = (double *)calloc(n, sizeof(double));
   for(i=0; i<n; i++) wgt[i] = 1.0;
//   idx = (int *)calloc(n, sizeof(int));
//   for(i=0; i<n; i++) idx[i] = i;

//   buf = (short int *)calloc(n, sizeof(short int));
//   for(i=0; i<n; i++) buf[i] = 32768 * (rand() / ((double)RAND_MAX) - 0.5);

//   med = quantile_s(n, buf, wgt, 0.5);

//   exit(0);



   double *std_d, *buf_d;
   std_d = (double *)calloc(n*niter, sizeof(double));
   buf_d = (double *)calloc(n, sizeof(double));

////////////////////////////////////////////////////////////////
// integer value double test
   for(i=0; i<n*niter; i++) std_d[i] = i % n;

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      dummy += k * buf_d[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("double intval n= %5d", n);

/* Test tsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      tsort_d(n, buf_d, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_d[(n-1)/2]+buf_d[n/2]);
   printf("  tsort= %8.3f %11.3f", 1e6*(time-time0)/niter, med);

/* Test insort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      insort2_d(n, buf_d, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_d[(n-1)/2]+buf_d[n/2]);
   printf("  insort= %8.3f %11.3f", 1e6*(time-time0)/niter, med);

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      med = quantile_d(n, buf_d, NULL, 0.25);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %11.3f", 1e6*(time-time0)/niter, med);

   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      med = quantile_d(n, buf_d, wgt, 0.25);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf(" wquantile= %8.3f %11.3f", 1e6*(time-time0)/niter, med);

   printf("\n");

// integer value double test
////////////////////////////////////////////////////////////////

//   exit(0);


   if( (which%10) == 1) {
////////////////////////////////////////////////////////////////
// double test
   for(i=0; i<n*niter; i++) std_d[i] = rand() / ((double)RAND_MAX);
   if(nasty == 3) { // nasty=3 -> reverse sorted
      for(k=0; k<niter; k++) {
	 for(j=0; j<n; j++) std_d[j+k*n] = -std_d[j+k*n];
	 tsort_d(n, std_d+k*n, NULL);
	 for(j=0; j<n; j++) std_d[j+k*n] = -std_d[j+k*n];
      }
   } else if(nasty) {
      j = n / 7;
      for(k=0; k<niter; k++) {
	 tsort_d(n, std_d+k*n, NULL); // nasty=1 -> sorted
	 if(nasty > 1) { // nasty=2 -> swap blocks around
	    memmove(buf_d, std_d+k*n+0*j, j*sizeof(double));
	    memmove(std_d+k*n+0*j, std_d+(k+1)*n-1-1*j, j*sizeof(double));
	    memmove(std_d+(k+1)*n-1-1*j, buf_d, j*sizeof(double));

	    memmove(buf_d, std_d+k*n+1*j, j*sizeof(double));
	    memmove(std_d+k*n+1*j, std_d+(k+1)*n-1-2*j, j*sizeof(double));
	    memmove(std_d+(k+1)*n-1-2*j, buf_d, j*sizeof(double));

	    memmove(buf_d, std_d+k*n+2*j, j*sizeof(double));
	    memmove(std_d+k*n+2*j, std_d+(k+1)*n-1-3*j, j*sizeof(double));
	    memmove(std_d+(k+1)*n-1-3*j, buf_d, j*sizeof(double));
	 }
      }
   }

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      dummy += k * buf_d[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("double        n= %5d", n);

/* Test tsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      tsort_d(n, buf_d, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_d[(n-1)/2]+buf_d[n/2]);
   printf("  tsort= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

/* Test insort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      insort2_d(n, buf_d, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_d[(n-1)/2]+buf_d[n/2]);
   printf("  insort= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      med = quantile_d(n, buf_d, NULL, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_d, std_d+k*n, n*sizeof(double));
      med = quantile_d(n, buf_d, wgt, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf(" wquantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

   printf("\n");

// double test
////////////////////////////////////////////////////////////////
   }



//   exit(0);

   if( ((which/10)%10) == 1) {
////////////////////////////////////////////////////////////////
// short int test
   short int *std_s, *buf_s;
   std_s = (short int *)calloc(n*niter, sizeof(short int));
   buf_s = (short int *)calloc(n, sizeof(short int));
   for(i=0; i<n*niter; i++) std_s[i] = 32768 * (rand() / ((double)RAND_MAX) - 0.5);
   if(nasty == 3) { // nasty=3 -> reverse sorted
      for(k=0; k<niter; k++) {
	 for(j=0; j<n; j++) std_s[j+k*n] = -std_s[j+k*n];
	 tsort_s(n, std_s+k*n, NULL);
	 for(j=0; j<n; j++) std_s[j+k*n] = -std_s[j+k*n];
      }
   } else if(nasty) {
      j = n / 7;
      for(k=0; k<niter; k++) {
	 tsort_s(n, std_s+k*n, NULL); // nasty=1 -> sorted
	 if(nasty > 1) { // nasty=2 -> swap blocks around
	    memmove(buf_s, std_s+k*n+0*j, j*sizeof(short int));
	    memmove(std_s+k*n+0*j, std_s+(k+1)*n-1-1*j, j*sizeof(short int));
	    memmove(std_s+(k+1)*n-1-1*j, buf_s, j*sizeof(short int));

	    memmove(buf_s, std_s+k*n+1*j, j*sizeof(short int));
	    memmove(std_s+k*n+1*j, std_s+(k+1)*n-1-2*j, j*sizeof(short int));
	    memmove(std_s+(k+1)*n-1-2*j, buf_s, j*sizeof(short int));

	    memmove(buf_s, std_s+k*n+2*j, j*sizeof(short int));
	    memmove(std_s+k*n+2*j, std_s+(k+1)*n-1-3*j, j*sizeof(short int));
	    memmove(std_s+(k+1)*n-1-3*j, buf_s, j*sizeof(short int));
	 }
      }
   }

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_s, std_s+k*n, n*sizeof(short int));
      dummy += k * buf_s[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("short int     n= %5d", n);

/* Test tsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_s, std_s+k*n, n*sizeof(short int));
      tsort_s(n, buf_s, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_s[(n-1)/2]+buf_s[n/2]);
   printf("  tsort= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

/* Test insort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_s, std_s+k*n, n*sizeof(short int));
      insort2_s(n, buf_s, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_s[(n-1)/2]+buf_s[n/2]);
   printf("  insort= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_s, std_s+k*n, n*sizeof(short int));
      med = quantile_s(n, buf_s, NULL, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_s, std_s+k*n, n*sizeof(short int));
      med = quantile_s(n, buf_s, wgt, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf(" wquantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

   printf("\n");
// short int test
////////////////////////////////////////////////////////////////
   }

   if( ((which/100)%10) == 1) {
////////////////////////////////////////////////////////////////
// unsigned short test
   unsigned short *std_u, *buf_u;
   std_u = (unsigned short *)calloc(n*niter, sizeof(unsigned short));
   buf_u = (unsigned short *)calloc(n, sizeof(unsigned short));
   for(i=0; i<n*niter; i++) std_u[i] = 65536.0 * rand() / ((double)RAND_MAX);
   if(nasty == 3) { // nasty=3 -> reverse sorted
      for(k=0; k<niter; k++) {
	 for(j=0; j<n; j++) std_u[j+k*n] = -std_u[j+k*n];
	 tsort_u(n, std_u+k*n, NULL);
	 for(j=0; j<n; j++) std_u[j+k*n] = -std_u[j+k*n];
      }
   } else if(nasty) {
      j = n / 7;
      for(k=0; k<niter; k++) {
	 tsort_u(n, std_u+k*n, NULL); // nasty=1 -> sorted
	 if(nasty > 1) { // nasty=2 -> swap blocks around
	    memmove(buf_u, std_u+k*n+0*j, j*sizeof(unsigned short));
	    memmove(std_u+k*n+0*j, std_u+(k+1)*n-1-1*j, j*sizeof(unsigned short));
	    memmove(std_u+(k+1)*n-1-1*j, buf_u, j*sizeof(unsigned short));

	    memmove(buf_u, std_u+k*n+1*j, j*sizeof(unsigned short));
	    memmove(std_u+k*n+1*j, std_u+(k+1)*n-1-2*j, j*sizeof(unsigned short));
	    memmove(std_u+(k+1)*n-1-2*j, buf_u, j*sizeof(unsigned short));

	    memmove(buf_u, std_u+k*n+2*j, j*sizeof(unsigned short));
	    memmove(std_u+k*n+2*j, std_u+(k+1)*n-1-3*j, j*sizeof(unsigned short));
	    memmove(std_u+(k+1)*n-1-3*j, buf_u, j*sizeof(unsigned short));
	 }
      }
   }

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_u, std_u+k*n, n*sizeof(unsigned short));
      dummy += k * buf_u[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("ushort        n= %5d", n);

/* Test tsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_u, std_u+k*n, n*sizeof(unsigned short));
      tsort_u(n, buf_u, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_u[(n-1)/2]+buf_u[n/2]);
   printf("  tsort= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

/* Test insort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_u, std_u+k*n, n*sizeof(unsigned short));
      insort2_u(n, buf_u, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_u[(n-1)/2]+buf_u[n/2]);
   printf("  insort= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_u, std_u+k*n, n*sizeof(unsigned short));
      med = quantile_u(n, buf_u, NULL, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_u, std_u+k*n, n*sizeof(unsigned short));
      med = quantile_u(n, buf_u, wgt, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf(" wquantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0);

   printf("\n");
// unsigned short test
////////////////////////////////////////////////////////////////
   }

   if( ((which/1000)%10) == 1) {
////////////////////////////////////////////////////////////////
// int test
   int *std_i, *buf_i;
   std_i = (int *)calloc(n*niter, sizeof(int));
   buf_i = (int *)calloc(n, sizeof(int));
   for(i=0; i<n*niter; i++) std_i[i] = rand() - RAND_MAX/2;
   if(nasty == 3) { // nasty=3 -> reverse sorted
      for(k=0; k<niter; k++) {
	 for(j=0; j<n; j++) std_i[j+k*n] = -std_i[j+k*n];
	 tsort_i(n, std_i+k*n, NULL);
	 for(j=0; j<n; j++) std_i[j+k*n] = -std_i[j+k*n];
      }
   } else if(nasty) {
      j = n / 7;
      for(k=0; k<niter; k++) {
	 tsort_i(n, std_i+k*n, NULL); // nasty=1 -> sorted
	 if(nasty > 1) { // nasty=2 -> swap blocks around
	    memmove(buf_i, std_i+k*n+0*j, j*sizeof(int));
	    memmove(std_i+k*n+0*j, std_i+(k+1)*n-1-1*j, j*sizeof(int));
	    memmove(std_i+(k+1)*n-1-1*j, buf_i, j*sizeof(int));

	    memmove(buf_i, std_i+k*n+1*j, j*sizeof(int));
	    memmove(std_i+k*n+1*j, std_i+(k+1)*n-1-2*j, j*sizeof(int));
	    memmove(std_i+(k+1)*n-1-2*j, buf_i, j*sizeof(int));

	    memmove(buf_i, std_i+k*n+2*j, j*sizeof(int));
	    memmove(std_i+k*n+2*j, std_i+(k+1)*n-1-3*j, j*sizeof(int));
	    memmove(std_i+(k+1)*n-1-3*j, buf_i, j*sizeof(int));
	 }
      }
   }
   buf_i = (int *)calloc(n, sizeof(int));

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_i, std_i+k*n, n*sizeof(int));
      dummy += k * buf_i[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("int           n= %5d", n);

/* Test tsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_i, std_i+k*n, n*sizeof(int));
      tsort_i(n, buf_i, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_i[(n-1)/2]+buf_i[n/2]);
   printf("  tsort= %8.3f %11.8f", 1e6*(time-time0)/niter, med/(65536.0*65536.0));

/* Test insort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_i, std_i+k*n, n*sizeof(int));
      insort2_i(n, buf_i, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_i[(n-1)/2]+buf_i[n/2]);
   printf("  insort= %8.3f %11.8f", 1e6*(time-time0)/niter, med/(65536.0*65536.0));

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_i, std_i+k*n, n*sizeof(int));
      med = quantile_i(n, buf_i, NULL, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med/(65536.0*65536.0));

   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_i, std_i+k*n, n*sizeof(int));
      med = quantile_i(n, buf_i, wgt, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf(" wquantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med/65536.0/65536.0);

   printf("\n");
// int test
////////////////////////////////////////////////////////////////
   }


   if( ((which/10000)%10) == 1) {
////////////////////////////////////////////////////////////////
// float test
   float *std_f, *buf_f;
   std_f = (float *)calloc(n*niter, sizeof(float));
   buf_f = (float *)calloc(n, sizeof(float));
   for(i=0; i<n*niter; i++) std_f[i] = rand() / ((float)RAND_MAX);
   if(nasty == 3) { // nasty=3 -> reverse sorted
      for(k=0; k<niter; k++) {
	 for(j=0; j<n; j++) std_f[j+k*n] = -std_f[j+k*n];
	 tsort_f(n, std_f+k*n, NULL);
	 for(j=0; j<n; j++) std_f[j+k*n] = -std_f[j+k*n];
      }
   } else if(nasty) {
      j = n / 7;
      for(k=0; k<niter; k++) {
	 tsort_f(n, std_f+k*n, NULL); // nasty=1 -> sorted
	 if(nasty > 1) { // nasty=2 -> swap blocks around
	    memmove(buf_f, std_f+k*n+0*j, j*sizeof(float));
	    memmove(std_f+k*n+0*j, std_f+(k+1)*n-1-1*j, j*sizeof(float));
	    memmove(std_f+(k+1)*n-1-1*j, buf_f, j*sizeof(float));

	    memmove(buf_f, std_f+k*n+1*j, j*sizeof(float));
	    memmove(std_f+k*n+1*j, std_f+(k+1)*n-1-2*j, j*sizeof(float));
	    memmove(std_f+(k+1)*n-1-2*j, buf_f, j*sizeof(float));

	    memmove(buf_f, std_f+k*n+2*j, j*sizeof(float));
	    memmove(std_f+k*n+2*j, std_f+(k+1)*n-1-3*j, j*sizeof(float));
	    memmove(std_f+(k+1)*n-1-3*j, buf_f, j*sizeof(float));
	 }
      }
   }

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_f, std_f+k*n, n*sizeof(float));
      dummy += k * buf_f[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("float         n= %5d", n);

/* Test tsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_f, std_f+k*n, n*sizeof(float));
      tsort_f(n, buf_f, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_f[(n-1)/2]+buf_f[n/2]);
   printf("  tsort= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

/* Test insort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_f, std_f+k*n, n*sizeof(float));
      insort2_f(n, buf_f, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_f[(n-1)/2]+buf_f[n/2]);
   printf("  insort= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_f, std_f+k*n, n*sizeof(float));
      med = quantile_f(n, buf_f, NULL, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_f, std_f+k*n, n*sizeof(float));
      med = quantile_f(n, buf_f, wgt, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf(" wquantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med);


   printf("\n");
// float test
////////////////////////////////////////////////////////////////
   }




   if( ((which/100000)%10) == 1) {
////////////////////////////////////////////////////////////////
// long double test
   long double *std_D, *buf_D;
   std_D = (long double *)calloc(n*niter, sizeof(long double));
   buf_D = (long double *)calloc(n, sizeof(long double));
   for(i=0; i<n*niter; i++) std_D[i] = rand() / ((double)RAND_MAX);
   if(nasty == 3) { // nasty=3 -> reverse sorted
      for(k=0; k<niter; k++) {
	 for(j=0; j<n; j++) std_D[j+k*n] = -std_D[j+k*n];
	 tsort_D(n, std_D+k*n, NULL);
	 for(j=0; j<n; j++) std_D[j+k*n] = -std_D[j+k*n];
      }
   } else if(nasty) {
      j = n / 7;
      for(k=0; k<niter; k++) {
	 tsort_D(n, std_D+k*n, NULL); // nasty=1 -> sorted
	 if(nasty > 1) { // nasty=2 -> swap blocks around
	    memmove(buf_D, std_D+k*n+0*j, j*sizeof(long double));
	    memmove(std_D+k*n+0*j, std_D+(k+1)*n-1-1*j, j*sizeof(long double));
	    memmove(std_D+(k+1)*n-1-1*j, buf_D, j*sizeof(long double));

	    memmove(buf_D, std_D+k*n+1*j, j*sizeof(long double));
	    memmove(std_D+k*n+1*j, std_D+(k+1)*n-1-2*j, j*sizeof(long double));
	    memmove(std_D+(k+1)*n-1-2*j, buf_D, j*sizeof(long double));

	    memmove(buf_D, std_D+k*n+2*j, j*sizeof(long double));
	    memmove(std_D+k*n+2*j, std_D+(k+1)*n-1-3*j, j*sizeof(long double));
	    memmove(std_D+(k+1)*n-1-3*j, buf_D, j*sizeof(long double));
	 }
      }
   }

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_D, std_D+k*n, n*sizeof(long double));
      dummy += k * buf_D[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("long double   n= %5d", n);

/* Test tsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_D, std_D+k*n, n*sizeof(long double));
      tsort_D(n, buf_D, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_D[(n-1)/2]+buf_D[n/2]);
   printf("  tsort= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

/* Test insort2 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_D, std_D+k*n, n*sizeof(long double));
      insort2_D(n, buf_D, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf_D[(n-1)/2]+buf_D[n/2]);
   printf("  insort= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

/* Test quantile */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_D, std_D+k*n, n*sizeof(long double));
      med = quantile_D(n, buf_D, NULL, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_D, std_D+k*n, n*sizeof(long double));
      med = quantile_D(n, buf_D, wgt, 0.5);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf(" wquantile= %8.3f %11.8f", 1e6*(time-time0)/niter, med);

   printf("\n");
// long double test
////////////////////////////////////////////////////////////////
   }


   if( ((which/1000000)%10) == 1) {
////////////////////////////////////////////////////////////////
// string test
   char **std_str, **buf_str;
   std_str = (char **)calloc(n*niter, sizeof(char *));
   buf_str = (char **)calloc(n, sizeof(char *));
   for(i=0; i<n*niter; i++) {
      j = 50 * (rand() / ((double)RAND_MAX));
      if(j < 3) j = 3;
      std_str[i] = malloc(j+1);
      for(k=0; k<j; k++) std_str[i][k] = (90-65) * (rand() / ((double)RAND_MAX)) + 65;
      std_str[i][j] = '\0';
   }

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf_str, std_str+k*n, n*sizeof(long double));
      dummy += k * buf_str[(k*n)/niter][2];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("string        n= %5d", n);

/* Test tsort */
#if 0
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_str, std_str+k*n, n*sizeof(long double));
      tsort_str_cpp(n, buf_str, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  tsort_C++= %8.3f `%s'", 1e6*(time-time0)/niter, buf_str[(n-1)/2]);
#endif

/* Test tsort, C version */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf_str, std_str+k*n, n*sizeof(long double));
      tsort_str(n, buf_str, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  qsort_C= %8.3f `%s'", 1e6*(time-time0)/niter, buf_str[(n-1)/2]);

   printf("\n");

// string test
////////////////////////////////////////////////////////////////
   }




   return(0);
}
