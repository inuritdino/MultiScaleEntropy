/* file: sampen.c	Doug Lake	2 August 2002
			Last revised:	9 March 2017 (by elias.potapov@gmail.com) 1.3
-------------------------------------------------------------------------------
sampen: calculate Sample Entropy
Copyright (C) 2002-2004 Doug Lake

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.  You may also view the agreement
at http://www.fsf.org/copyleft/gpl.html.

You may contact the author via electronic mail (dlake@virginia.edu).  For
updates to this software, please visit PhysioNet (http://www.physionet.org/).

_______________________________________________________________________________

Revision history:
  1.0 (2 August 2002, Doug Lake)	Original version
  1.1 (6 January 2004, George Moody)	Removed limits on input series length
  1.2 (1 November 2004, George Moody)	Merged bug fixes from DL (normalize
					by standard deviation, detect and
					avoid divide by zero); changed code to
					use double precision, to avoid loss of
					precision for small m and large N
  1.3 (9 March, 2017, Ilya Potapov)     Made a library to use with external tools+
                                        +Multi-Scale Entropy analysis from mse.c by M. Costa


Compile this library using any standard C compiler, linking with the standard C
math library.  For example, if your compiler is gcc, use:
    gcc -shared -o libsampen.so -O -Wall -fPIC sampen.c -lm

Additional information and the executable program are available at:
    http://www.physionet.org/physiotools/sampen/.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double *readdata(char *filenm, int *filelen);
void sampen(double *y, int mm, double r, unsigned long n, double *SampEn);
void sampen2(double *y, int mm, double r, unsigned long n, double *SampEn, double *SampEnSD);
void normalize(double *data, unsigned long n);

double *readdata(char *filenm, int *filelen)
{
    FILE *ifile;
    long maxdat = 0L, npts = 0L;
    double *data = NULL, y = 0.0;

    if (strcmp(filenm, "-") == 0) {
	filenm = "standard input";
	ifile = stdin;
    }
    else if ((ifile = fopen(filenm, "rt")) == NULL) {
	fprintf(stderr, "sampen:  Could not open %s \n", filenm);
	exit(1);
    }

    while (fscanf(ifile, "%lf", &y) == 1) {
	if (++npts >= maxdat) {
	    double *s;

	    maxdat += 5000;	/* allow the input buffer to grow (the
				   increment is arbitrary) */
	    if ((s = realloc(data, maxdat * sizeof(double))) == NULL) {
		fprintf(stderr,
		"sampen: insufficient memory, truncating input at row %ld\n",
			npts);
		break;
	    }
	    data = s;
	}
	data[npts - 1] = y;
    }

    fclose(ifile);

    if (npts < 1) {
	fprintf(stderr, "sampen: %s contains no data\n", filenm);
	exit(1);
    }

    *filelen = npts;
    return (data);
}

/* This function subtracts the mean from data, then divides the data by their
   standard deviation. */
void normalize(double *data, unsigned long n)
{
    int i;
    double mean = 0;
    double var = 0;
    for (i = 0; i < n; i++)
	mean += data[i];
    mean = mean / n;
    for (i = 0; i < n; i++)
	data[i] = data[i] - mean;
    for (i = 0; i < n; i++)
	var += data[i] * data[i];
    var = sqrt(var / n);
    for (i = 0; i < n; i++)
	data[i] = data[i] / var;
}

/* sampen2 calculates an estimate of sample entropy and the variance of the
   estimate. */
void sampen2(double *y, int mm, double r, unsigned long n, double *SampEn, double *SampEnSD)
{
    double *p = NULL;
    double *v1 = NULL, *v2 = NULL, *s1 = NULL, dv;
    int *R1 = NULL, *R2 = NULL, *F2 = NULL, *F1 = NULL, *F = NULL, FF;
    int *run = NULL, *run1 = NULL;
    double *A = NULL, *B = NULL;
    double *K = NULL, *n1 = NULL, *n2 = NULL;
    int MM;
    int m, m1, i, j, nj, jj, d, d2, i1, i2, dd;
    int nm1, nm2, nm3, nm4;
    double y1;

    mm++;
    MM = 2 * mm;

    if ((run = (int *) calloc(n, sizeof(int))) == NULL)
	exit(1);
    if ((run1 = (int *) calloc(n, sizeof(int))) == NULL)
	exit(1);
    if ((R1 = (int *) calloc(n * MM, sizeof(int))) == NULL)
	exit(1);
    if ((R2 = (int *) calloc(n * MM, sizeof(int))) == NULL)
	exit(1);
    if ((F = (int *) calloc(n * MM, sizeof(int))) == NULL)
	exit(1);
    if ((F1 = (int *) calloc(n * mm, sizeof(int))) == NULL)
	exit(1);
    if ((F2 = (int *) calloc(n * mm, sizeof(int))) == NULL)
	exit(1);
    if ((K = (double *) calloc((mm + 1) * mm, sizeof(double))) == NULL)
	exit(1);
    if ((A = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((B = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((p = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((v1 = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((v2 = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((s1 = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((n1 = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((n2 = (double *) calloc(mm, sizeof(double))) == NULL)
	exit(1);

    for (i = 0; i < n - 1; i++) {
	nj = n - i - 1;
	y1 = y[i];
	for (jj = 0; jj < nj; jj++) {
	    j = jj + i + 1;
	    if (((y[j] - y1) < r) && ((y1 - y[j]) < r)) {
		run[jj] = run1[jj] + 1;
		m1 = (mm < run[jj]) ? mm : run[jj];
		for (m = 0; m < m1; m++) {
		    A[m]++;
		    if (j < n - 1)
			B[m]++;
		    F1[i + m * n]++;
		    F[i + n * m]++;
		    F[j + n * m]++;
		}
	    }
	    else
		run[jj] = 0;
	}			/* for jj */

	for (j = 0; j < MM; j++) {
	    run1[j] = run[j];
	    R1[i + n * j] = run[j];

	}
	if (nj > MM - 1)
	    for (j = MM; j < nj; j++)
		run1[j] = run[j];
    }				/* for i */

    for (i = 1; i < MM; i++)
	for (j = 0; j < i - 1; j++)
	    R2[i + n * j] = R1[i - j - 1 + n * j];
    for (i = MM; i < n; i++)
	for (j = 0; j < MM; j++)
	    R2[i + n * j] = R1[i - j - 1 + n * j];
    for (i = 0; i < n; i++)
	for (m = 0; m < mm; m++) {
	    FF = F[i + n * m];
	    F2[i + n * m] = FF - F1[i + n * m];
	    K[(mm + 1) * m] += FF * (FF - 1);
	}

    for (m = mm - 1; m > 0; m--)
	B[m] = B[m - 1];
    B[0] = (double) n *(n - 1) / 2;
    for (m = 0; m < mm; m++) {
	p[m] = (double) A[m] / B[m];
	v2[m] = p[m] * (1 - p[m]) / B[m];
    }
    dd = 1;
    for (m = 0; m < mm; m++) {
	d2 = m + 1 < mm - 1 ? m + 1 : mm - 1;
	for (d = 0; d < d2 + 1; d++) {
	    for (i1 = d + 1; i1 < n; i1++) {
		i2 = i1 - d - 1;
		nm1 = F1[i1 + n * m];
		nm3 = F1[i2 + n * m];
		nm2 = F2[i1 + n * m];
		nm4 = F2[i2 + n * m];
		for (j = 0; j < (dd - 1); j++) {
		    if (R1[i1 + n * j] >= m + 1)
			nm1--;
		    if (R2[i1 + n * j] >= m + 1)
			nm4--;
		}
		for (j = 0; j < 2 * (d + 1); j++)
		    if (R2[i1 + n * j] >= m + 1)
			nm2--;
		for (j = 0; j < (2 * d + 1); j++)
		    if (R1[i2 + n * j] >= m + 1)
			nm3--;
		K[d + 1 + (mm + 1) * m] +=
		    (double) 2 *(nm1 + nm2) * (nm3 + nm4);
	    }
	}
    }

    n1[0] = (double) n *(n - 1) * (n - 2);
    for (m = 0; m < mm - 1; m++)
	for (j = 0; j < m + 2; j++)
	    n1[m + 1] += K[j + (mm + 1) * m];
    for (m = 0; m < mm; m++) {
	for (j = 0; j < m + 1; j++)
	    n2[m] += K[j + (mm + 1) * m];
    }

    for (m = 0; m < mm; m++) {
	v1[m] = v2[m];
	dv = (n2[m] - n1[m] * p[m] * p[m]) / (B[m] * B[m]);
	if (dv > 0)
	    v1[m] += dv;
	s1[m] = (double) sqrt((double) (v1[m]));
    }

    for (m = 0; m < mm; m++) {
      if (p[m] == 0){
	    SampEn[m] = -1.;//-1 indicates Inf
	    SampEnSD[m] = -1.;
      }
      else{
	    SampEn[m] = -log(p[m]);
	    SampEnSD[m] = s1[m];
      }
    }

    free(A);
    free(B);
    free(p);
    free(run);
    free(run1);
    free(s1);
    free(K);
    free(n1);
    free(R1);
    free(R2);
    free(v1);
    free(v2);
    free(F);
    free(F1);
    free(F2);
}

/* sampen() calculates an estimate of sample entropy but does NOT calculate
   the variance of the estimate */
void sampen(double *y, int M, double r, unsigned long n, double * SampEn)
{
    double *p = NULL;
    long *run = NULL, *lastrun = NULL, N;
    double *A = NULL, *B = NULL;
    int M1, m;
    unsigned long i, j, jj, nj;
    double y1;

    M++;
    if ((run = (long *) calloc(n, sizeof(long))) == NULL)
	exit(1);
    if ((lastrun = (long *) calloc(n, sizeof(long))) == NULL)
	exit(1);
    if ((A = (double *) calloc(M, sizeof(double))) == NULL)
	exit(1);
    if ((B = (double *) calloc(M, sizeof(double))) == NULL)
	exit(1);
    if ((p = (double *) calloc(M, sizeof(double))) == NULL)
	exit(1);

    /* start running */
    for (i = 0; i < n - 1; i++) {
	nj = n - i - 1;
	y1 = y[i];
	for (jj = 0; jj < nj; jj++) {
	    j = jj + i + 1;
	    if (((y[j] - y1) < r) && ((y1 - y[j]) < r)) {
		run[jj] = lastrun[jj] + 1;
		M1 = M < run[jj] ? M : run[jj];
		for (m = 0; m < M1; m++) {
		    A[m]++;
		    if (j < n - 1)
			B[m]++;
		}
	    }
	    else
		run[jj] = 0;
	}			/* for jj */
	for (j = 0; j < nj; j++)
	    lastrun[j] = run[j];
    }				/* for i */

    N = (long) (n * (n - 1) / 2);
    p[0] = A[0] / N;

    SampEn[0] = -log(p[0]);
    for (m = 1; m < M; m++) {
	p[m] = A[m] / B[m - 1];
	if (p[m] == 0)
	    SampEn[m] = -1.;
	else
	    SampEn[m] = -log(p[m]);
    }

    free(A);
    free(B);
    free(p);
    free(run);
    free(lastrun);
}

//Another method for SampEn from mse.c file
#define M_MAX 20
double SampleEntropy(double *y, unsigned long int nlin, int m, double r, int scale)
{
    int j = scale;
    int m_max = m;
    unsigned long int i, k, l, nlin_j; 
    unsigned long int cont[M_MAX+1];
    double r_new;
    double SE = -1.0;
    
    nlin_j = (unsigned long) ((nlin/j) - m_max); 
    r_new = r;              

    for (i = 0; i <= M_MAX; i++)
	cont[i]=0;

    for (i = 0; i < nlin_j; ++i) {
	for (l = i+1; l < nlin_j; ++l) { /*self-matches are not counted*/
	    k = 0;
	    while (k < m_max && fabs(y[i+k] - y[l+k]) <= r_new)
		cont[++k]++;
	    if (k == m_max && fabs(y[i+m_max] - y[l+m_max]) <= r_new)
		cont[m_max+1]++;
	} 
    }     

    for (i = 1; i <= m_max; i++)
	if (cont[i+1] == 0 || cont[i] == 0)
	    SE = -log((double)1/((nlin_j)*(nlin_j-1)));
	else
	    SE = -log((double)cont[i+1]/cont[i]);

    if (SE < 0.0){
      printf("Error: SampEn is negative: %lf\n",SE);
      printf("Report:\n");
      printf("nlin=%lu, j=%d, nlin_j=%lu, m=%d, r=%lf\n",nlin,j,nlin_j,m,r);
      for (i = 0; i < M_MAX; i++)
	printf("cont[%lu] = %lu, cont[%lu] = %lu\n",i,cont[i],i+1,cont[i+1]);
    }
    return SE;
}

// Coarse Graining for MSE analysis
void CoarseGraining(double *y, unsigned long int nlin, int j, double *data)
{
    int i, k;

    for (i = 0; i < (unsigned long)nlin/j; i++) {
	y[i] = 0;
	for (k = 0; k < j; k++)
	    y[i] += data[i*j+k];
	y[i] /= j; 
    }
}

void MSE(double *data, unsigned long int n, int m, double r, int max_scale, double *SE)
{
  int i;
  
  //Allocate space for coarse grain time series
  double *y = (double *)malloc((size_t)n*sizeof(double));
  if(y == NULL){
    fprintf(stderr,"Error: memory allocation problem\n");
    return ;
  }
  
  //Perform normalization on the original data
  normalize(data, n);

  //Perform the loop over the scales
  for(i = 1; i <= max_scale; i++){
    CoarseGraining(y, n, i, data);
    SE[i-1] = SampleEntropy(y, n, m, r, i);
  }

  free(y);
}

void MSESD(double *data, unsigned long int n, int m, double r, int max_scale, double *SE, double *SESD)
{
  int i;
  
  //Allocate space for coarse grain time series
  double *y = (double *)malloc((size_t)n*sizeof(double));
  if(y == NULL){
    fprintf(stderr,"Error: memory allocation problem\n");
    return ;
  }
  
  //Perform normalization on the original data
  normalize(data, n);

  //Temporary variables to store SE and SESD values
  double *se_tmp = (double *)malloc((size_t)(m+1)*sizeof(double));
  double *sesd_tmp = (double *)malloc((size_t)(m+1)*sizeof(double));
  
  //Perform the loop over the scales
  for(i = 1; i <= max_scale; i++){
    CoarseGraining(y, n, i, data);
    sampen2(y, m, r, (unsigned long)n/i, se_tmp, sesd_tmp);
    SE[i-1] = se_tmp[m];
    SESD[i-1] = sesd_tmp[m];
  }

  free(y);
  free(se_tmp);
  free(sesd_tmp);
}
