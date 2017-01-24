#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef double probv; /* probability type */
void allocate_all ( double ***, double ***, double ***, double **,
		    double **, double **, double **, double **,
		    double ***, double ***, int, int, int );
int query_user ( int *, char *, int *, int *,
		 char *, int *, int *, char * );
double read_fpts ( char *, double **, double *, int, int, int );
void build_k ( probv *, probv *, double *, double *, double, int, int );
void write_k ( probv **, probv **, int, int, char * );
probv evolveqc ( probv **, probv **, probv **, probv **, int, int, int );
void write_p ( probv **, int, int, char * );
void markovize ( double *, double *, double *, double, int );
void shallmfpt ( double *, double *, int, int );
void solve_mat ( double *, double *, double *, double *, int );
/*double p2mfpt ( probv **, int, int );*/
void swap ( probv **, probv **, double *, double *, int, int );
void swap_k ( probv **, probv **, int, int );
void print_equilibrium ( probv *, double *, double *, int );
void action ( probv **, probv **, probv **, probv **, double *, double *,
	      int, int, int, int, int, int, char * );





int main() {
  char fpt_name[128];
  char k_name[128];
  char p_name[128];
  probv ** p, ** q, ** kp, ** km;
  int nmlst, b, nq, n, i, j, nbins, init_mlst, ab;
  double ** fpt, max_fpt, * mfpts, * frac_p, * rp, * rm, * peq;


  /* NOTE: ab is the mode (0 = mode 1, 1 = mode 2, 2 = mode 3) */

  ab = query_user ( & nmlst,
		    fpt_name,
		    & n,
		    & b,
		    k_name,
		    & init_mlst,
		    & nq,
		    p_name );


  allocate_all ( & fpt, & kp, & km, & mfpts, & frac_p,
		 & rp, & rm, & peq, & p, & q, nmlst, n, nq );

  /* read fpts */
  max_fpt = read_fpts ( fpt_name, fpt, frac_p, nmlst, n, ab );

  nbins = ceil ( max_fpt / ((double) b) );
  printf ( "Max incubation time is %f, "
	   "so will use %d bin%s per histogram\n\n",
 	   max_fpt, nbins, nbins == 1 ? "" : "s" );

  /* construct and write histograms */
  for ( i = 0; i < nmlst; i++ ) {
    build_k ( kp[i], km[i], fpt[i], mfpts + i, b, n, nbins );
    markovize ( frac_p + i, rp + i, rm + i, mfpts[i], n );
  }
  write_k ( kp, km, nmlst, nbins, k_name );

  /* done with fpts */
  for ( i = 0; i < nmlst; i++ ) free ( fpt[i] );
  free ( fpt );

  if ( !ab ) print_equilibrium ( peq, rp, rm, nmlst );
  
  /* go! */
  action ( kp, km, p, q, frac_p, mfpts, nmlst, init_mlst,
	   nbins, nq, ab, b, p_name );

  if ( ab < 2 ) goto finish;

  printf ( "\n---- REVERSE RUN...\n\n" );

  swap ( kp, km, frac_p, mfpts, nmlst, nbins );
  /*  write_k ( kp, km, nmlst, nbins, "k_backward" );*/


  /* reset & go again */
  for ( i = 0; i < nmlst; i++ ) {
    for ( j = 0; j < nq + 2; j++ ) {
      p[i][j] = 0.0;
      q[i][j] = 0.0;
    }
  }


  action ( kp, km, p, q, frac_p, mfpts, nmlst,
	   nmlst - init_mlst + 1,
	   nbins, nq, ab, b, strcat(p_name, "_reverse") );

 finish:
  return 0;
}




void action ( probv ** kp, probv ** km, probv ** p, probv ** q,
	      double * fp, double * mfpts, int nmlst, int mi,
	      int nbins, int nq, int ab, int b, char * p_name )
{
  probv prob_end;

  mi--;

  if ( mi == nmlst - 1 ) {
    printf ( "MFPT = 0 (skipping integration)\n" );
    return;
  }

  if ( ab ) shallmfpt ( fp, mfpts, nmlst, mi );

  /* initialize */
  q[mi][1] = 1.0;
  p[mi][1] = 1.0;

  prob_end = evolveqc ( kp, km, p, q,
			nmlst - (ab ? 1 : 0),
			nbins, nq + 1 );

  printf ( "Probability remaining after integration is %f\n",
	   prob_end );

  write_p ( p, nmlst, nq, p_name );

}


/* Shalloway MFPT */
void shallmfpt ( double * fp, double * mfpts, int nmlst, int mi ) {
  int i;
  double * fm, * x, * bvec, mfpt;

  x = (double *) calloc ( nmlst - 1, sizeof(double) );
  bvec = (double *) calloc ( nmlst - 1, sizeof(double) );
  fm = (double *) calloc ( nmlst, sizeof(double) );

  for ( i = 0; i < nmlst; i++ )
    fm[i] = 1.0 - fp[i];

  bvec[mi] = 1.0;
  solve_mat ( fp, fm, x, bvec, nmlst - 1 );

  mfpt = 0.0;
  for ( i = 0; i < nmlst; i++ )
    mfpt += -mfpts[i] * x[i];

  printf ( "MFPT = %f\n", mfpt );

  free ( x );
  free ( bvec );
  free ( fm );
}





void allocate_all ( double *** fpt, double *** kp, double *** km,
		    double ** mfpts, double ** frac_p, double ** rp,
		    double ** rm, double ** peq, double *** p,
		    double *** q, int nmlst, int n, int nq )
{
  int i;

  * fpt = (double **) calloc ( nmlst, sizeof(double *) );
  * kp = (probv **) calloc ( nmlst, sizeof(probv *) );
  * km = (double **) calloc ( nmlst, sizeof(probv *) );
  * mfpts = (double *) calloc ( nmlst, sizeof(double) );
  * frac_p = (double *) calloc ( nmlst, sizeof(double) );
  * rp = (double *) calloc ( nmlst, sizeof(double) );
  * rm = (double *) calloc ( nmlst, sizeof(double) );
  * peq = (double *) calloc ( nmlst, sizeof(double) );
  * q = (probv **) calloc ( nmlst, sizeof(probv *) );
  * p = (probv **) calloc ( nmlst, sizeof(probv *) );

  for ( i = 0; i < nmlst; i++ ) {
    (*fpt)[i] = (double *) calloc ( n + 1, sizeof(double) );
    (*kp)[i] = (probv *) calloc ( n + 1, sizeof(probv) );
    (*km)[i] = (probv *) calloc ( n + 1, sizeof(probv) );
    (*p)[i] = (probv *) calloc ( nq + 2, sizeof(probv) );
    (*q)[i] = (probv *) calloc ( nq + 2, sizeof(probv) );
  }
}





void print_equilibrium ( probv * peq, double * rp,
			 double * rm, int nmlst ) {
  int i;
  double totp = 1.0;

  peq[0] = 1.0;
  for ( i = 1; i < nmlst; i++ ) {
    peq[i] = peq[i-1] * rp[i-1] / rm[i];
    totp += peq[i];
  }

  printf ( "peq: " );
  for ( i = 0; i < nmlst; i++ ) {
    peq[i] /= totp;
    printf ( "%f\t ", peq[i] );
  }
  printf ( "\n\n" );

}



void swap ( probv ** kp, probv ** km,
	    double * frac_p, double * mfpts,
	    int nmlst, int nbins )
{
  int i, j;
  double tmp;

  /* swap LFPTDs */
  swap_k ( kp, km, nmlst, nbins );

  /* swap frac_p and mfpts */
  for ( i = 0; i < ceil ( nmlst / 2 ); i++ ) {
    j = nmlst - i - 1;
    tmp = frac_p[i];
    frac_p[i] = 1.0 - frac_p[j];
    frac_p[j] = 1.0 - tmp;

    tmp = mfpts[i];
    mfpts[i] = mfpts[j];
    mfpts[j] = tmp;
  }
}



void swap_k ( probv ** kp, probv ** km, int nmlst, int nbins ) {
  int i, j, s;
  probv ** kmm, ** kpp;

  /* copy */
  kmm = (probv **) calloc ( nmlst, sizeof(probv *) );
  kpp = (probv **) calloc ( nmlst, sizeof(probv *) );
  for ( i = 0; i < nmlst; i++ ) {
    kpp[i] = (probv *) calloc ( nbins + 1, sizeof(probv) );
    kmm[i] = (probv *) calloc ( nbins + 1, sizeof(probv) );
    for ( j = 1; j <= nbins; j++ ) {
      kpp[i][j] = kp[i][j];
      kmm[i][j] = km[i][j];
    }
  }

  /* swap */
  for ( i = 0; i < nmlst; i++ ) {
    for ( j = 1; j <= nbins; j++ ) {
      s = nmlst - i - 1;
      kp[i][j] = kmm[s][j];
      km[i][j] = kpp[s][j];
    }
  }


  for ( i = 0; i < nmlst; i++ ) {
    free ( kpp[i] );
    free ( kmm[i] );
  }
  free ( kpp );
  free ( kmm );
}



/* Old way of computing MFPT */
/*double p2mfpt ( probv ** p, int nmlst, int n ) {
  int i, j;
  double mfpt;

  mfpt = 0.0;
  for ( i = 0; i < nmlst - 1; i++ )
    for ( j = 1; j <= n; j++ )
      mfpt += p[i][j];

  return mfpt;
  }*/




probv evolveqc ( probv ** kp, probv ** km, probv ** p, probv ** q,
		 int nmlst, int nbins, int nf )
{
  int t, i, j, u, z = nmlst - 1;
  probv ** cpt, totp;

  /* conditional prob of no transition. */
  cpt = (probv **) calloc ( nmlst, sizeof(probv *) );
  for ( i = 0; i < nmlst; i++ )
    cpt[i] = (probv *) calloc ( nf + 2, sizeof(probv) );


  for ( i = 0; i < nmlst; i++ ) {
    cpt[i][1] = 1.0;
    for ( j = 1; j <= nf; j++ ) {
      cpt[i][j+1] = cpt[i][j];
      if ( j <= nbins ) cpt[i][j+1] -= ( kp[i][j] + km[i][j] );
    }
  }


  /* go! */
  for ( t = 2; t <= nf; t++ ) {
    for ( i = 0; i < nmlst; i++ ) {

      u = MAX ( 1 , t - nbins );

      for ( j = u; j <= t - 1; j++ ) {
	if ( i > 0 ) q[i][t] += ( q[i-1][j] * kp[i-1][t-j] );
	if ( i < z ) q[i][t] += ( q[i+1][j] * km[i+1][t-j] );
      }

      for ( j = 1; j <= t; j++ )
	p[i][t] += ( q[i][j] * cpt[i][t+1-j] );

    }
  }

  for ( i = 0; i < nmlst; i++ ) free ( cpt[i] );
  free ( cpt );

  totp = 0.0;
  for ( i = 0; i < nmlst; i++ ) totp += p[i][nf];

  return totp;
}
 




void markovize ( double * frac_p, double * rp, double * rm,
		 double mfpt, int n )
{
  * rp = * frac_p / mfpt;
  * rm = (1.0 - * frac_p) / mfpt;
}



void build_k ( probv * kp, probv * km, double * fpt,
	       double * mfpt, double b, int n, int nbins ) {
  int j, bin;
  double fptp, tau = 0.0;

  for ( j = 1; j <= n; j++ ) {
    fptp = fabs ( fpt[j] );
    tau += fptp;
    bin = ceil ( fptp / ((double) b) );
    ( fpt[j] >= 0.0 ) ? kp[bin]++ : km[bin]++;
  }

  * mfpt = tau / (double) n;

  for ( j = 1; j <= nbins; j++ ) {
    kp[j] /= (probv) n;
    km[j] /= (probv) n;
  }

}



void write_k ( probv ** kp, probv ** km, int nmlst, int nbins, char * f ) {
  FILE * fp, * fm;
  char fname[128];
  int i, j;

  strcpy ( fname, f );
  strcat ( fname, "_p" );
  if ( (fp = fopen ( fname, "w" )) == NULL ) {
      printf ( "ERROR: can't write to file %s\n", fname );
      exit ( 1 );
  }
  strcpy ( fname, f );
  strcat ( fname, "_m" );
  if ( (fm = fopen ( fname, "w" )) == NULL ) {
      printf ( "ERROR: can't write to file %s\n", fname );
      exit ( 1 );
  }

  for ( i = 0; i < nmlst; i++ ) {
    for ( j = 1; j <= nbins; j++ ) {
      fprintf ( fp, "%f\t", kp[i][j] );
      fprintf ( fm, "%f\t", km[i][j] );
    }
    fprintf ( fp, "\n" );
    fprintf ( fm, "\n" );
  }

  fclose ( fp );
  fclose ( fm );
}




double read_fpts ( char * fpt_name, double ** fpt,
		   double * frac_p, int nmlst, int n, int ab )
{
  FILE * fp;
  char f[128], suffix[32];
  int i, j, forw_flag, back_flag;
  double max;

  max = 0.0;

  for ( i = 0; i < nmlst - (ab == 1 ? 1 : 0); i++ ) {
    forw_flag = 0;
    back_flag = 0;
    frac_p[i] = 0.0;
    strcpy ( f, fpt_name );
    sprintf ( suffix, "_%d", i + 1 );
    strcat ( f, suffix );
    if ( (fp = fopen ( f, "r" )) == NULL ) {
      printf ( "ERROR: can't open file %s\n", f );
      exit ( 1 );
    }
    printf ("Reading %s...", f );
    for ( j = 1; j <= n; j++ ) {
      if ( fscanf(fp, "%lf\n", fpt[i] + j) == EOF ) {
	printf ( "ERROR: only %d FPTs in file %s\n", j - 1, f );
	exit ( 1 );
      }
      if ( (i == 0) && (fpt[i][j] < 0) ) {
	printf ( "ERROR: first milestone leaks probability\n" );
	exit ( 1 );
      }
      if ( (i == nmlst - 1) && (fpt[i][j] > 0) ) {
	printf ( "ERROR: last milestone leaks probability\n" );
	exit ( 1 );
      }
      if ( fabs(fpt[i][j]) > max ) max = fabs ( fpt[i][j] );
      if ( fpt[i][j] > 0.0 ) {
	forw_flag = 1;
	frac_p[i]++;
      }
      else { back_flag = 1; }
    }
    if ( (i < nmlst-1) && !forw_flag ) {
      printf ( "ERROR: milestone %d has no forward probability\n", i+1 );
      exit ( 1 );
    }
    if ( (i > 0) && !back_flag && (ab == 2) ) {
      printf ( "ERROR: milestone %d has no backward probability\n", i+1 );
      exit ( 1 );
    }
    frac_p[i] /= ((double) n);
    fclose ( fp );
    printf ( "DONE\n" );
  }

  return max;
}




int query_user ( int * nmlst,
		 char * fpt_name,
		 int * n,
		 int * b,
		 char * k_name,
		 int * init_mlst,
		 int * nq,
		 char * p_name )
{
  int ab;
  char c;

  printf ( "Number of milestones: " );
  scanf ( "%d", nmlst );

  printf ( "Filename base of fpts: " );
  scanf ( "%s", fpt_name );

  printf ( "Number of trajectories per milestone: " );
  scanf ( "%d", n );

  printf ( "Desired number of steps per histogram bin: " );
  scanf ( "%d", b );

  printf ( "Desired filename base of histograms (kp/km): " );
  scanf ( "%s", k_name );

  printf ( "Initial milestone: " );
  scanf ( "%d", init_mlst );
  if ( (*init_mlst < 1)  ||  (*init_mlst > *nmlst) ) {
    printf ( "ERROR: initial milestone number must be "
	     "between 1 and %d\n", * nmlst );
    exit ( 1 );
  }

  printf ( "Desired number of integration steps: " );
  scanf ( "%d", nq );

  printf ( "Desired filename of system probability p(t): " );
  scanf ( "%s", p_name );



  /* now ask for boundary conditions */
  ab = 0;
  while ( 1 ) {
    printf ( "Absorbing final milestone? (y/n) " );

    if ( (c = getc(stdin)) == EOF ) exit(1);
    if ( c == 'n' || c == 'y' ) break;
  }
  ab += ( c == 'y' );
  if ( ab ) {
    while ( 1 ) {
      printf ( "Reverse-absorbing? (y/n) " );
      
      if ( (c = getc(stdin)) == EOF ) exit(1);
      if ( c == 'n' || c == 'y' ) break;
    }
    ab += ( c == 'y' );
  }


  printf ( "\n" );

  printf ( "%d milestones\n", * nmlst );
  printf ( "fpt filenames begin with: %s\n", fpt_name );
  printf ( "%d trajectories per milestone\n",  * n );
  printf ( "%d steps per histogram bin\n", * b );
  printf ( "Histogram filenames will begin with %s\n", k_name );
  printf ( "All probability begins in milestone %d\n",  * init_mlst );
  printf ( "%d integration steps\n",  * nq );
  printf ( "System distribution filename is %s\n", p_name );
  if ( ab ) printf ( "Absorbing final milestone\n" );
  if ( ab == 2 ) printf ( "Reverse-absorbing\n" );

  printf ( "\n" );

  return ab;
}



void write_p ( probv ** p, int nmlst, int nq, char * f ) {
  FILE * fp;
  int i, j;

  if ( (fp = fopen ( f, "w" )) == NULL ) {
      printf ( "ERROR: can't write to file %s\n", f );
      exit ( 1 );
  }
  
  for ( j = 1; j <= nq + 1; j++ ) {
    for ( i = 0; i < nmlst; i++ )
      fprintf ( fp, "%f\t", p[i][j] );
    fprintf ( fp, "\n" );
  }
    
  fclose ( fp );

}



void solve_mat ( double * fp, double * fm,
		 double * x, double * b, int n )
{
  int i;
  double * d;

  d = (double *) calloc ( n, sizeof(double) );

  d[n-1] = -1.0;

  /* make lower triangular */
  for ( i = n-2; i >= 0; i-- ) {
    d[i] = -1.0 - fp[i] * fm[i+1] / d[i+1];
    b[i] -= fm[i+1] * b[i+1] / d[i+1];
  }

  /* solve */
  x[0] = b[0] / d[0];
  for ( i = 1; i < n; i++ )
    x[i] = (b[i] - fp[i-1] * x[i-1]) / d[i];

  free ( d );
}
