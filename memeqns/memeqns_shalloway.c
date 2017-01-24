#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
enum modes { EQUIL, MFPT1, MFPT2 };

typedef double probv; /* probability type */
enum modes query_user ( int *, char *, int * );
double rshallmfpt ( probv *, double *, int, int );
double shallmfpt ( probv *, double *, int, int );
void solve_mat ( double *, double *, double *, double *, int );
void read_fpts ( char *, double *, probv *, int, enum modes );
void check_fpts ( probv *, int, enum modes );
void print_equilibrium (  probv *, double *, int );



int main() {
  char fpt_name[128];
  int nmlst, init_mlst;
  double * mfpts;
  probv * frac_p;
  enum modes mode;

  mode = query_user ( & nmlst, fpt_name, & init_mlst );

  mfpts = (double *) calloc ( nmlst, sizeof(double) );
  frac_p = (probv *) calloc ( nmlst, sizeof(probv) );

  read_fpts ( fpt_name, mfpts, frac_p, nmlst, mode );
  check_fpts ( frac_p, nmlst, mode );

  if ( mode == EQUIL )
    print_equilibrium ( frac_p, mfpts, nmlst );
  else {
    printf ("MFPT = %lf\n",
	    shallmfpt ( frac_p, mfpts, nmlst - 1, init_mlst - 1 ) );
    if ( mode == MFPT2 ) 
      printf ( "Reverse MFPT = %lf\n",
	       rshallmfpt ( frac_p, mfpts, nmlst, init_mlst ) );
  }

  free ( mfpts );
  free ( frac_p );
  return 0;
}





void check_fpts ( probv frac_p[], int nmlst, enum modes mode ) {
  int i;
  if ( frac_p[0] != 1 ) {
    printf ( "ERROR: first milestone leaks\n" );
    exit ( 1 );
  }
  if ( mode != MFPT1 && frac_p[nmlst-1] != 0 ) {
    printf ( "ERROR: last milestone leaks\n" );
    exit ( 1 );
  }

  for ( i = 0; i < nmlst - 1; i++ )
    if ( frac_p[i] == 0 ) {
      printf ( "ERROR: no forward flux at milestone %d\n", i+1 );
      exit ( 1 );
    }
}




double rshallmfpt ( double * frac_p, double * mfpts, int nmlst, int mi ) {
  double * rmfpts = (double *) calloc ( nmlst, sizeof(double) );
  probv * rfrac_p = (probv *) calloc ( nmlst, sizeof(probv) );
  int i;
  double mfpt;

  for ( i = 0; i < nmlst; i++ ) {
    rmfpts[i] = mfpts[nmlst-1-i];
    rfrac_p[i] = 1 - frac_p[nmlst-1-i];
  }

  mfpt = shallmfpt ( rfrac_p, rmfpts, nmlst - 1, nmlst - mi  );

  free ( rmfpts );
  free ( rfrac_p );

  return mfpt;
}




void print_equilibrium ( double * frac_p,
			 double * mfpts,
			 int nmlst )
{
  int i;
  double totp = 1.0;
  double rp; /* forward markov rate of mlst i-1 */
  double rm; /* backward markov rate of mlst i */
  double * peq;

  peq = (probv *) calloc ( nmlst, sizeof(probv) );

  peq[0] = 1.0;
  for ( i = 1; i < nmlst; i++ ) {

    rp = frac_p[i-1] / mfpts[i-1];
    rm = (1.0 - frac_p[i]) / mfpts[i];

    peq[i] = peq[i-1] * rp / rm;
    totp += peq[i];
  }


  printf ( "peq: " );
  for ( i = 0; i < nmlst; i++ ) {
    peq[i] /= totp;
    printf ( "%f\t ", peq[i] );
  }
  printf ( "\n\n" );

  free ( peq );
}



/* Shalloway MFPT */
double shallmfpt ( double * frac_p, double * mfpts, int nmlst, int mi ) {
  int i;
  double * frac_m, * x, * bvec, mfpt;

  x = (double *) calloc ( nmlst, sizeof(double) );
  bvec = (double *) calloc ( nmlst, sizeof(double) );
  frac_m = (double *) calloc ( nmlst, sizeof(double) );

  for ( i = 0; i < nmlst; i++ )
    frac_m[i] = 1.0 - frac_p[i];

  bvec[mi] = 1.0;
  solve_mat ( frac_p, frac_m, x, bvec, nmlst );

  mfpt = 0.0;
  for ( i = 0; i < nmlst; i++ )
    mfpt += -mfpts[i] * x[i];

  free ( x );
  free ( bvec );
  free ( frac_m );

  return mfpt;
}






void read_fpts ( char * fpt_name,
		 double * mfpts,
		 probv * frac_p,
		 int nmlst,
		 enum modes mode )
{
  FILE * fpl, * fp;
  char f[128];
  double fpt;
  int i, n, nm;

  if ( (fpl = fopen ( fpt_name, "r" )) == NULL ) {
    printf ( "ERROR: can't open file %s\n", fpt_name );
    exit ( 1 );
  }

  nm = mode == MFPT1 ? nmlst - 1 : nmlst;

  for ( i = 0; i < nm; i++ ) {
    fscanf ( fpl, "%s\n", f );
    if ( (fp = fopen ( f, "r" )) == NULL ) {
      printf ( "ERROR: can't open file %s\n", f );
      exit ( 1 );
    }

    printf ("Reading %s...", f );
    n = 0;
    mfpts[i] = 0.0;
    frac_p[i] = 0.0;
    while ( fscanf(fp, "%lf\n", & fpt ) != EOF ) {
      n++;
      if ( fpt > 0.0 ) frac_p[i]++;
      mfpts[i] += abs(fpt);
    }

    frac_p[i] /= ((double) n);
    mfpts[i] /= ((double) n);
    fclose ( fp );
    printf ( "DONE\n" );
  }

  fclose ( fpl );
  printf ( "\n" );
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






enum modes query_user ( int * nmlst, char * fpt_name, int * init_mlst ) {
  char c;
  enum modes mode;

  printf ( "Number of milestones: " );
  scanf ( "%d", nmlst );

  printf ( "Filename of fpt file list: " );
  scanf ( "%s", fpt_name );

  printf ( "Equilibrium run? (y/n) " );
  scanf ( "%s", & c );
  if ( c == 'y') mode = EQUIL;
  else {
    printf ( "Initial milestone: " );
    scanf ( "%d", init_mlst );
    if ( (*init_mlst < 1)  ||  (*init_mlst > *nmlst) ) {
      printf ( "ERROR: initial milestone number must be "
	       "between 1 and %d\n", * nmlst );
      exit ( 1 );
    }

    printf ( "Two-way MFPT? (y/n) " );
    scanf ( "%s", & c );
    mode = ( c == 'y' ) ? MFPT2 : MFPT1;
  }
  
  printf ( "\n" );
  switch ( mode ) {
  case EQUIL:
    printf ( "Equilibrium run\n" );
    break;
  case MFPT1:
    printf ( "MFPT run\n" );
    break;
  case MFPT2:
    printf ( "Two-way MFPT run\n" );
    break;
  default:
    printf ( "ERROR: bad mode in query_user\n" );
    exit ( 1 );
  }

  printf ( "%d milestones\n", * nmlst );
  printf ( "fpt filenames are listed in %s\n", fpt_name );
  if ( mode != EQUIL )
    printf ( "All probability begins in milestone %d\n",
	     * init_mlst );

  printf ( "\n" );

  return mode;
}
