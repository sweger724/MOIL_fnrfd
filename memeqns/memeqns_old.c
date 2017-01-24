#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef float probv;

int query_user ( int *, char *, int *, int *,
		 char *, int *, int *, char * );
int read_fpts ( char *, int **, int, int, int );
void build_k ( int, int, int **, int, probv **, probv **, int );
void write_k ( probv **, probv **, int, int, char * );
probv evolveqc ( probv **, probv **,
		 probv **, probv **,
		 int, int, int );
void write_probv_mat ( probv **, int, int, char * );
double compute_mfpt ( probv **, int, int );



int main() {
  char fpt_name[128] /*= "fpfpt"*/;
  char k_name[128] /*"k"*/;
  char p_name[128] /*="p"*/;
  double mfpt;
  int nmlst, b, nq, n, i, ** fpt, mafpt, nbins, init_mlst, ab;
  probv ** kp, ** km, ** p, ** q, prob_end;


  ab = query_user ( & nmlst,
		    fpt_name,
		    & n,
		    & b,
		    k_name,
		    & init_mlst,
		    & nq,
		    p_name );
  
  /*  nmlst = 11;
  n = 2247;
  b = 20;
  init_mlst = 3;
  nq = 10000;
  ab = 1;

  printf ( "%d milestones\n", nmlst );
  printf ( "fpt filenames begin with: %s\n", fpt_name );
  printf ( "%d trajectories per milestone\n",  n );
  printf ( "%d steps per histogram bin\n", b );
  printf ( "Histogram filenames will begin with %s\n", k_name );
  if ( init_mlst < 1 || init_mlst > nmlst ) {
    printf ( "ERROR: initial milestone number must be "
	     "between 1 and %d\n", nmlst );
    exit ( 1 );
  }
  printf ( "Initial milestone number is %d\n", init_mlst );
  printf ( "%d integration steps\n", nq );
  printf ( "System distribution filename is %s\n", p_name ); */

  init_mlst--;


  /* allocate and read fpts */
  fpt = (int **) calloc ( nmlst, sizeof(int *) );
  for ( i = 0; i < nmlst; i++ )
    fpt[i] = (int *) calloc ( n + 1, sizeof(int) );

  mafpt = read_fpts ( fpt_name, fpt, nmlst, n, ab );
  nbins = ceil ( mafpt / ((double) b) );
  printf ( "Max incubation time is %d, "
	   "so will use %d bins per histogram\n",
	   mafpt, nbins );


  /* construct and write histograms */
  kp = (probv **) calloc ( nmlst, sizeof(probv *) );
  for ( i = 0; i < nmlst; i++ )
    kp[i] = (probv *) calloc ( nbins + 1, sizeof(probv) );

  km = (probv **) calloc ( nmlst, sizeof(probv *) );
  for ( i = 0; i < nmlst; i++ )
    km[i] = (probv *) calloc ( nbins + 1, sizeof(probv) );

  build_k ( nmlst, n, fpt, b, kp, km, nbins );
  write_k ( kp, km, nmlst, nbins, k_name );

  /* done with fpts */
  for ( i = 0; i < nmlst; i++ ) free ( fpt[i] );
  free ( fpt );


  p = (probv **) calloc ( nmlst, sizeof(probv *) );
  for ( i = 0; i < nmlst; i++ )
    p[i] = (probv *) calloc ( nq + 2, sizeof(probv) );

  q = (probv **) calloc ( nmlst, sizeof(probv *) );
  for ( i = 0; i < nmlst; i++ )
    q[i] = (probv *) calloc ( nq + 2, sizeof(probv) );
  

  /* initial condition */
  q[init_mlst][1] = 1.0;
  p[init_mlst][1] = 1.0;

  /* go! */
  prob_end = evolveqc ( q, p, kp, km, nmlst, nbins, nq );

  printf ( "End of integration. Final total probability: %f\n",
	   prob_end );

  /* write p */
  write_probv_mat ( p, nq, nmlst, p_name );
  printf ( "Wrote system prob dist'n to %s\n", p_name );


  if ( ab ) {
    mfpt = compute_mfpt ( p, nmlst, nq );
    printf ( "MFPT = %f\n", mfpt );
  }
  

  /* cleanup */
  for ( i = 0; i < nmlst; i++ ) {
    free ( kp[i] );
    free ( km[i] );
    free ( q[i] );
    free ( p[i] );
  }

  free ( kp );
  free ( km );
  free ( q );
  free ( p );

  return 0;
}





double compute_mfpt ( probv ** p, int nmlst, int n ) {
  int i, j;
  double mfpt;

  mfpt = 0.0;
  for ( i = 0; i < nmlst; i++ )
    for ( j = 1; j <= n; j++ ) mfpt += p[i][j];

  return mfpt;
}




probv evolveqc ( probv ** q, probv ** p,
		 probv ** kp, probv ** km,
		 int nmlst, int nbins, int nf )
{
  int t, i, j, u, m = nmlst - 1;
  probv ** cpt, totp;

  /* conditional prob of no transition. */
  cpt = (probv **) calloc ( nmlst, sizeof(probv *) );
  for ( i = 0; i < nmlst; i++ )
    cpt[i] = (probv *) calloc ( nf + 1, sizeof(probv) );

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
	if ( i < m ) q[i][t] += ( q[i+1][j] * km[i+1][t-j] );
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




void build_k ( int nmlst, int n, int ** fpt,
	       int b, probv ** kp, probv ** km, int nbins )
{
  int i, j, bin;

  for ( i = 0; i < nmlst; i++ )
    for ( j = 1; j <= n; j++ ) {
      bin = ceil ( (double) abs ( fpt[i][j] ) / ((double) b) );
      if ( fpt[i][j] > 0 ) {
	kp[i][bin]++;
	kp[i][0]++;
      }
      else {
	km[i][bin]++;
	km[i][0]++;
      }
    }


  /* normalize */
  for ( i = 0; i < nmlst; i++ )
    for ( j = 1; j <= nbins; j++ ) {
      kp[i][j] /= (probv) n;
      km[i][j] /= (probv) n;
    }

}



void write_k ( probv ** kp, probv ** km,
	       int nmlst, int nbins, char * f )
{
  char fname[128], suffix[2];

  /* kp */
  strcpy ( fname, f );
  sprintf ( suffix, "_p" );
  strcat ( fname, suffix );
  write_probv_mat ( kp, nbins, nmlst, fname );
  printf ( "Wrote kp histogram to %s\n", fname );
  fflush ( stdout );

  /* km */
  strcpy ( fname, f );
  sprintf ( suffix, "_m" );
  strcat ( fname, suffix );
  write_probv_mat ( km, nbins, nmlst, fname );
  printf ( "Wrote kp histogram to %s\n", fname );
  fflush ( stdout );
}




int read_fpts ( char * fpt_name, int ** fpt, int nmlst, int n, int ab ) {
  FILE * fp;
  char f[128], suffix[32];
  int i, j, max;

  max = 0;
  
  for ( i = 0; i < nmlst; i++ ) {
    strcpy ( f, fpt_name );
    sprintf ( suffix, "_%d", i + 1 );
    strcat ( f, suffix );
    if ( (fp = fopen ( f, "r" )) == NULL ) {
      printf ( "ERROR: can't open file %s\n", f );
      exit ( 1 );
    }
    for ( j = 1; j <= n; j++ ) {
      if ( fscanf(fp, "%d\n", fpt[i] + j) == EOF ) {
	printf ( "ERROR: only %d elements in file %s\n", j - 1, f );
	exit ( 1 );
      }
      if ( (i == 0) && (fpt[i][j] < 0) ) {
	printf ( "ERROR: first milestone leaks probability\n" );
	exit ( 1 );
      }
      if ( !ab && (i == nmlst - 1) && (fpt[i][j] > 0) ) {
	printf ( "ERROR: last milestone leaks probability\n" );
	exit ( 1 );
      }
      if ( abs(fpt[i][j]) > max ) max = abs ( fpt[i][j] );
    }
    fclose ( fp );
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

  while ( 1 ) {
    printf ( "Absorbing final milestone? (y/n) " );

    if ( (c = getc(stdin)) == EOF ) exit(1);
    if ( c == 'n' || c == 'y' ) break;
  }
  ab = ( c == 'y' );

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

  return ab;
}




void write_probv_mat ( probv ** m, int rows, int cols, char * fname ) {
  FILE * fp;
  int i, j;
  
  if ( (fp = fopen ( fname, "w" )) == NULL ) {
      printf ( "ERROR: can't write to file %s\n", fname );
      exit ( 1 );
  }
  for ( j = 1; j <= rows; j++ ) {
    for ( i = 0; i < cols; i++ )
      fprintf ( fp, "%f\t", m[i][j] );
    fprintf ( fp, "\n" );
  }
  fclose ( fp );
}
