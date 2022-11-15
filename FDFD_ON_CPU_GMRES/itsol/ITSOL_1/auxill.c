#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include "LIB/globheads.h"
#include "LIB/protos.h"
#include "ios.h"


#define MAX_NUM_LEV 5        /* maximum number of levels for arms    */ 

int readhb_c(int *NN, double **AA, int **JA, int **IA, io_t *pio, 
	     double **rhs, double **guess, int *rsa)
{
    int job, ncol, nrow, nrhs, ierr;
    char guesol[3], title[73], key[9], type[4];
    int *ia = NULL, *ja = NULL, *Tia = NULL, *Tja = NULL;
    double *Ta = NULL, *a = NULL; 
    
    int n, nnz, tmp1, tmp2;


/* find out size of Harwell-Boeing matrix ---------------------------*/
    *rsa = 0;
    job = 0;
    tmp1 = tmp2 = 1;
    readmtc( &tmp1, &tmp2, &job, pio->Fname, Ta, Tja, Tia, *rhs, &nrhs,
             guesol, &nrow, &ncol, &nnz, title, key, type, &ierr );
    if( ierr != 0 ) {
        fprintf( stderr, "readhb: err in read matrix header = %d\n", ierr );
        return ierr;
    }
/* some consistency checks ------------------------------------------*/
    pio->ndim = n = ncol;
    if( nrow != ncol ) {
        fprintf( stderr, "readhb: matrix is not square\n" );
        return -1;
    }
    if( type[1] == 'S' || type[1] == 's' ) *rsa = 1;
/* allocate space ---------------------------------------------------*/
    Tia    = (int *)Malloc( sizeof(int)*(n+1), "readhb" );
    Tja    = (int *)Malloc( sizeof(int)*nnz, "readhb" );
    Ta     = (double *)Malloc( sizeof(double)*nnz, "readhb" );
    *rhs   = (double *)Malloc( sizeof(double)*n, "readhb" );
    *guess = (double *)Malloc( sizeof(double)*n, "readhb" );
/* read matrix ------------------------------------------------------*/
    job = 2;
    tmp1 = n+1;
    tmp2 = nnz;
    readmtc( &tmp1, &tmp2, &job, pio->Fname, Ta, Tja, Tia, *rhs, &nrhs,
             guesol, &nrow, &ncol, &nnz, title, key, type, &ierr );
    if( ierr != 0 ) {
        fprintf( stderr, "readhb: err in read matrix data = %d\n", ierr );
        return ierr;
    }
    tmp1 = tmp2 = 1;
    ia     = (int *)Malloc( sizeof(int)*(n+1), "readhb" );
    ja     = (int *)Malloc( sizeof(int)*nnz, "readhb" );
    a      = (double *)Malloc( sizeof(double)*nnz, "readhb" );
    csrcsc( &n, &tmp1, &tmp2, Ta, Tja, Tia, a, ja, ia); 
/*-------------------- io struct update */
    pio->ndim = n;
    pio->nnz  = nnz;
    if (*rsa == 1) 
          pio->nnz  = nnz+nnz-n;
      
    strncpy( pio->type, type, 3 );
    pio->type[3] = '\0';
/*-------------------- copy pointers */
    *AA=a; 
    *JA=ja;
    *IA = ia; 
    *NN    = n; 
/*-------------------- free temps */
    free( Ta );  Ta  = NULL;
    free( Tja ); Tja = NULL;
    free( Tia ); Tia = NULL;

    return 0;
}


int read_inputs( char *in_file, io_t *pio )
{
  FILE *finputs;
  char line[MAX_LINE], *p1, *p2;
  if( NULL == ( finputs = fopen( in_file, "r" ) ) )
      return -1;
  memset( line, 0, MAX_LINE );
  /*-------------------- line 1 : Number of params = number of tests*/
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->nparam = atoi( p1 );
  /*-------------------- line 2 : krylov subspace dimension */  
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->im = atoi( p1 );
  /*-------------------- Line 3 : maximum number of iterations */
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->maxits = atoi( p1 );
  /*-------------------- Line 4 : tolerance for stopping */  
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->tol = atof( p1 );
  
  /* this is now set manually in main 
     memset( line, 0, MAX_LINE );
     fgets( line, MAX_LINE, finputs );
     for( p1 = line; ' ' == *p1; p1++ );
     for( p2 = p1; ' ' != *p2; p2++ );
     *p2 = '\0';
     pio->eps = atof( p1 );
  */
  
/*-------------------- Line 5 : initial lfil for iluk, ilut etc */
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->lfil0 = atoi( p1 );
  /*-------------------- Line 6 : increment for lfil for consecutive tests */
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->lfilInc = atoi( p1 );
  /*-------------------- Line 7 : initial drop tol for ilut etc.. */  
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->tol0 = atof( p1 );
  /*-------------------- Line 8: multiplier for droptop for consecutive tests */
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->tolMul = atof( p1 );
/*-------------------- Line 9: fill-level used by ILUK ONLY */
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->fill_lev = atoi( p1 );
/*-------------------- Fill-lev increment is set to one */
  pio->fill_lev_inc = 1; 
/*-------------------- Line 10: Ind. set perms or PQ perms [ARMS only] */
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->perm_type = atoi( p1 );
/*-------------------- Line 11: Block size --  */
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, finputs );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  pio->Bsize = atoi( p1 );

/*-------------------- DONE -------------------- */
  fclose( finputs );
    
  return 0;
}

int get_matrix_info( FILE *fmat, io_t *pio )
{
    char line[MAX_LINE], *p1, *p2;

    memset( line, 0, MAX_LINE );
    fgets( line, MAX_LINE, fmat );
    for( p1 = line; '\'' != *p1; p1++ );
    p1++;
    for( p2 = p1; '\'' != *p2; p2++ );
    *p2 = '\0';
    strcpy( pio->Fname, p1 );
    
    for( p1 = p2+1; '\'' != *p1; p1++ );
    p1++;
    for( p2 = p1; '\'' != *p2; p2++ );
    *p2 = '\0';
    strcpy( pio->HBnameF, p1 );
    
    return 0;
}


void output_blocks( int nBlock, int *nB, FILE *f )
{
    int i;
    fprintf( f, "\nBlocks:\n" );
    for( i = 0; i < nBlock; i++ ) {
        fprintf( f, "%2d ", nB[i] );
        if( (i+1) % 25 == 0 ) fprintf( f, "\n" );
    }
    fprintf( f, "\n" );
    fflush( f );
}

void output_perm( int n, int *perm, FILE *f )
{
    int i;
    fprintf( f, "\nPermutation array:\n" );
    for( i = 0; i < n; i++ ) {
        fprintf( f, "%6d ", perm[i] );
        if( (i+1) % 10 == 0 ) fprintf( f, "\n" );
    }
    fprintf( f, "\n" );
    fflush( f );
}

int read_coo(double **VAL, int **COL, int **ROW, io_t *pio, 
	      double **rhs, double **sol) 
{
/*-------------------- reads a matrix in coordinate format. 
!  arrays VAL, COL, ROW are allocated and created 
!  for rhs: memory allocation done + artificial rhs created.
!  various other things are filled in pio  
!------------------------------------------------------------*/
  FILE *matf = NULL;
  double *aa;
  int *ii, *jj;
  int k, n, nnz;
  char *p1, *p2;
  char line[MAX_LINE];
/*-------------------- start */
  matf = fopen(pio->Fname,"r");

/*-------------------- mtx format .. in some cases n, 
                       nnz are read separately per line */
  fscanf(matf," %d %d %d\n", &n, &k, &nnz); 
  if (n != k) {
    fprintf(stdout,"This is not a square matrix -- stopping \n");
    return(1);
    } 
  /* separate reads for n and nnz 
     fscanf(matf," %d\n", &n); 
     fscanf(matf," %d\n", &nnz);  
  */
  pio->ndim = n; 
  pio->nnz  = nnz;
  pio->type[3] = '\0';
 
/*-------------------- allocate memory for matrix and rhs */
  *rhs = (double *)Malloc( n*sizeof(double), "read_coo:1" );
  *sol = (double *)Malloc( n*sizeof(double), "read_coo:2" );
  aa  = (double *)Malloc( nnz*sizeof(double),"read_coo:3" );
  jj  = (int *)Malloc( nnz*sizeof(int), "read_coo:4" );
  ii  = (int *)Malloc( nnz*sizeof(int), "read_coo:5" );

  for (k=0; k<nnz; k++) {
      fgets(line, MAX_LINE, matf);
      for( p1 = line; ' ' == *p1; p1++ );
      
/*-------------------- 1st entry - row index*/
      for( p2 = p1; ' ' != *p2; p2++ ); 
      *p2 = '\0';
      ii[k] = atoi(p1);
      
/*-------------------- 2nd entry - column index*/
      for( p1 = p2+1; ' ' == *p1; p1++ );
      for( p2 = p1; ' ' != *p2; p2++ );
      *p2 = '\0';
      jj[k]  = atoi(p1); 
      
/*-------------------- 3rd entry - nonzero entry */
      p1 = p2+1;
      aa[k] = atof(p1); 
      
  }
  *ROW = ii;
  *COL = jj;
  *VAL = aa;
  
  return 0;
}

void output_header( io_t *pio )
{
    FILE *f = pio->fout;
    fprintf( f, "\n \n"); 
    fprintf( f, " ======================================================\n" );
    fprintf( f, "   MATRIX TESTED  MATRIX %-15s  \n",pio->HBnameF);
    fprintf( f, " ------------------------------------------------------\n" ); 
    fprintf( f, "   SIZE = %-12d   NonZeros = %-12d \n",pio->ndim,pio->nnz );
    fprintf( f, "   PRECONDITIONER =   %s \n",pio->PrecMeth);
    fprintf( f, " ======================================================\n" );
    fprintf( f, "\n" );
    fprintf( f, " -------------------------------------------------------------------------\n");
    fprintf( f, "| lf  |   tol   |   P-T   |   I-T   | nz(LU)/nz | Its | ErrNorm | RsdNorm |\n");
    fprintf( f, " -------------------------------------------------------------------------\n");
    fflush( f );
}

void output_header_vb( io_t *pio ){
    FILE *f = pio->fout;
    fprintf( f, "\n \n"); 
    fprintf( f, " ======================================================\n" );
    fprintf( f, "   MATRIX TESTED  MATRIX %-15s  \n",pio->HBnameF);
    fprintf( f, " ------------------------------------------------------\n" ); 
    fprintf( f, "   SIZE = %-12d   NonZeros = %-12d \n",pio->ndim,pio->nnz );
    fprintf( f, "   PRECONDITIONER =   %s \n",pio->PrecMeth);
    fprintf( f, "|  with angle tolerance  : %-25.2f |\n", pio->eps );
  fprintf( f, "======================================================\n" );
  fprintf( f, "\n" );
  fprintf( f, "=============== PRECONDITIONER STATS =================\n" ); 
  fprintf( f, "|  Time of Hash  |  Time of Angle |  Time of Total   |\n" );
  fprintf( f, "------------------------------------------------------\n" ); 
  fprintf( f, "|  %-12.3f  |  %-13.3f |  %-13.3f   |\n",
	   pio->tm_h, pio->tm_a, pio->tm_b );
  fprintf( f, "------------------------------------------------------\n" ); 
  fprintf( f, "|  V cmpr. rate  |  E cmpr. rate  | Cmpr. Effc. (%%)  |\n" );
  fprintf( f, "------------------------------------------------------\n" ); 
  fprintf( f, "|  %-12.2f  |  %-13.2f |  %-13.2f   |\n",
	   pio->rt_v, pio->rt_e, pio->ceff );
  fprintf( f, "======================================================\n" );
  fprintf( f, "\n" );
    fprintf( f, " -------------------------------------------------------------------------\n");
    fprintf( f, "| lf  |   tol   |   P-T   |   I-T   | nz(LU)/nz | Its | ErrNorm | RsdNorm |\n");
    fprintf( f, " -------------------------------------------------------------------------\n");
    fflush( f );
}

void output_result( int lfil, io_t *pio, int iparam ){
  FILE *f = pio->fout;
  int i;
  double tol = pio->tol0;
  for( i = 1; i < iparam; i++ ) tol *= pio->tolMul;
  
  fprintf( f,"| %3d |%8.2g |%8.3f |%8.3f |%10.3f | %3d |%8.2g |%8.2g |\n", lfil, tol,
	   pio->tm_p, pio->tm_i, pio->fillfact, pio->its, pio->enorm, pio->rnorm );
  fprintf(f," -------------------------------------------------------------------------\n"); 
  fflush( f );
}


void set_arms_pars(io_t* io, int Dscale, int *ipar, 
		   double *dropcoef, int *lfil){
/*-------------------------------------------------*/
/* sets parameters required by arms preconditioner */
/* input io_t, Dscale                            */
/* output ipar tolcoef, lfil                       */
/*-------------------- trigger an error if not set */
  int j;
  for (j=0; j< 17; j++) 
    ipar[j] = 0; 
/*-------------------- */
  ipar[0]  = MAX_NUM_LEV;     /* max number of levels allowed */
  fprintf(stdout," %d maxlev \n",ipar[0]);
  ipar[1]  = io->perm_type;   /* Indset (0) / PQ (1)    permutation   */
                              /* note that these refer to completely  */
                              /* different methods for reordering A   */
                              /* 0 = standard ARMS independent sets   */
                              /* 1 = arms with ddPQ ordering          */
  ipar[2]  = io->Bsize;           /* smallest size allowed for last schur comp. */
  ipar[3]  = 1;             /* whether or not to print statistics */
/*-------------------- interlevel methods */  
  ipar[10] = 0;         /* Always do permutations - currently not used  */    
  ipar[11] = 0;         /* ILUT or ILUTP - currently only ILUT is implemented */
  ipar[12] = Dscale;    /* diagonal row scaling before PILUT 0:no 1:yes */
  ipar[13] = Dscale;    /* diagonal column scaling before PILUT 0:no 1:yes */
/*-------------------- last level methods */  
  ipar[14] = 0;         /* Always do permutations at last level */
  ipar[15] = 1;         /* -> ILUTP for last level (0 = ILUT at last level) */
  ipar[16] = Dscale;    /* diagonal row scaling  0:no 1:yes */
  ipar[17] = Dscale;    /* diagonal column scaling  0:no 1:yes */
/*-------------------- set lfil */ 
  for (j=0; j<7; j++){
    lfil[j] = io->lfil0;
  }
/*--------- dropcoef (droptol[k] = tol0*dropcoef[k]) ----- */
    dropcoef[0] = 1.6;  /* dropcoef for L of B block */
    dropcoef[1] = 1.6;   /* dropcoef for U of B block */
    dropcoef[2] = 1.6;   /* dropcoef for L\ inv F */
    dropcoef[3] = 1.6;   /* dropcoef for U\ inv E */
    dropcoef[4] = 0.004;   /* dropcoef for forming schur comple. */
    dropcoef[5] = 0.004;   /* dropcoef for last level L */
    dropcoef[6] = 0.004;   /* dropcoef for last level U */
}

void randvec (double *v, int n) {
  /* fills v with random values */
  double x;
  int k, seed = 4321;
  srand(seed);
  for (k=0; k<n; k++) {
    x = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
    v[k] = x;
  }
}




