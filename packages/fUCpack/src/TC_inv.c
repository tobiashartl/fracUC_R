#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>


SEXP TC_inverse(SEXP S0t, SEXP S0c, SEXP nu, SEXP y, SEXP n, SEXP ma, SEXP ma_c){

    /* initialize integer indices */
    int t, i, j, k, ii;

    /* copy relevant variables and data pointers */
    int N = *INTEGER(n);
    double * S0tP = REAL(S0t);
    double * S0cP = REAL(S0c);
    double * yP = REAL(y);
    double * ma_cP = REAL(ma_c);
    double * maP = REAL(ma);
    double a = REAL(nu)[0];

    /* initialize result object */
    SEXP result = PROTECT(allocMatrix(REALSXP, N+1, 2));
    double *resultP = REAL(result);

    /*first two entries for t = 1;*/
    resultP[0] = 0;
    resultP[N+1] = 0;

    /* allocate space for intermediate results */
    double * Xacc, *X, *Bacc, *Yt, *X1y;

    /* in Xacc, we iteratively store BB' + nu SS' */
    Xacc = (double *) R_alloc(N * N, sizeof(double));

    /* initialize it to zero */
    for (i = 0; i< N*N; i++) Xacc[i] = 0.0;

    /* in Xacc, we iteratively store BB' + nu SS' */
    X = (double *) R_alloc(N * N, sizeof(double));

    /* initialize it to zero */
    for (i = 0; i< N*N; i++) X[i] = 0.0;

    /* in Xacc, we iteratively store BB' + nu SS' */
    Bacc = (double *) R_alloc(N * N, sizeof(double));

    /* initialize it to zero */
    for (i = 0; i< N*N; i++) Bacc[i] = 0.0;

    /* in Xacc, we iteratively store BB' + nu SS' */
    Yt = (double *) R_alloc(N * 1, sizeof(double));

    /* initialize it to zero */
    for (i = 0; i< N; i++) Yt[i] = 0.0;

    /* in Xacc, we iteratively store BB' + nu SS' */
    X1y = (double *) R_alloc(N, sizeof(double));

    /* initialize it to zero */
    for (i = 0; i< N; i++) X1y[i] = 0.0;



    /* recursively compute X = A (Aa + Bb)' + B(Ab + Bc)' using all the information available:
     with A is the n-t+2:n block of S0t where the order of columns are reverted
     and same for B and S0c,
     a = Q[1,1], B = Q[1,2], c = Q[2,2]

     - S0t and S0c are lower triangular matrices
     - symmetry of X
     - to feed DPOSV, we only need the upper part of X
     */


    /*
     loop over t,
     NOTE: t starts at zero not 1 as in the r code
     */
    double tempA = 0.0;
    double tempB = 0.0;
    double tempYx = 0.0;
    double tempX = 0.0;
    double tempC = 0.0;


    /* int tCheck = 6; */

    for (t = 1; t < N+1; t++){

        /* compute Xacc using previous results */
        k = t;            /* dimension of the current iteration */
        /* i,j th entry of the lower block of Xacc
         that is, we have to start at m*N + m and add m whenever j/k is integer

         in order to revert the columns of the corresponding block of A and B
         instead of i'th or j'th block, we take the k-i'th or k-j'th column
         respectively
         */
        /* if (t == tCheck)Rprintf("\n------\nt = %i\n\n", t); */

        for (i = 0; i < k; i++){
            for (j = i; j < k; j++){

                tempA = 0.0;
                tempB = 0.0;

                if (i == 0){

                    /* in this case we have to compute the full row */
                    tempA += S0tP[ i ] * S0tP[ j ];

                    tempB += S0cP[ i ] * S0cP[ j ];



                    /* using triangular structure */
                    Xacc[ j * N + i] = tempA * a + tempB;
                    Bacc[ j * N + i] = tempB;


                    /* if (t<=tCheck) { */
                    /*   Rprintf("AA'[%i,%i]: %f\nIndex: %i\n", i,j, temp * a, */
                    /*           m * (N + 1) + j*N + i); */
                    /* } */

                    /* if (t==tCheck) { */
                    /*   Rprintf("Xacc[%i,%i]: %f\n", i,j, Xacc[m * (N + 1) + j * N + i]); */
                    /* } */

                } else {

                    /* in this case, we only have to add the new stuff to the old one*/

                    //if (t==2) for (ii = 0; ii < k; ii++) Rprintf("A[%i,%i]: %f", i, ii, S0tP[m * (N + 1 + (k - ii - 1)) + i]);


                    /* if (t<=tCheck) { */
                    /*   Rprintf("AApre'[%i,%i]: %f\nIndex: %i\n", i,j, Xacc[m * (N + 1) + j*N + i], */
                    /*           m * (N + 1) + j*N + i); */
                    /* } */


                    Xacc[ j * N + i] = Xacc[ (j-1) * N + i - 1] + S0tP[i] * S0tP[j] * a + S0cP[i] * S0cP[j];
                    Bacc[ j * N + i] = Bacc[ (j-1) * N + i - 1] + S0cP[i] * S0cP[j];



                    /* if (t==tCheck) { */
                    /*   Rprintf("Xacc[%i,%i]: %f\n", i,j, Xacc[m * (N + 1) + j*N + i]); */
                    /* } */

                }
            }
        }




        for (i = 0; i < k; i++){
            for (j = i; j < k; j++){
                X[j * k + i] = Xacc[j * N + i];
            }
        }




        /* Filter x_t|t and c_t|t*/
        for ( i = 0; i < k; i++){
            tempYx= 0.0;

            for (j = 0; j < k; j++){
                if(i > j){
                    tempYx += Bacc[i * N + j] * yP[k - j - 1];
                }else{
                    tempYx += Bacc[j * N + i] * yP[k - j - 1];
                }
            }
            Yt[i] = tempYx;
        }




        int nrh = 1;
        int info;

        for (ii = 0; ii < k; ii++){
            X1y[ii] = Yt[ii];
            /*X1y[ii + k] = Yc[ii];*/
        }

        F77_CALL(dposv)("U", &k, &nrh, X, &k, X1y, &k, &info);

        if (info) {
            error("Error ocurred during cholesky decomposition for x, Info = %i!", info);
        }




        /* compute the final results */
        /* compute the final results */
        tempX = 0.0;
        tempC= 0.0;

        for(i=0; i < k; i++){
            tempX += -maP[i] * X1y[i];
            tempC += -ma_cP[i] * (yP[k-i-1] - X1y[i]);

        }
        resultP[t] = tempX;
        /*resultP[t+N+1] = tempC;*/
        resultP[t+N+1] = tempC;

    }

    UNPROTECT(1);
    return result;
}









