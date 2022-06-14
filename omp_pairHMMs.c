/* Pairwise alignment using HMMs and accelerated by openMP
 * @author: zhzj
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <math.h>

#include <omp.h>

#define MAX_LEN 1000

int main() {
    double delta;
    delta = 0.9;

    double epsilon;
    epsilon = 0.1;

    double alpha = 1 - (2*delta);
    double beta = 1 - epsilon;
    double gamma = 1 - epsilon;

    char seq1[MAX_LEN] = "ACGTC";
    char seq2[MAX_LEN] = "ACGAA";

    int m = strlen(seq1);
    int n = strlen(seq2);

    // declare the matrix
    double M[m][n];     // match
    double I[m][n];     // insertion
    double D[m][n];     // deletion

    // initialization:
    M[0][0] = 1;
    D[0][0] = I[0][0] = 0;
    
    #pragma omp parallel for
    for (int i = 1; i < m; i++) {
        M[i][0] = 0;
        D[i][0] = I[i][0] = 0;
    }

    #pragma omp parallel for
    for ( int j = 1; j < n; j++) {
        M[0][j] = 0;
        D[0][j] = I[0][j] = 0;
    }

    double prior = 0;
    int q = 10;
    double Q = pow( 10, -q/10 );

    #pragma omp parallel
    
    /* non uso "for" perché vi è la dipendenza tra le varie antidiagonali:
     * intertask: barriera, perché i task devono essere completati prima di passare alla prossima iterazione
     * intratask: computare in parallelo le matrici
     */

    /* pseudocode:
     * compute numbers of antidiagonals
     *   for each antidiaagonal
     *     compute numbers of elements
     *      for each element
     *       compute matrix
     */

    int nAntidiagonals = m-1;
    int nElements = (m+n+1);

    for ( int k = 0; k < nElements; k++ ) {
        if ( m - 1 >= k ) {
            #pragma omp parallel for
            for ( int i = 0; i <= d; i++ ) {
                
            }
        }
    }
}