/* Pairwise alignment using HMMs
 * @author: zhzj
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main( ) {
    // INPUT :
    // R[], H[] : read bases and haplotype bases
    // Q[] : quality scores (MANCANTE) (necessario per calcolo degli emission probabilities)
    // T[][] : transition probabilities
    
    // parameters (given by GATK HaplotypeCaller)
    double delta;           // indel start probability
    delta = 0.9;
    double epsilon;         // indel continuation probability
    epsilon = 0.1;

    // define matrix of state transistion probabilities
    double T[3][3] = {
        // poiché sono utilizzati separatamente, posso pensare di definire tre vettori anziché una matrice
        { 1-2*delta, delta, delta },
        { 1-epsilon, epsilon, 0 },
        { 1-epsilon, 0, epsilon }
    };

    // read in the sequence
    char * seq1 = NULL;
    char * seq2 = NULL;

    printf("Please input the first sequence: ");
    scanf("%s", seq1);
    printf("Please input the second sequence: ");
    scanf("%s", seq2);

    int len1 = strlen(seq1);
    int len2 = strlen(seq2);

    // initialize the matrix
    double M[len1][len2];     // match
    double I[len1][len2];     // insertion
    double D[len1][len2];     // deletion

    // initialization:
    for ( int i = 0; i < len1; i++ ) {
        M[i][0] = I[i][0] = D[i][0] = 0;
    }
    for ( int j = 0; j < len2; j++ ) {
        M[0][j] = I[0][j] = 0;
        D[0][j] = 1/len2;
    }

    // define the matrix of emission probabilities
    double lambda[len1][len2];
    double Q[len1];
    for ( int i = 0; i < len1; i++ ) {
        for ( int j = 0; j < len2; j++ ) {
            if ( seq1[i] == seq2[j] )
                // match
                lambda[i][j] = Q[i]/3;
            else
                // mismatch
                lambda[i][j] = 1-Q[i];
            
            // fill in the rest of the matrix
            M[i][j] = lambda[i][j] * ( T[0][0] * M[i-1][j-1] + T[1][0] * I[i-1][j-1] + T[2][0] * D[i-1][j-1] );
            I[i][j] = M[i-1][j] * T[0][1] + I[i-1][j] * T[1][1];
            D[i][j] = M[i][j-1] * T[0][2] + D[i][j-1] * T[2][2];

            /* TODO:
             * Le matrici M, I, D hanno dipendenza solo sugli elementi che si trovano sulla loro sinistra, alto-sinistra o alto,
             * questo implica il fatto che gli elementi sulla stessa anti-diagonale non hanno dipendenza.
             * -> parallelizzare questo calcolo
             */
        }
    }

    // return the final score
    double finalScore = 0;
    for ( int j = 0; j < len2; j++ ) {
        finalScore += M[len1][j];
    }
     
}