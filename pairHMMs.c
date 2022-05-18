/* Pairwise alignment using HMMs
 * @author: zhzj
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <math.h>

#define MAX_LEN 1000

int main( ) {
    // INPUT :
    // R[], H[] : read bases and haplotype bases
    // Q[] : quality scores
    // T[][] : transition probabilities
    
    // parameters (given by GATK HaplotypeCaller)
    double delta;           // indel start probability
    delta = 0.9;
    double epsilon;         // indel continuation probability
    epsilon = 0.1;

    // define matrix of state transistion probabilities
    /*
    double T[3][3] = {
        { 1-(2*delta), delta, delta },
        { 1-epsilon, epsilon, 0 },
        { 1-epsilon, 0, epsilon }
    };
    */

    // poiché sono utilizzati separatamente, posso pensare di definire tre vettori anziché una matrice
    double alpha = 1 - (2*delta);
    double beta = 1 - epsilon;
    double gamma = 1 - epsilon;

    // read in the sequence
    char seq1[MAX_LEN] = "ACGTC";
    char seq2[MAX_LEN] = "ACGAA";

    // printf("Please input the first sequence: ");
    // scanf("%s", seq1);
    // printf("Please input the second sequence: ");
    // scanf("%s", seq2);

    int len1 = strlen(seq1);
    int len2 = strlen(seq2);

    // declare the matrix
    double M[len1][len2];     // match
    double I[len1][len2];     // insertion
    double D[len1][len2];     // deletion

    // initialization:
    M[0][0] = 1;
    D[0][0] = I[0][0] = 0;
    for ( int i = 1; i < len1; i++ ) {
        M[i][0] = I[i][0] = D[i][0] = 0;
    }
    
    for ( int j = 1; j < len2; j++ ) {
        M[0][j] = I[0][j] = 0;
        D[0][j] = 1/len2;               // To allow the alignment of the haplotype to start anywhere on the read
                                        // without penalty, we need to initialize the entire first row of the deletion
                                        // matrix with the normalized factor 1/len2
    }

    // define the matrix of emission probabilities
    double prior = 0;                   // probability of emitting an aligned pair of symbols
    int q = 10;                         // constant default phred-scaled indel start quality
    double Q = pow( 10, -q/10 );        // quality score
    
    for ( int i = 1; i < len1; i++ ) {
        for ( int j = 1; j < len2; j++ ) {
            if ( seq1[i] == seq2[j] )
                // match
                prior = 1 - Q;
            else
                // mismatch
                prior = Q / 3;
            

            // M[i][j] = prior * ( T[0][0] * M[i-1][j-1] + T[1][0] * I[i-1][j-1] + T[2][0] * D[i-1][j-1] );
            // I[i][j] = T[0][1] * M[i-1][j] + T[1][1] * I[i-1][j];
            // D[i][j] = T[0][2] * M[i][j-1] + T[2][2] * D[i][j-1];

            M[i][j] = prior * ( alpha * M[i-1][j-1] + beta * I[i-1][j-1] + gamma * D[i-1][j-1] );
            I[i][j] = delta * M[i-1][j] + epsilon * I[i-1][j];
            D[i][j] = delta * M[i][j-1] + epsilon * D[i][j-1];

            /* TODO:
             * Le matrici M, I, D hanno dipendenza solo sugli elementi che si trovano sulla loro sinistra, alto-sinistra o alto,
             * questo implica il fatto che gli elementi sulla stessa anti-diagonale non hanno dipendenza.
             * -> parallelizzare questo calcolo
             */

        }
    }

    // return the final score
    double finalScore = 0.0;

    // for ( int j = 0; j < len2; j++ ) 
    //    finalScore += ( M[len1-1][j] + I[len1-1][j] );

    finalScore = M[len1-1][len2-1] + I[len1-1][len2-1];

    printf("Final score: %f\n", finalScore);

}