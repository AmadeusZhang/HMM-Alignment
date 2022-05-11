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

    int len_max = len1 > len2 ? len1 : len2;
    int len_min = len1 < len2 ? len1 : len2;

    // initialize the matrix
    double M[len_max][len_max];     // match
    double I[len_max][len_max];     // insertion
    double D[len_max][len_max];     // deletion

    // initialization:
    for ( int i = 0; i < len_max; i++ ) {
        M[i][0] = I[i][0] = D[i][0] = 0;
    }
    for ( int j = 0; j < len_min; j++ ) {
        M[0][j] = I[0][j] = 0;
        D[0][j] = 1/len_min;
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
        }
    }

    // return the final score
    double finalScore = 0;
    for ( int j = 0; j < len2; j++ ) {
        finalScore += M[len1][j];
    }
     
}