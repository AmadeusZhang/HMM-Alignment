// Find the most probable path through a pair HMM given sequences x and y

#include <stdio.h>
#include <stdlib.h>

find_max( int a, int b, int c ) {
  int max;
  if ( a > b )
    max = a;
  else
    max = b;
  if ( max < c )
    max = c;
  return max;
}

int main() {

    int n, m;

    // costante: delta, tau, epsilon

    // inizializzazione delle matrici vm, vx, vy
    int vm[n][m] = 0;
    vm[0][0] = 1;

    int vx[n][m] = 0;
    int vy[n][m] = 0;
    vx[0][0] = vy[0][0] = 0;

    // definire matrice emiss_prob
    
    for ( i = 0; i < n; i++ ) {
        for ( j = 0; j < m; j++ ) {
            vm[i][j] = emiss_prob[i][j] * find_max((1-2*delta-tau)*vm[i-1][j-1], (1-epsilon-tau)*vx[i-1][j-1], (1-epsilon-tau)*vy[i-1,j-1])
            vx[i][j] = emiss_x[i] * find_max( delta*vm[i-1][j], epsilon*vx[i-1][j] );
            vy[i][j] = emiss_y[j] * find_max( delta*vm[i][j-1], epsilon*vy[i][j-1] );
        }
    }

    v_end = tau * find_max(vm[n][m], vx[n][m], vy[n][m]);

    return v_end;
}