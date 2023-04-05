#include <stdio.h>
#include <math.h>

void gaussSolver(int n, double A[n][n], double b[n], double X[n]) {
    int i, j, k, l, m;
  
    //ESCALONAMENTO
    for (k = 0; k < n - 1; k++) {
        double max = fabs(A[k][k]);
        int maxIndex = k;
        
        for (i = k + 1; i < n; i++) {
            if (max < fabs(A[i][k])) {
                max = fabs(A[i][k]);
                maxIndex = i;
            }
        }
        if (maxIndex != k) {
            
            for (j = 0; j < n; j++) {
                double temp = A[k][j];
                A[k][j] = A[maxIndex][j];
                A[maxIndex][j] = temp;
            }
            double temp = b[k];
            b[k] = b[maxIndex];
            b[maxIndex] = temp;
        }

      
        if (A[k][k] == 0) {
            printf("A matriz dos coeficientes é singular\n");
            return;
        } else {
            
            for (m = k + 1; m < n; m++) {
                double F = -A[m][k] / A[k][k];
                A[m][k] = 0; //evita uma iteração
                b[m] = b[m] + F * b[k];
                for (l = k + 1; l < n; l++) {
                    A[m][l] = A[m][l] + F * A[k][l];
                }
            }
        }
    }
  
    //RESOLUÇÃO DO SISTEMA Ax=b
    for (i = n - 1; i >= 0; i--) {
        X[i] = b[i];
        for (j = i + 1; j < n; j++) {
            X[i] = X[i] - X[j] * A[i][j];
        }
        X[i] = X[i] / A[i][i];
    }
}

int main(void) {

  int n;
  int i,j;
  
  printf("*********************************\n");
  printf("* MÉTODO DE ELIMINAÇÃO DE GAUSS *\n");
  printf("*********************************\n\n");



  printf("Qual o tamanho da matriz?");
  scanf("%d", &n);
  printf("Será utilizado o método para uma matiz quadrada de lado %d", n);

  double matrizA[n][n], b[n], x[n];

  printf("\nInforme os coeficientes [A] \n");
  for(i=0; i<n; i++){
    printf("\nLinha %d\n", i+1);
    for(j=0; j<n; j++){
      printf("A[%d][%d]: ", i+1, j+1);
      scanf("%lf", &matrizA[i][j]);
    }
  }

  printf("\n***** MATRIZ A *****\n");
  for(i=0; i<n; i++){

    for(j=0; j<n; j++){
      printf("[%.2lf]", matrizA[i][j]);
    }
    printf("\n");
  }

  printf("\n***** B *****\n");
  for(i=0; i<n; i++){
    printf("B[%d]: ", i+1);
    scanf("%lf", &b[i]);
  }

  gaussSolver(n, matrizA, b, x);

  printf("\n***** RESULTADOS *****\n");
  for(i=0; i<n; i++){
    printf("X%d = %.2lf\n", i+1, x[i]);
  }
  return 0;
}