#include "mpi.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define tamanho 2048
#define iteracoes 2000

int getNeighbors(float **grid, int i, int j) {
  int neighbors = 0;
  if (grid[((i - 1) + tamanho) % tamanho][j] == 1.0) // vizinho superior
    neighbors++;
  if (grid[((i - 1) + tamanho) % tamanho][((j - 1) + tamanho) % tamanho] ==
      1.0) // vizinho superior esquerdo
    neighbors++;
  if (grid[((i - 1) + tamanho) % tamanho][(j + 1) % tamanho] ==
      1.0) // vizinho superior direito
    neighbors++;
  if (grid[i][((j - 1) + tamanho) % tamanho] == 1.0) // vizinho esquerda
    neighbors++;
  if (grid[i][(j + 1) % tamanho] == 1.0) // viziho direito
    neighbors++;
  if (grid[(i + 1) % tamanho][((j - 1) + tamanho) % tamanho] ==
      1.0) // vizinho inferior esquerdo
    neighbors++;
  if (grid[(i + 1) % tamanho][j] == 1.0) // vizinho inferior
    neighbors++;
  if (grid[(i + 1) % tamanho][(j + 1) % tamanho] ==
      1.0) // vizinho inferior direito
    neighbors++;
  return neighbors;
}

int main(int argc, char *argv[]) {

  int rank, size;
  int nameSize;
  int start = 0, end = 0;

  float media, soma;

  float **grid;
  float **newGrid;

  float recvBuf1[tamanho], recvBuf2[tamanho], sendBuf1[tamanho],
      sendBuf2[tamanho];
  int anteProc, proxProc;

  struct timeval inicio, final;
  long long tmili;

  MPI_Status status;

  gettimeofday(&inicio, NULL);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int tamanhoLocal = (int)ceil(((float)tamanho / (float)size));
  int cel_vivas = 0;

  start = rank * tamanhoLocal;
  tamanhoLocal = start + tamanhoLocal < tamanho ? tamanhoLocal : tamanho - start + 1;

  anteProc = ((rank + size - 1) % size);
  proxProc = ((rank + 1) % size);

  grid = (float **)malloc(sizeof(float *) * (tamanhoLocal + 2));
  newGrid = (float **)malloc(sizeof(float *) * (tamanhoLocal + 2));

  for (int i = 0; i < tamanhoLocal + 2; i++) {
    grid[i] = (float *)malloc(sizeof(float) * tamanho + 2);
    newGrid[i] = (float *)malloc(sizeof(float) * tamanho + 2);
  }

  if (rank == 0) {
    int lin = 1;
    int col = 1;
    // GLIDER
    grid[lin][col + 1] = 1.0;
    grid[lin + 1][col + 2] = 1.0;
    grid[lin + 2][col] = 1.0;
    grid[lin + 2][col + 1] = 1.0;
    grid[lin + 2][col + 2] = 1.0;

    lin = 10;
    col = 30;
    // R-pentomino
    grid[lin][col + 1] = 1.0;
    grid[lin][col + 2] = 1.0;
    grid[lin + 1][col] = 1.0;
    grid[lin + 1][col + 1] = 1.0;
    grid[lin + 2][col + 1] = 1.0;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int k = 0; k < iteracoes; k++) {
    for (int j = 0; j < tamanho; j++) {
      sendBuf1[j] = grid[1][j + 1];
      sendBuf2[j] = grid[tamanhoLocal][j + 1];
    }

    MPI_Sendrecv(sendBuf2, tamanho, MPI_INTEGER, proxProc, 1, recvBuf1, tamanho,
                 MPI_INTEGER, anteProc, 1, MPI_COMM_WORLD, &status);

    MPI_Sendrecv(sendBuf1, tamanho, MPI_INTEGER, anteProc, 1, recvBuf2, tamanho,
                 MPI_INTEGER, proxProc, 1, MPI_COMM_WORLD, &status);

    for (int j = 0; j < tamanho; j++) {
      grid[0][j + 1] = recvBuf1[j];
      grid[tamanhoLocal + 1][j + 1] = recvBuf2[j];
    }

    for (int i = 0; i < tamanhoLocal + 2; i++) {
      grid[i][tamanho + 1] = grid[i][1];
      grid[i][0] = grid[i][tamanho];
    }

    for (int i = 1; i < tamanhoLocal + 1; i++) {
      for (int j = 1; j < tamanho + 1; j++) {
        soma = getNeighbors(grid, i, j);
        if (grid[i][j] == 1 && (soma == 2 || soma == 3)) {
          media = soma / 8.0;
          newGrid[i][j] = ceilf(media);
        } else if (grid[i][j] == 0.0 && soma == 3) {
          media = soma / 8.0;
          newGrid[i][j] = ceilf(media);
        } else
          newGrid[i][j] = 0.0;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cel_vivas = 0;
    for (int i = 1; i < tamanhoLocal + 1; i++) {
      for (int j = 1; j < tamanho + 1; j++) {
        grid[i][j] = newGrid[i][j];
        cel_vivas += newGrid[i][j];
      }
    }
    if(rank==0 && k<5){
      int cel_vivas50 = 0;
      int i;
      int j;
      for (i = 0; i < 50; i++) {
        for (j = 0; j < 50; j++) {
          printf("%.2f ", grid[i][j]);
          if (grid[i][j] == 1.0) {
           cel_vivas50++;
          }
        }
      }
      printf("\nGeração %d: %d\n\n", k + 1, cel_vivas50);
    }
  }

  if (rank == 0) {
    int results[size], sum = 0;
    results[0] = cel_vivas;

    for (int i = 1; i < size; i++) {
      MPI_Recv(&results[i], 1, MPI_INTEGER, i, 10, MPI_COMM_WORLD, &status);
    }

    for (int i = 0; i < size; i++) {
      sum += results[i];
    }
    printf("\nTotal: %d\n", sum);
  } else {
    MPI_Send(&cel_vivas, 1, MPI_INTEGER, 0, 10, MPI_COMM_WORLD);
  }

  MPI_Finalize();

  gettimeofday(&final, NULL);
  tmili = (int)(1000 * (final.tv_sec - inicio.tv_sec) +
                (final.tv_usec - inicio.tv_usec) / 1000);
  printf("Tempo decorrido: %lld ms\n", tmili);

  return 0;
}
