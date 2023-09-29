#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MAX_ITER 300
#define EPS 0.0001

typedef struct PointType {
    double *data;
    int length;
} Point;

typedef struct ListType {
    Point *pointsArr;
    int length;
} PointList;

void handleError(void);
double distance(Point p1, Point p2);
void addPointToList(PointList*, Point);
PointList readInput(char*);
double** createMatrixDynamically(int rows, int columns);
void printMatrix(double **a, int rows, int columns);
double** mulMatrices(double** a, double** b, int rowsA, int columnsARowsB, int columnsB);
double** subMatrices(double** a, double** b, int rows, int columns);
double** transpose(double** a, int rows, int columns);
double frobeniusNorm(double** a, int rows, int columns);
double** sym(PointList *pointList);
double** ddg(PointList *pointList);
double** norm(PointList *pointList);
double** symnmf(double** h0, PointList* pointList);
double** converge(double** h0, double** w, int rows, int columns);
int isConverged(double** h0, double** h1, int rows, int columns);
double** updateH(double** h, double** w, int rows, int columns);

