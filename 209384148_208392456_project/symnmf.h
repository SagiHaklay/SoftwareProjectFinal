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

double** sym(PointList *pointList);
double** ddg(PointList *pointList);
double** norm(PointList *pointList);
double** symnmf(double** h0, double** w, int rows, int columns);