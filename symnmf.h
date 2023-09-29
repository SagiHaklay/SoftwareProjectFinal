#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MAX_ITER 300
#define EPS 0.0001

typedef struct PointType Point;
typedef struct ListType PointList;

double** sym(PointList *pointList);
double** ddg(PointList *pointList);
double** norm(PointList *pointList);
