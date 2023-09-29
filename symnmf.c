#include "symnmf.h"

void handleError(void) {
    printf("An Error Has Occurred\n");
    exit(1);
}

double distance(Point p1, Point p2) {
    double sum = 0.0;
    int i;
    for (i = 0; i < p1.length; i++){
        sum += pow(p1.data[i] - p2.data[i], 2);
    }
    return sqrt(sum);
}

void addPointToList(PointList *list, Point point) {
    if (list->length == 0) {
        list->pointsArr = (Point *)malloc(sizeof(Point));
    } else {
        list->pointsArr = (Point *)realloc(list->pointsArr, (list->length + 1)*sizeof(Point));
    }
    if (list->pointsArr == NULL) {
        handleError();
    }
    list->pointsArr[list->length] = point;
    list->length++;
}

PointList readInput(char* fileName) {
    FILE *ifp = fopen( fileName,"r");
    PointList list = {NULL, 0};
    Point p = {NULL, 0};
    char line[100], *ptr, currNum[10];
    int dLength = 0, currLen = 0;
    double *data = NULL;
    while (fscanf(ifp,"%s\n", line) > 0) {
        ptr = line;
        currLen = 0;
        data = NULL;
        dLength = 0;
        while (ptr[currLen] != '\0') {
            while (ptr[currLen] != ',' && ptr[currLen] != '\0') {
                currLen++;
            }
            strncpy(currNum, ptr, currLen);
            currNum[currLen] = '\0';
            if (dLength == 0) {
                data = (double*)malloc(sizeof(double));
            } else {
                data = (double*)realloc(data, (dLength + 1)*sizeof(double));
            }
            if (data == NULL) {
                handleError();
            }
            data[dLength] = atof(currNum);
            dLength++;
            if (ptr[currLen] == ','){
                ptr += currLen + 1;
                currLen = 0;
            }

        }
        p.data = data;
        p.length = dLength;
        addPointToList(&list, p);
    }
    fclose(ifp);
    return list;
}

double** createMatrixDynamically(int rows, int columns) {
    double *p;
    double **a;
    p = calloc(rows*columns, sizeof(double));
    a = calloc(rows,sizeof(double *));
    for (int i=0 ; i<rows; i++)
        a[i] = p+i*columns;
    return a;
}

void printMatrix(double **a, int rows, int columns) {
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            printf("%.4f", a[i][j]);
            if (j != columns-1)
                printf(",");
        }
        printf("\n");
    }
}

double** mulMatrices(double** a, double** b, int rowsA, int columnsARowsB, int columnsB) {
    double **mul = createMatrixDynamically(rowsA, columnsB);
    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < columnsB; j++){
            for (int k = 0; k < columnsARowsB; k++) {
                mul[i][j] += a[i][k]*b[k][j];
            }
        }
    }
    return mul;
}

double** subMatrices(double** a, double** b, int rows, int columns) {
    double** sub = createMatrixDynamically(rows,columns);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            sub[i][j] = a[i][j] - b[i][j];
        }
    }
    return sub;
}

double** transpose(double** a, int rows, int columns) {
    double **transposed = createMatrixDynamically(columns, rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            transposed[j][i] = a[i][j];
        }
    }
    return transposed;
}

double frobeniusNorm(double** a, int rows, int columns) {
    double sum = 0.0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            sum += pow(fabs(a[i][j]),2);
            sum = sqrt(sum);
        }
    }
    return sum;
}


double** sym(PointList *pointList) {
    int n = pointList->length;
    double **a = createMatrixDynamically(n,n);
    for (int i=0 ; i<n; i++) {
        for (int j=0 ; j<n; j++) {
            if (i != j)
                a[i][j] = exp(-0.5 * pow(distance(pointList->pointsArr[i], pointList->pointsArr[j]), 2));
            else
                a[i][j] = 0.0;
        }
    }
    return a;
}

double** ddg(PointList *pointList) {
    int n = pointList->length;
    double **ddg = createMatrixDynamically(n,n);
    double **symMatrix = sym(pointList);
    double sumRow = 0.0;
    for (int i=0 ; i<n; i++) {
        for (int j=0 ; j<n; j++) {
            sumRow += symMatrix[i][j];
            if (j == n-1) {
                ddg[i][i] = sumRow;
                sumRow = 0.0;
            }
            else
                ddg[i][j] = 0.0;
        }
    }
    return ddg;
}

double** norm(PointList *pointList) {
    int n = pointList->length;
    double **ddgMatrix = ddg(pointList);
    double **symMatrix = sym(pointList);
    for (int i=0 ; i<n; i++) {
        ddgMatrix[i][i] = 1/sqrt(ddgMatrix[i][i]);
    }
    for (int i=0 ; i<n; i++) {
        for (int j=0 ; j<n; j++) {
            symMatrix[i][j] = ddgMatrix[i][i]*symMatrix[i][j]*ddgMatrix[j][j];
        }
    }
    return symMatrix;
}

double** symnmf(double** h0, PointList* pointList) { // NEED TO PUT K INSTEAD OF 2
    double** normMatrix = norm(pointList);
    double** symnmf = converge(h0, normMatrix, pointList->length, 2); //dim of h0: rows x columns. dim of w: rows x rows.
    return symnmf;
};

double** converge(double** h0, double** w, int rows, int columns) {
    int i = 0;
    double** h1;
    while (i < MAX_ITER) {
        h1 = updateH(h0, w, rows, columns);
        if (isConverged(h0, h1, rows, columns)){
            return h1;
        }
        h0 = h1;
        i++;
    }
    return h1;
}

int isConverged(double** h0, double** h1, int rows, int columns) {
    double** sub = subMatrices(h1, h0, rows, columns);
    double a = pow(frobeniusNorm(sub, rows, columns),2);
    if (a < EPS)
        return 1;
    else
        return 0;
}


double** updateH(double** h, double** w, int rows, int columns) { //dim of h: rows x columns. dim of w: rows x rows.
    double** updated = createMatrixDynamically(rows,columns);
    double** WH = mulMatrices(w, h, rows, rows, columns);
    double** Ht = transpose(h, rows, columns);
    double** HTimesHt = mulMatrices(h, Ht, rows, columns, rows);
    double** HTimesHtTimesH = mulMatrices(HTimesHt, h, rows, rows, columns);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            updated[i][j] = h[i][j]*(1-0.5+0.5*WH[i][j]/HTimesHtTimesH[i][j]);
        }
    }
    return updated;
}

int check() {
    PointList points = readInput("input_1.txt");
    PointList *pointsPtr = &points;
    double **matrix1 = sym(pointsPtr);
    double **matrix2 = ddg(pointsPtr);
    double **matrix3 = norm(pointsPtr);

    printMatrix(matrix1, points.length, points.length);
    printf("\n\n\n");
    printMatrix(matrix2, points.length, points.length);
    printf("\n\n\n");
    printMatrix(matrix3, points.length, points.length);
    printf("\n\n\n");

    return 0;
}


int main(int argc, char *argv[]) {
    PointList pointList = readInput(argv[2]);
    PointList* pointListPtr = &pointList;
    double** matrix;
    if (strcmp(argv[1], "sym")) {
        matrix = sym(pointListPtr);
    } else if (strcmp(argv[1], "ddg")) {
        matrix = ddg(pointListPtr);
    } else if (strcmp(argv[1], "norm")) {
        matrix = norm(pointListPtr);
    } else {
        handleError();
    }
    printMatrix(matrix, pointList.length, pointList.length);
    return 0;
}
