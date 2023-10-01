#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

double** matrixPytonToC(PyObject* pyMatrix, int* rowNum, int* colNum, PointList* points){
    // convert python nested list to C double 2D array or PointList
    double **cMatrix = NULL, data, *rowData = NULL;
    PyObject *row, *item;
    int i, j;
    // get number of rows
    *rowNum = PyObject_Length(pyMatrix);
    if (*rowNum < 0) {
        return NULL;
    }
    // allocate row array
    if (points != NULL) {
        points->pointsArr = (Point*)calloc(*rowNum, sizeof(Point));
        if (points->pointsArr == NULL) {
            printf("Memory allocation failed!\n");
            return NULL;
        }
        points->length = *rowNum;
    } else {
        cMatrix = (double **)calloc(*rowNum, sizeof(double *));
        if (cMatrix == NULL) {
            printf("Memory allocation failed!\n");
            return NULL;
        }
    }
    
    
    // get number of columns
    row = PyList_GetItem(pyMatrix, 0);
    *colNum = PyObject_Length(row);
    // iterate over rows
    for (i = 0; i < *rowNum; i++) {
        // allocate row
        rowData = (double *)calloc(*colNum, sizeof(double));
        if (rowData == NULL) {
            printf("Memory allocation failed!\n");
            return NULL;
        }
        if (points != NULL) {
            points->pointsArr[i].data = rowData;
            points->pointsArr[i].length = *colNum;
        } else {
            cMatrix[i] = rowData;
        }
        // copy data to row
        row = PyList_GetItem(pyMatrix, i);
        for (j = 0; j < *colNum; j++) {
            item = PyList_GetItem(row, j);
            data = PyFloat_AsDouble(item);
            if (points != NULL) {
                points->pointsArr[i].data[j] = data;
            } else {
                cMatrix[i][j] = data;
            }
            
        }
    }
    return cMatrix;
}

PyObject* matrixCToPython(double** cMatrix, int rowNum, int colNum){
    // convert C double 2D array to python nested list
    PyObject *pyMatrix, *row, *item;
    int i, j;
    // initialize row list
    pyMatrix = PyList_New(rowNum);
    // iterate over rows
    for (i = 0; i < rowNum; i++) {
        // initialize row
        row = PyList_New(colNum);
        // copy data to row
        for (j = 0; j < colNum; j++) {
            item = PyFloat_FromDouble(cMatrix[i][j]);
            PyList_SetItem(row, j, item);
        }
        PyList_SetItem(pyMatrix, i, row);
    }
    return pyMatrix;
}

void freePointList(PointList list) {
    int i;
    for (i = 0; i < list.length; i++) {
        free(list.pointsArr[i].data);
    }
    free(list.pointsArr);
}
void freeMatrix(double **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++) {
        free(mat[i]);
    }
    free(mat);
}

static PyObject* symnmf_method(PyObject* self, PyObject* args) {
    PyObject *pyH, *pyW, *result;
    double **cResult, **cH, **cW;
    int rows, cols;
    if (!PyArg_ParseTuple(args, "OO", &pyH, &pyW)) {
        return NULL;
    }
    cH = matrixPytonToC(pyH, &rows, &cols, NULL);
    cW = matrixPytonToC(pyW, &rows, &rows, NULL);
    cResult = symnmf(cH, cW, rows, cols);
    result = matrixCToPython(cResult, rows, cols);
    freeMatrix(cH, rows);
    freeMatrix(cW, rows);
    //freeMatrix(cResult, rows);
    return result;
}

static PyObject* sym_method(PyObject* self, PyObject* args) {
    PointList points = {pointsArr: NULL, length: 0};
    PyObject *datapoints, *result;
    double** cResult;
    int rowNum, colNum;
    if (!PyArg_ParseTuple(args, "O", &datapoints)) {
        return NULL;
    }
    matrixPytonToC(datapoints, &rowNum, &colNum, &points);
    cResult = sym(&points);
    result = matrixCToPython(cResult, rowNum, rowNum);
    freePointList(points);
    //freeMatrix(cResult, rowNum);
    return result;
}

static PyObject* ddg_method(PyObject* self, PyObject* args) {
    PointList points = {pointsArr: NULL, length: 0};
    PyObject *datapoints, *result;
    double** cResult;
    int rowNum, colNum;
    if (!PyArg_ParseTuple(args, "O", &datapoints)) {
        return NULL;
    }
    matrixPytonToC(datapoints, &rowNum, &colNum, &points);
    cResult = ddg(&points);
    result = matrixCToPython(cResult, rowNum, rowNum);
    freePointList(points);
    //freeMatrix(cResult, rowNum);
    return result;
}

static PyObject* norm_method(PyObject* self, PyObject* args) {
    PointList points;
    PyObject *datapoints, *result;
    double** cResult;
    int rowNum, colNum;
    if (!PyArg_ParseTuple(args, "O", &datapoints)) {
        return NULL;
    }
    matrixPytonToC(datapoints, &rowNum, &colNum, &points);
    cResult = norm(&points);
    result = matrixCToPython(cResult, rowNum, rowNum);
    freePointList(points);
    //freeMatrix(cResult, rowNum);
    return result;
}

static PyMethodDef symnmfMethods[] = {
    {"symnmf",
    (PyCFunction)symnmf_method,
    METH_VARARGS,
    PyDoc_STR("Receives initial decomposition H and normalized similarity W and optimizes H")},
    {"sym",
    (PyCFunction)sym_method,
    METH_VARARGS,
    PyDoc_STR("Receives Datapoints and returns similarity matrix")},
    {"ddg",
    (PyCFunction)ddg_method,
    METH_VARARGS,
    PyDoc_STR("Receives datapoints and returns diagonal degree matrix")},
    {"norm",
    (PyCFunction)norm_method,
    METH_VARARGS,
    PyDoc_STR("Receives datapoints and returns normalized similarity matrix")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfModule = {
    PyModuleDef_HEAD_INIT,
    "symNMF",
    NULL,
    -1,
    symnmfMethods
};

PyMODINIT_FUNC PyInit_symNMF(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfModule);
    if (!m) {
        return NULL;
    }
    return m;
}