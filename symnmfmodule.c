#define PY_SSIZE_T_CLEAN
#include <Python.h>

double** matrixPytonToC(PyObject* pyMatrix){
    // convert python nested list to C double 2D array
    double **cMatrix, data;
    PyObject *row, *item;
    int rowNum, colNum, i, j;
    // get number of rows
    rowNum = PyObject_Length(pyMatrix);
    if (rowNum < 0) {
        return NULL;
    }
    // allocate row array
    cMatrix = (double **)calloc(rowNum, sizeof(double *));
    if (cMatrix == NULL) {
        printf("Memory allocation failed!\n");
        return NULL;
    }
    // get number of columns
    row = PyList_GetItem(pyMatrix, 0);
    colNum = PyObject_Length(row);
    // iterate over rows
    for (i = 0; i < rowNum; i++) {
        // allocate row
        cMatrix[i] = (double *)calloc(colNum, sizeof(double));
        if (cMatrix[i] == NULL) {
            printf("Memory allocation failed!\n");
            return NULL;
        }
        // copy data to row
        row = PyList_GetItem(pyMatrix, i);
        for (j = 0; j < colNum; j++) {
            item = PyList_GetItem(row, j);
            data = PyFloat_AsDouble(item);
            cMatrix[i][j] = data;
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

static PyObject* symnmf(PyObject* self, PyObject* args) {

}

static PyObject* sym(PyObject* self, PyObject* args) {

}

static PyObject* ddg(PyObject* self, PyObject* args) {

}

static PyObject* norm(PyObject* self, PyObject* args) {

}

static PyMethodDef symnmfMethods[] = {
    {"symnmf",
    (PyCFunction)symnmf,
    METH_VARARGS,
    PyDoc_STR("Receives initial decomposition H and normalized similarity W and optimizes H")},
    {"sym",
    (PyCFunction)sym,
    METH_VARARGS,
    PyDoc_STR("Receives Datapoints and returns similarity matrix")},
    {"ddg",
    (PyCFunction)ddg,
    METH_VARARGS,
    PyDoc_STR("Receives datapoints and returns diagonal degree matrix")},
    {"norm",
    (PyCFunction)norm,
    METH_VARARGS,
    PyDoc_STR("Receives datapoints and returns normalized similarity matrix")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfModule = {
    PyModuleDef_HEAD_INIT,
    "symnmf",
    NULL,
    -1,
    symnmfMethods
};

PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfModule);
    if (!m) {
        return NULL;
    }
    return m;
}