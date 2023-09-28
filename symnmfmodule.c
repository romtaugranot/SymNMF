#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "symnmf.h"

static int N = 0;
static int d = 0;

/** declarations **/

static PyObject* sym_module_imp(PyObject *self, PyObject *args);

static PyObject* ddg_module_imp(PyObject *self, PyObject *args);

static PyObject* norm_module_imp(PyObject *self, PyObject *args);

static PyObject* symnmf_module_imp(PyObject *self, PyObject *args);

static PyObject* convert_to_python_list_of_lists(double** matrix, int n, int m);

static PyObject* convert_to_python_list(double* array, int n);

static double** convert_to_c_2d_array(PyObject *list_of_lists);

static struct vector* convert_from_python_to_c(PyObject *list_of_lists);

/** implementations **/

static PyObject* sym_module_imp(PyObject *self, PyObject *args) {  
    struct vector* data_points = NULL;
    PyObject* data;
    int n;
    if(!PyArg_ParseTuple(args, "O", &data)) {
        return NULL;
    }

    data_points = convert_from_python_to_c(data);

    n = N;

    double** sym_matrix = compute_similarity_matrix(data_points, n, d);

    PyObject* returned = convert_to_python_list_of_lists(sym_matrix, n, n);

    return Py_BuildValue("O", returned);
}

static PyObject* ddg_module_imp(PyObject *self, PyObject *args)
{
    struct vector* data_points = NULL;
    PyObject* data;
    int n;
    if(!PyArg_ParseTuple(args, "O", &data)) {
        return NULL;
    }
    data_points = convert_from_python_to_c(data);
    n = N;
    double** sym_matrix = compute_similarity_matrix(data_points, n, d);
    double* ddg_list = compute_degree_matrix(sym_matrix);
    return Py_BuildValue("O", convert_to_python_list(ddg_list, n));
}

static PyObject* norm_module_imp(PyObject *self, PyObject *args)
{
    struct vector* data_points = NULL;
    PyObject* data;
    int n;
    if(!PyArg_ParseTuple(args, "O", &data)) {
        return NULL;
    }
    data_points = convert_from_python_to_c(data);
    n = N;
    double** sym_matrix = compute_similarity_matrix(data_points, n, d);
    double* ddg_list = compute_degree_matrix(sym_matrix);
    double** laplacian = compute_laplacian_matrix(sym_matrix, ddg_list);
    return Py_BuildValue("O", convert_to_python_list_of_lists(laplacian, n, n));
}

static PyObject* symnmf_module_imp(PyObject *self, PyObject *args)
{

    double** W;
    double** H;
    PyObject* W_P;
    PyObject* H_P;
    int n, k;

    if(!PyArg_ParseTuple(args, "OOii", &W_P, &H_P, &n, &k)) {
        return NULL;
    }
    W = convert_to_c_2d_array(W_P);
    H = convert_to_c_2d_array(H_P);

    H = optimize_H(W, H, n, k);
    return Py_BuildValue("O", convert_to_python_list_of_lists(H, n, k));
}

static PyMethodDef symnmfMethods[] = {
    {"sym", sym_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the similarity matrix.")},
    {"ddg", ddg_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the Diagonal Degree Matrix.")},
    {"norm", norm_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the normalized similarity matrix.")},
    {"symnmf", symnmf_module_imp, METH_VARARGS, PyDoc_STR("An implementation of symnmf.")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule", /* name of module exposed to Python */
    "SymNMF Python wrapper for C extension library", /* module documentation */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    symnmfMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&module);
    if (!m) {
        return NULL;
    }
    return m;
}


/** conversion functions **/
static PyObject* convert_to_python_list_of_lists(double** matrix, int n, int m){
    PyObject *list_of_lists;

    list_of_lists = PyList_New(n);
    
    int i,j;
    for (i = 0; i < n; i++) {
        PyList_SetItem(list_of_lists, i, PyList_New(m));
        for (j = 0; j < m; j++) {
            PyObject* python_double = Py_BuildValue("d", matrix[i][j]);
            PyList_SetItem(PyList_GetItem(list_of_lists, i), j, python_double);
        }
    }

    return list_of_lists;
}

static PyObject* convert_to_python_list(double* array, int n){
    PyObject* list;

    list = PyList_New(n);
    
    int i;
    for (i = 0; i < n; i++) {
        PyList_SetItem(list, i, Py_BuildValue("d", array[i]));
    }

    return list;
}

static double** convert_to_c_2d_array(PyObject* list_of_lists) {
    PyObject* row;
    PyObject* entry;

    int n = PyObject_Length(list_of_lists);
    int m = PyObject_Length(PyList_GetItem(list_of_lists, 0));

    double** matrix = calloc(n, sizeof(double*));

    int i, j;
    for (i = 0; i < n; i++) {
        row = PyList_GetItem(list_of_lists, i);
        matrix[i] = calloc(m, sizeof(double));
        for (j = 0; j < m; j++) {
            entry = PyList_GetItem(row, j);
            matrix[i][j] = PyFloat_AsDouble(entry);
        }
    }
    return matrix;
}

static struct vector* convert_from_python_to_c(PyObject *list_of_lists) {
    PyObject *list;
    PyObject *item;

    N = PyObject_Length(list_of_lists);
    d = PyObject_Length(PyList_GetItem(list_of_lists, 0));

    struct vector *head_vec, *curr_vec;
    struct entry *head_entry, *curr_entry;

    head_entry = malloc(sizeof(struct entry));
    if (head_entry == NULL) {   /* Memory allocation failed */
        print_error();
    }
    curr_entry = head_entry;
    curr_entry->next = NULL;

    curr_vec = malloc(sizeof(struct vector));
    if (curr_vec == NULL) {   /* Memory allocation failed */
        print_error();
    }
    curr_vec->next = NULL;
    head_vec = curr_vec;

    int i,j;
    for (i = 0; i < N; i++) {
        list = PyList_GetItem(list_of_lists, i);
        for (j = 0; j < d; j++) {
            item = PyList_GetItem(list, j);
            double num = PyFloat_AsDouble(item);
            curr_entry->value = num;
            if (j+1 < d){
                curr_entry->next = malloc(sizeof(struct entry));
                if (curr_entry == NULL) {   /* Memory allocation failed */
                    print_error();
                }
                curr_entry = curr_entry->next;
            }
            curr_entry->next = NULL;
        }

        curr_vec->entries = head_entry;
        if (i+1<N){
            curr_vec->next = calloc(1, sizeof(struct vector));
            if (curr_vec->next == NULL) {   /* Memory allocation failed */
                print_error();
            }
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
        }
        head_entry = malloc(sizeof(struct entry));
        if (head_entry == NULL) {   /* Memory allocation failed */
            print_error();
        }
        curr_entry = head_entry;
        curr_entry->next = NULL;
    }

    free_entries(head_entry);

    return head_vec;
}
