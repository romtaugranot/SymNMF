#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "symnmf.h"

static PyObject* sym_module_imp(PyObject *self, PyObject *args) {  
    struct vector* data_points;
    char* file_path;
    if(!PyArg_ParseTuple(args, "s", &file_path)) {
        return NULL;
    }
    printf("path is: %s\n", file_path);
    data_points = read_data_points(file_path);
    

    return Py_BuildValue("O", data_points);
}

static PyObject* ddg_module_imp(PyObject *self, PyObject *args)
{

    if(!PyArg_ParseTuple(args, "", NULL)) {
        return NULL;
    }

    return Py_BuildValue("O", NULL);
}

static PyObject* norm_module_imp(PyObject *self, PyObject *args)
{

    if(!PyArg_ParseTuple(args, "", NULL)) {
        return NULL;
    }

    return Py_BuildValue("O", NULL);
}

static PyObject* symnmf_module_imp(PyObject *self, PyObject *args)
{

    if(!PyArg_ParseTuple(args, "", NULL)) {
        return NULL;
    }

    return Py_BuildValue("O", NULL);
}

static PyObject* get_list(PyObject* self, PyObject* args)
{
    int N,r;
    PyObject* python_val;
    PyObject* python_int;
    if (!PyArg_ParseTuple(args, "i", &N)) {
        return NULL;
    }
    python_val = PyList_New(N);
    for (int i = 0; i < N; ++i)
    {
        r = i;
        python_int = Py_BuildValue("i", r);
        PyList_SetItem(python_val, i, python_int);
    }
    return python_val;
}

static PyMethodDef symnmfMethods[] = {
    {"sym", sym_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the similarity matrix.")},
    {"ddg", ddg_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the Diagonal Degree Matrix.")},
    {"norm", norm_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the normalized similarity matrix.")},
    {"symnmf", symnmf_module_imp, METH_VARARGS, PyDoc_STR("An implementation of symnmf.")},
    {"get_list", get_list, METH_VARARGS, PyDoc_STR("get list of ints")},
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
