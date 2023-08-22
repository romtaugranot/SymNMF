static PyObject* sym_module_imp(PyObject *self, PyObject *args)
{

    if(!PyArg_ParseTuple(args, "", ...)) {
        return NULL;
    }

    return Py_BuildValue("O", NULL);
}

static PyObject* ddg_module_imp(PyObject *self, PyObject *args)
{

    if(!PyArg_ParseTuple(args, "", ...)) {
        return NULL;
    }

    return Py_BuildValue("O", NULL);
}

static PyObject* norm_module_imp(PyObject *self, PyObject *args)
{

    if(!PyArg_ParseTuple(args, "", ...)) {
        return NULL;
    }

    return Py_BuildValue("O", NULL);
}

static PyObject* symnmf_module_imp(PyObject *self, PyObject *args)
{

    if(!PyArg_ParseTuple(args, "", ...)) {
        return NULL;
    }

    return Py_BuildValue("O", NULL);
}

static PyMethodDef symnmfMethods[] = {
    {"sym",(PyCFunction) sym_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the similarity matrix.")},
    {"ddg",(PyCFunction) ddg_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the Diagonal Degree Matrix.")},
    {"norm",(PyCFunction) norm_module_imp, METH_VARARGS, PyDoc_STR("Calculate and output the normalized similarity matrix.")},
    {"symnmf",(PyCFunction) symnmf_module_imp, METH_VARARGS, PyDoc_STR("An implementation of symnmf.")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    symnmfMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}