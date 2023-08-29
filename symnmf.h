#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

static int k = 0;
static int N = 0;
static int d = 1;
static struct vector *backup_vectors;

const double beta = 0.5;
const int max_iter = 1000;
const double eps = 0.0001;

/* Structs definitions */
struct entry {
    double value;
    struct entry *next;
};

struct vector {
    struct vector *next;
    struct entry *entries;
};

/* Goals */
double** compute_similarity_matrix(struct vector* vectors);
double* compute_degree_matrix(double** A);
double** compute_laplacian_matrix(double** A, double* D);
double** initialize_H(double** W);
double** optimize_H(double** W);

/* Argument reading and processing functions */
struct vector* read_data_points(char* file_name);
void mem_error();
void free_entries(struct entry *head);
void free_matrix(double** matrix, int rows);
void print_matrix(double** matrix, int rows, int cols);
void print_diagonal_matrix(double* diagonal, int size);

/* Computational functions */
double squared_dist(struct vector u, struct vector v);
double** matrix_multiply(double** A, double** B, int m, int n, int p);
double** transpose(double** matrix, int rows, int cols);


/* Functions implemented in symnmfmodule.c */
static PyObject* sym_module_imp(PyObject *self, PyObject *args);

static PyObject* ddg_module_imp(PyObject *self, PyObject *args);

static PyObject* norm_module_imp(PyObject *self, PyObject *args);

static PyObject* symnmf_module_imp(PyObject *self, PyObject *args);

PyMODINIT_FUNC PyInit_symnmf(void);
