#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "symnmf.h"

/* dimension of data points */
static int d = 1;
/* number of data points */
static int N = 0;
/* input argument */
static int k = 0;

/* constants for the SymNMF algorithm */
static const double beta = 0.5;
static const int max_iter = 1000;
static const double eps = 0.0001;

static struct vector *backup_vectors;

/* Computational functions */
double squared_dist(struct vector u, struct vector v);
double** matrix_multiply(double** A, double** B, int m, int n, int p);
double** transpose(double** matrix, int rows, int cols);

void print_vector(struct vector vector);

/*note in 2.9.2 in project_NMF.pdf we are told we don't need to validate args.*/
int main(int argc, char *argv[]) {

    char* file_path = argv[3];
    char* goal = argv[2];
    struct vector *data_points;
    double **A;
    double *D;
    double **W;

    if (argc != 4) {
        printf("Usage: %s <path_to_txt_file> <k>\n", argv[0]);
        return 1;
    }

    /* Read k from command line arguments */
    k = atoi(argv[1]);
    if (k <= 0) {
        printf("Invalid value for k. It should be a positive integer.\n");
        return 1;
    }

    /* Read data points */
    data_points = read_data_points(file_path);

    /* Compute matrices */
    A = compute_similarity_matrix(data_points, N, d);
    D = compute_degree_matrix(A);
    W = compute_laplacian_matrix(A, D);

    if (!strcmp(goal, "sym")) /*why we need to negate?*/
        print_matrix(A, N, N);

    if (!strcmp(goal, "ddg"))
        print_diagonal_matrix(D, N);

    if (!strcmp(goal, "norm"))
        print_matrix(W, N, N);

    /* Free allocated memory */
    free_matrix(A, N);
    free_matrix(W, N);
    free(D);
    free_vectors(data_points);
    free_vectors(backup_vectors);

    return 0;
}

struct vector* read_data_points(char* file_name){

    FILE *ifp = NULL;
    struct vector *head_vec, *curr_vec;
    struct entry *head_entry, *curr_entry;
    double n;
    char c;

    ifp = fopen(file_name, "r" );

    head_entry = malloc(sizeof(struct entry));
    if (head_entry == NULL)     /* Memory allocation failed */
        mem_error();
    curr_entry = head_entry;
    curr_entry->next = NULL;

    curr_vec = malloc(sizeof(struct vector));
    if (curr_vec == NULL)       /* Memory allocation failed */
        mem_error();
    curr_vec->next = NULL;
    head_vec = curr_vec;
    backup_vectors = head_vec;

    while (fscanf(ifp, "%lf%c", &n, &c) == 2) {

        /* We have read all the entries for the current vector */
        if (c == '\n') {

            curr_entry->value = n;
            curr_vec->entries = head_entry;
            curr_vec->next = calloc(1, sizeof(struct vector));
            if (curr_vec->next == NULL)     /* Memory allocation failed */
                mem_error();
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_entry = malloc(sizeof(struct entry));
            if (head_entry == NULL)         /* Memory allocation failed */
                mem_error();
            curr_entry = head_entry;
            curr_entry->next = NULL;

            /* Count the number of vectors N */
            N++;
            continue;
        }

        /* Read the next entry of the current vector */
        curr_entry->value = n;
        curr_entry->next = malloc(sizeof(struct entry));
        if (curr_entry == NULL)         /* Memory allocation failed */
            mem_error();
        curr_entry = curr_entry->next;
        curr_entry->next = NULL;

        /* Count the dimension d */
        if (N == 0)
            d++;
    }

    free_entries(head_entry);
    fclose(ifp);

    return head_vec;
}

double** compute_similarity_matrix(struct vector* vectors, int N_, int d_) {
    int i = 0, j = 0;
    struct vector *vec_i = vectors, *vec_j = vectors;
    double **A;
    N = N_;
    d = d_;
    A = (double **) malloc(N * sizeof(double *));
    if (A == NULL)      /* Memory allocation failed */
        mem_error();

    for (; i < N; vec_i = vec_i->next, i++) {
        A[i] = (double *) malloc(N * sizeof(double));
        for (j = 0; j < N; vec_j = vec_j->next, j++) {
            if (i != j){
                A[i][j] = exp(-squared_dist(*vec_i, *vec_j) / 2); 
            }
            else
                A[i][j] = 0.0;
        }
        vec_j = vectors;
    }
    return A;
}

double* compute_degree_matrix(double** A) {
    int i = 0, j = 0;
    double *D = (double *)malloc(N * sizeof(double));
    if (D == NULL)      /* Memory allocation failed */
        mem_error();

    for (; i < N; i++) {
        D[i] = 0.0;
        for (j = 0; j < N; j++)
            D[i] += A[i][j];
    }
    return D;
}

double** compute_laplacian_matrix(double** A, double* D) {
    int i = 0;

    double **intermediate;
    double **W;

    /* Compute D^{-1/2} as a diagonal matrix */
    double** D_inv_sqrt = (double**)malloc(N * sizeof(double*));
    if (D_inv_sqrt == NULL)      /* Memory allocation failed */
        mem_error();

    for (; i < N; i++) {
        D_inv_sqrt[i] = (double*)calloc(N, sizeof(double));
        if (D_inv_sqrt[i] == NULL)      /* Memory allocation failed */
            mem_error();
        D_inv_sqrt[i][i] = 1.0 / sqrt(D[i]);
    }

    /* Compute intermediate matrix: D^{-1/2} * A */
    intermediate = matrix_multiply(D_inv_sqrt, A, N, N, N);

    /* Compute W = (D^{-1/2} * A) * D^{-1/2} */
    W = matrix_multiply(intermediate, D_inv_sqrt, N, N, N);

    /* Free allocated memory */
    free_matrix(D_inv_sqrt, N);
    free_matrix(intermediate, N);

    return W;
}

double** optimize_H(double** W, double** H, int n, int k) {
    int i = 0, j = 0;
    double** H_transpose = transpose(H, n, k);
    double** temp1 = NULL;
    double** temp2 = NULL;
    double** temp3 = NULL;
    double** H_prev = (double**)malloc(n * sizeof(double*));
    int iteration;
    double diff;
    if (H_prev == NULL)
        mem_error();

    for (; i < n; i++) {
        H_prev[i] = (double*)malloc(k * sizeof(double));
        if (H_prev[i] == NULL)
            mem_error();
    }

    iteration = 0;
    while (iteration < max_iter) {
        /* Copy current H to H_prev */
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++)
                H_prev[i][j] = H[i][j];
        }

        /* Update H using the provided rule */
        temp1 = matrix_multiply(W, H, n, n, k);
        temp2 = matrix_multiply(H, H_transpose, n, k, n);
        temp3 = matrix_multiply(temp2, H, n, n, k);

        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++)
                H[i][j] *= (1 - beta + (beta * temp1[i][j] / temp3[i][j]));
        }

        /* Check for convergence */
        diff = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++)
                diff += pow(H[i][j] - H_prev[i][j], 2);
        }

        /* Free temp matrices */
        free_matrix(temp1, n);
        free_matrix(temp2, n);
        free_matrix(temp3, n);

        if (diff < eps)
            break;

        H_transpose = transpose(H, n, k);
        iteration++;

    }

    /* Free allocated memory for H_prev */
    free_matrix(H_prev, n);

    return H;
}

void mem_error(void) {
    printf("An Error Has Occurred\n");
    exit(1);
}

void free_vectors(struct vector *head) {
    if (head != NULL){
        free_entries(head->entries);
        free_vectors(head->next);
        free(head);
    }
    head = NULL;
    backup_vectors = head;
}

void free_entries(struct entry *head) {
    if (head != NULL){
        free_entries(head->next);
        free(head);
    }
    head = NULL;
}

void free_matrix(double **matrix, int rows) {
    int i = 0;
    for (; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

void print_matrix(double **matrix, int rows, int cols) {
    int i = 0, j = 0;
    for (; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j != cols - 1)
                printf(",");
        }
        printf("\n");
    }
}

void print_diagonal_matrix(double *diagonal, int size) {
    int i = 0, j = 0;
    for (; i < size; i++) {
        for (j = 0; j < size; j++) {
            if (i == j)
                printf("%.4f", diagonal[i]);
            else
                printf("0.0000");
            if (j != size - 1)
                printf(",");
        }
        printf("\n");
    }
}

double squared_dist(struct vector u, struct vector v) {
    int i = 0;

    struct entry *u_entry = u.entries;
    struct entry *v_entry = v.entries;
    double sum = 0;

    for(; i < d; i++) {
        sum += pow(u_entry->value - v_entry->value, 2);
        u_entry = u_entry->next;
        v_entry = v_entry->next;
    }

    return sum;
}

double** matrix_multiply(double** A, double** B, int m, int n, int p) {
    int i = 0, j = 0, l = 0;
    double **C = (double **) malloc(m * sizeof(double*));
    if (C == NULL)         /* Memory allocation failed */
        mem_error();

    for (; i < m; i++) {
        C[i] = (double *) malloc(m * sizeof(double));
        if (C[i] == NULL)         /* Memory allocation failed */
            mem_error();

        for (j = 0; j < p; j++) {
            C[i][j] = 0.0;
            for (l = 0; l < n; l++)
                C[i][j] += A[i][l] * B[l][j];
        }
    }
    return C;
}

double** transpose(double** matrix, int rows, int cols) {
    int i = 0, j = 0;
    double **transposed = (double**)malloc(cols * sizeof(double*));
    if (transposed == NULL)
        mem_error();

    for (; i < cols; i++) {
        transposed[i] = (double*)malloc(rows * sizeof(double));
        if (transposed[i] == NULL)
            mem_error();

        for (j = 0; j < rows; j++)
            transposed[i][j] = matrix[j][i];
    }

    return transposed;
}

void print_vectors(struct vector *vectors) {
    struct vector* curr_vec = vectors;

    while (curr_vec != NULL) {
        struct entry* entry = curr_vec->entries;
        while (entry != NULL) {
            printf("%.4f ", entry->value);
            entry = entry->next;
        }
        printf("\n");
        fflush(stdout); 
        curr_vec = curr_vec->next;
    }
}

void print_vector(struct vector vector) {
    struct entry* entry = vector.entries;
    while (entry != NULL) {
        printf("%.4f ", entry->value);
        entry = entry->next;
    }
    printf("\n");
    fflush(stdout); 
}
