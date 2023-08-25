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

/* Functions declarations */

int main(int argc, char *argv[]);

/* Goals */
double** compute_similarity_matrix(struct vector* vectors);
double* compute_degree_matrix(double** A);
double** compute_laplacian_matrix(double** A, double* D);
double** initialize_H(double** W);
double** optimize_H(double** W);

/* Argument reading and processing functions */
struct vector* read_data_points();
void mem_error();
void free_entries(struct entry *head);
void free_matrix(double** matrix, int rows);
void print_matrix(double** matrix, int rows, int cols);
void print_diagonal_matrix(double* diagonal, int size);

/* Computational functions */
double squared_dist(struct vector u, struct vector v);
double** matrix_multiply(double** A, double** B, int m, int n, int p);
double** transpose(double** matrix, int rows, int cols);

/* Algorithm */
int main(int argc, char *argv[]) {

    if (argc != 3) {
        printf("Usage: %s <path_to_txt_file> <k>\n", argv[0]);
        return 1;
    }

    /* Read k from command line arguments */
    k = atoi(argv[2]);
    if (k <= 0) {
        printf("Invalid value for k. It should be a positive integer.\n");
        return 1;
    }

    /* Redirect stdin to read from the file */
    freopen(argv[1], "r", stdin);
    if (!stdin) {
        printf("Failed to open file: %s\n", argv[1]);
        return 1;
    }

    /* Read data points */
    struct vector *data_points = read_data_points();

    /* Compute matrices */
    double **A = compute_similarity_matrix(data_points);
    double *D = compute_degree_matrix(A);
    double **W = compute_laplacian_matrix(A, D);

    /* Optimize H */
    double **H = optimize_H(W);

    /* Print matrices */
    printf("A:\n");
    print_matrix(A, N, N);
    printf("D:\n");
    print_diagonal_matrix(D, N);
    printf("W:\n");
    print_matrix(W, N, N);
    printf("H:\n");
    print_matrix(H, N, k);

    /* Free allocated memory */
    free_matrix(A, N);
    free_matrix(W, N);
    free(D);
    free_matrix(H, N);

    return 0;
}

double** compute_similarity_matrix(struct vector* vectors) {
    int i = 0, j = 0;
    struct vector *vec_i = vectors, *vec_j = vectors;
    double **A = (double **) malloc(N * sizeof(double *));
    if (A == NULL)      /* Memory allocation failed */
        mem_error();

    for (; i < N; vec_i = vec_i->next, i++) {
        A[i] = (double *) malloc(N * sizeof(double));
        for (j = 0; j < N; vec_j = vec_j->next, j++) {
            if (i != j)
                A[i][j] = exp(-squared_dist(*vec_i, *vec_j) / 2); /* MIGHT BE PROBLEMATIC */
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
    double **intermediate = matrix_multiply(D_inv_sqrt, A, N, N, N);

    /* Compute W = (D^{-1/2} * A) * D^{-1/2} */
    double **W = matrix_multiply(intermediate, D_inv_sqrt, N, N, N);

    /* Free allocated memory */
    free_matrix(D_inv_sqrt, N);
    free_matrix(intermediate, N);

    return W;
}

double** initialize_H(double** W) {
    int i = 0, j = 0;
    double m = 0.0;

    /* Calculate the average of all entries of W */
    for (; i < N; i++) {
        for (j = 0; j < N; j++)
            m += W[i][j];
    }
    m /= pow(N, 2);

    /* Allocate memory for H */
    double **H = (double**)malloc(N * sizeof(double*));
    if (H == NULL)
        mem_error();

    for (i = 0; i < N; i++) {
        H[i] = (double*)malloc(k * sizeof(double));
        if (H[i] == NULL)
            mem_error();

        /* Randomly initialize H with values from the interval [0, 2 * sqrt(m/k)] */
        for (j = 0; j < k; j++)
            H[i][j] = ((double)rand() / RAND_MAX) * 2 * sqrt(m / k);
    }
    return H;
}

double** optimize_H(double** W) {
    int i = 0, j = 0;
    double** H = initialize_H(W);
    double** H_transpose = transpose(H, N, k);
    double** temp1 = NULL;
    double** temp2 = NULL;
    double** temp3 = NULL;
    double** H_prev = (double**)malloc(N * sizeof(double*));
    if (H_prev == NULL)
        mem_error();

    for (; i < N; i++) {
        H_prev[i] = (double*)malloc(k * sizeof(double));
        if (H_prev[i] == NULL)
            mem_error();
    }

    int iteration = 0;
    while (iteration < max_iter) {
        /* Copy current H to H_prev */
        for (i = 0; i < N; i++) {
            for (j = 0; j < k; j++)
                H_prev[i][j] = H[i][j];
        }

        /* Update H using the provided rule */
        temp1 = matrix_multiply(W, H, N, N, k);
        temp2 = matrix_multiply(H, H_transpose, N, k, N);
        temp3 = matrix_multiply(temp2, H, N, N, k);

        for (i = 0; i < N; i++) {
            for (j = 0; j < k; j++)
                H[i][j] *= (1 - beta + (beta * temp1[i][j] / temp3[i][j]));
        }

        /* Check for convergence */
        double diff = 0.0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < k; j++)
                diff += pow(H[i][j] - H_prev[i][j], 2);
        }

        /* Free temp matrices */
        free_matrix(temp1, N);
        free_matrix(temp2, N);
        free_matrix(temp3, N);

        if (diff < eps)
            break;

        H_transpose = transpose(H, N, k);
        iteration++;
    }

    /* Free allocated memory for H_prev */
    free_matrix(H_prev, N);

    return H;
}

struct vector* read_data_points(){

    struct vector *head_vec, *curr_vec;
    struct entry *head_entry, *curr_entry;
    double n;
    char c;

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

    while (scanf("%lf%c", &n, &c) == 2) {

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

    return head_vec;
}

void mem_error(){
    printf("An Error Has Occurred\n");

    exit(1);
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
