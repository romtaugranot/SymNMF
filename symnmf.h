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
double** optimize_H(double** W, double** H, int n, int k);

/* Argument reading and processing functions */
struct vector* read_data_points(char* file_name);

void print_matrix(double** matrix, int rows, int cols);
void print_diagonal_matrix(double* diagonal, int size);
void print_vectors(struct vector *vectors);

int getN();
int getK();
