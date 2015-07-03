#ifndef SIMPLEX_H
#define SIMPLEX_H 1

typedef unsigned int uint;

typedef enum {
	less_than = '<',
	greater_than = '>',
	equal = '='
} Restriction;

typedef struct {
	double* data;
	uint cols;
	uint lines;
} Matrix;

typedef struct {
	double* data;
	uint size;
} Vector;

typedef struct {
	Matrix A;
	Vector c;
	Vector b;
	Restriction* restrictions_type;
	uint slack_variables;
	int* identity_cols_indexes;
}PPL;

typedef struct {
	Matrix table;
	uint n_variables_in_base;
	uint* variables_in_base;
} Simplex_table;

void set_matrix_value(Matrix*, uint, uint, double);
void set_vector_value(Vector*, uint, double);

double get_matrix_value(Matrix*, uint, uint);
double get_vector_value(Vector*, uint);

double inner_product(Vector*, Vector*);

void get_vector_from_matrix_line(Matrix*, Vector*, uint);
void get_vector_from_matrix_col_with_offset(Matrix*, Vector*, uint, uint);

void print_vector_as_coefficients(Vector*);
void print_PPL(PPL*);

int get_PPL_from_file(FILE*, PPL*);

void expand_PPL(PPL*);

//TODO: This function can return an int to indicate whether the 2-phase method is needed
int get_first_phase_table(Simplex_table*, PPL*);
void print_simplex_table(Simplex_table*);
#endif