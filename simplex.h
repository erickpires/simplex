#ifndef SIMPLEX_H
#define SIMPLEX_H 1

typedef unsigned int uint;

typedef enum {
	less_than = '<',
	greater_than = '>',
	equal = '='
} Restriction;

typedef enum {
	single_solution,
	mutiple_solution,
	unbounded,
	unfeasible,
	feasible
}LPP_type;

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
	Restriction* restrictions_type;
	int* identity_cols_indexes;
	Matrix A;
	Vector c;
	Vector b;
	uint slack_variables;
}PPL;

typedef struct {
	uint* variables_in_base;
	Matrix table;
	Vector costs;
	uint n_variables_in_base;
} Simplex_table;

typedef struct {
	uint* values;
	uint* tmp_base;
	uint variables_per_base;
	uint n_bases;
	uint capacity_in_bases;
} Remembered_bases;

void set_matrix_value(Matrix*, uint, uint, double);
void set_vector_value(Vector*, uint, double);

double get_matrix_value(Matrix*, uint, uint);
double get_vector_value(Vector*, uint);

void copy_matrix_with_line_offset(Matrix*, Matrix*, uint, uint);
void copy_vector(Vector*, Vector*);
double inner_product(Vector*, Vector*);
void sum_vector_with_multiplier(Vector*, Vector*, double);
void divide_vector_by_scalar(Vector*, double);

void get_vector_from_matrix_line(Matrix*, Vector*, uint);
void get_vector_from_matrix_col_with_offset(Matrix*, Vector*, uint, uint);

void print_vector_as_coefficients(Vector*);
void print_PPL(PPL*);

int get_PPL_from_file(FILE*, PPL*);

void expand_PPL(PPL*);

void get_c_b(Simplex_table*, Vector*);
void fill_simplex_table_z_line(Simplex_table*, Vector*);

int get_first_phase_table(Simplex_table*, PPL*);

LPP_type run_simplex(Simplex_table*, int, int);
LPP_type lpp_type_from_solved_table(Simplex_table*);
LPP_type get_second_phase_table(Simplex_table*, PPL*);

void print_other_solutions_from_base_solution(Simplex_table*, Remembered_bases*);

int is_variable_in_base(Simplex_table*, uint);

void put_in_base(Simplex_table*, uint, uint);

void print_simplex_table(Simplex_table*);

int is_same_base(uint*, uint*, uint);

void alloc_remembered_bases(Remembered_bases*, uint);
void realloc_remembered_bases(Remembered_bases*);

void sort_tmp_base(Remembered_bases*);
void remember_tmp_base(Remembered_bases*);
void copy_to_tmp_base(Remembered_bases*, uint*);
int is_base_remembered(Remembered_bases*, uint*, uint, uint);

#endif
