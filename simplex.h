#ifndef SIMPLEX_H
#define SIMPLEX_H 1

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

//TODO: A PPL probably needs to now how many loosen (and virtual) variables it has
typedef struct {
	Matrix A;
	Vector c;
	Vector b;
	Restriction type;
}PPL;

void set_matrix_value(Matrix*, uint, uint, double);
void set_vector_value(Vector*, uint, double);

double get_matrix_value(Matrix*, uint, uint);
double get_vector_value(Vector*, uint);

void inner_product(Vector*, Vector*, Vector*);

void get_vector_from_matrix_line(Matrix*, Vector*, uint);

void print_vector_as_coefficients(Vector*);
void print_PPL(PPL*);

int get_PPL_from_file(FILE*, PPL*);

void expand_PPL(PPL*);
#endif