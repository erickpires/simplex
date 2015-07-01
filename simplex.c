#include <stdio.h>
#include <stdlib.h>
#include "simplex.h"

#define abs(n) (n < 0 ? -n : n)

int get_PPL_from_file(FILE* input, PPL* ppl) {
	uint n_variables;
	uint n_restrictions;

	printf("# of variables: ");
	fscanf(input, "%d", &n_variables);

	printf("\n# of restrictions: ");
	fscanf(input, "%d", &n_restrictions);

	printf("\nRestriction type ( < | > | = ): ");
	char tmp_char = 0;
	do {
		fscanf(input, "%c", &tmp_char);
	} while(tmp_char != '<' && tmp_char != '>' && tmp_char != '=');

	ppl->type = (int) tmp_char;

	ppl->A.cols = n_variables;
	ppl->A.lines = n_restrictions;
	ppl->c.size = n_variables;
	ppl->b.size = n_restrictions;

	ppl->A.data = (double*) malloc(n_variables * n_restrictions * sizeof(double));
	ppl->c.data = (double*) malloc(n_variables * sizeof(double));
	ppl->b.data = (double*) malloc(n_restrictions * sizeof(double));

	printf("Enter the function coefficients\n");
	for(uint i = 0; i < n_variables; i++) {
		double tmp_double;
		fscanf(input, "%lf", &tmp_double);
		set_vector_value(&(ppl->c), i, tmp_double);
	}

	printf("Enter the \"A\" matrix values (per line)\n");
	for(uint i = 0; i < n_restrictions; i++) {
		for(uint j = 0; j < n_variables; j++) {
			double tmp_double;
			fscanf(input, "%lf", &tmp_double);
			set_matrix_value(&(ppl->A), i, j, tmp_double);
		}
	}

	printf("Enter the \"b\" vector values\n");
	for(uint i = 0; i < n_restrictions; i++) {
		double tmp_double;
		fscanf(input, "%lf", &tmp_double);
		set_vector_value(&(ppl->b), i, tmp_double);
	}

	return 1;
}

void print_PPL(PPL* ppl) {
	Vector* c = &(ppl->c);
	Vector* b = &(ppl->b);
	Matrix* A = &(ppl->A);

	printf("Min\n    ");
	print_vector_as_coefficients(c);

	printf("\nsubject to\n");
	for(uint i = 0; i < A->lines; i++) {
		Vector tmp_vector;
		get_vector_from_matrix_line(A, &tmp_vector, i);
		printf("    ");
		print_vector_as_coefficients(&tmp_vector);	

		printf(" %c %.2lf\n", ppl->type, get_vector_value(b, i));
	}
}

void print_vector_as_coefficients(Vector* vector) {
	if(vector->size == 0) return;

	printf("%.2lf x%d", get_vector_value(vector, 0), 0);
	for(uint i = 1; i < vector->size; i++) {
		double value = get_vector_value(vector, i);
		if(value < 0) 
			printf(" - ");
		else
			printf(" + ");

		printf("%.2lf x%d", abs(value), i);
	}
}

inline void set_matrix_value(Matrix* matrix, uint i, uint j, double value) {
	matrix->data[i * matrix->cols + j] = value;
}

inline void set_vector_value(Vector* vector, uint i, double value) {
	vector->data[i] = value;
}

inline double get_matrix_value(Matrix* matrix, uint i, uint j) {
	return matrix->data[i * matrix->cols + j];
}

inline double get_vector_value(Vector* vector, uint i) {
	return vector->data[i];
}

void get_vector_from_matrix_line(Matrix* matrix, Vector* vec, uint line) {
	vec->size = matrix->cols;
	vec->data = matrix->data + (line * matrix->cols);
}

void inner_product(Vector* a, Vector* b, Vector* out) {
	if(a->size != b->size || a->size != out->size) return;

	for(uint i = 0; i < a->size; i++) {
		//TODO: Use pointer arithmetic here
		out->data[i] = a->data[i] * b->data[i];
	}
}

void expand_PPL(PPL* ppl) {
	// If the problem already has '=' restrictions, there's nothing to be done
	if(ppl->type == equal) 
		return;

	//TODO: A (and probably c) should be extract in this local context
	//TODO: Calculate the diff from the value below
	int new_n_varibles;
	if (ppl->type == less_than)
		new_n_varibles = ppl->A.cols + ppl->A.lines;
	else // greater_than
		new_n_varibles = ppl->A.cols + 2 * ppl->A.lines;

	ppl->c.data = realloc(ppl->c.data, new_n_varibles * sizeof(double));
	for(uint i = ppl->c.size; i < new_n_varibles; i++) {
		set_vector_value(&(ppl->c), i, 0.0);
	}
	ppl->c.size = new_n_varibles;

	Matrix new_matrix;
	new_matrix.cols = new_n_varibles;
	new_matrix.lines = ppl->A.lines;
	new_matrix.data = (double*) malloc(ppl->A.lines * new_n_varibles * sizeof(double));

	for(uint i = 0; i < ppl->A.lines; i++) {
		for(uint j = 0; j < ppl->A.cols; j++) {
			set_matrix_value(&new_matrix, i, j, get_matrix_value(&(ppl->A), i, j));
		}
	}

	for(uint i = 0; i < ppl->A.lines; i++) {
		for(uint j = ppl->A.cols; j < new_n_varibles; j++) {
			uint diff = j - ppl->A.cols;
			double new_value;
			if(diff == i) {
				if(ppl->type == less_than)
					new_value = 1.0;
				else
					new_value = -1.0;
			}
			else if(diff == i + ppl->A.lines)
				new_value = 1.0;
			else
				new_value = 0.0;
			
			set_matrix_value(&new_matrix, i, j, new_value);
		}
	}

	free(ppl->A.data);
	ppl->A = new_matrix;
	ppl->type = equal;
}

int main(int argc, char** argv){

	PPL ppl;

	get_PPL_from_file(stdin, &ppl);

	printf("\n\n\n\n\n");
	print_PPL(&ppl);

	//TODO: Will probably need to save tha PPL type before expanding it
	PPL expanded_ppl;
	expand_PPL(&ppl);

	print_PPL(&ppl);

    return 0;
}
