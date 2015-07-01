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
	ppl->slack_variables = 0;
	ppl->identity_cols_indexes = (int*) malloc(n_restrictions * sizeof(uint));

	for(uint i = 0; i < n_restrictions; i++) {
		ppl->identity_cols_indexes[i] = -1;
	}

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

void print_simplex_table(Simplex_table* simplex_table) {
	printf("    | z | ");
	for(uint i = 0; i < simplex_table->table.cols - 1; i++)
		printf("  x%d       ", i);
	printf("| RHS\n");

	for(uint i = 0; i < simplex_table->table.cols; i++)
		printf("============", i);
	printf("\n");

	printf(" z  | 1 | ");
	{
		uint i;
		for(i = 0; i < simplex_table->table.cols - 1; i++){
			printf(" %-9.2lf ", get_matrix_value(&(simplex_table->table), 0, i));
		}
		printf("| %-9.2lf\n", get_matrix_value(&(simplex_table->table), 0, i));
	}

	for(uint i = 0; i < simplex_table->table.cols; i++)
		printf("============", i);
	printf("\n");

	for(uint i = 0; i < simplex_table->n_variables_in_base; i++) {
		printf(" x%d | 0 | ", simplex_table->variables_in_base[i]);
		uint j;
		for(j = 0; j < simplex_table->table.cols - 1; j++) {
			printf(" %-9.2lf ", get_matrix_value(&(simplex_table->table), i + 1, j));
		}
		printf("| %-9.2lf\n", get_matrix_value(&(simplex_table->table), i + 1, j));
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

	Matrix* A = &(ppl->A);
	Vector* c = &(ppl->c);
	
	int new_n_varibles = A->cols + A->lines;
	ppl->slack_variables = A->lines;

	c->data = realloc(c->data, new_n_varibles * sizeof(double));
	for(uint i = c->size; i < new_n_varibles; i++) {
		set_vector_value(c, i, 0.0);
	}

	c->size = new_n_varibles;

	Matrix new_matrix;
	new_matrix.cols = new_n_varibles;
	new_matrix.lines = A->lines;
	new_matrix.data = (double*) malloc(A->lines * new_n_varibles * sizeof(double));

	for(uint i = 0; i < A->lines; i++) {
		for(uint j = 0; j < A->cols; j++) {
			set_matrix_value(&new_matrix, i, j, get_matrix_value(A, i, j));
		}
	}

	for(uint i = 0; i < A->lines; i++) {
		for(uint j = A->cols; j < new_n_varibles; j++) {
			uint diff = j - A->cols;
			double new_value;
			if(diff == i) {
				if(ppl->type == less_than) {
					new_value = 1.0;
					ppl->identity_cols_indexes[i] = j;
				}
				else
					new_value = -1.0;
			}
			else
				new_value = 0.0;
			
			set_matrix_value(&new_matrix, i, j, new_value);
		}
	}

	free(A->data);
	ppl->A = new_matrix;
	ppl->type = equal;
}

void fill_simplex_table(Simplex_table* simplex_table, PPL* ppl) {
	simplex_table->n_variables_in_base = ppl->b.size;
	simplex_table->variables_in_base = (uint*) malloc(simplex_table->n_variables_in_base * sizeof(uint));

	simplex_table->table.cols = ppl->A.cols + 1;
	simplex_table->table.lines = ppl->A.lines + 1;

	size_t table_size_in_bites = simplex_table->table.cols * simplex_table->table.lines * sizeof(double);
	simplex_table->table.data = (double*) malloc(table_size_in_bites);

	for(uint i = 0; i < simplex_table->n_variables_in_base; i++) {
		if(ppl->identity_cols_indexes[i] == -1) {
			// TODO: An artificial variable is needed
		}
		else {
			simplex_table->variables_in_base[i] = ppl->identity_cols_indexes[i];
		}
	}

	for(uint i = 0; i < ppl->c.size; i++) {
		//TODO: calculate c_B^t * b^(-1) * a_j - c_j
		set_matrix_value(&(simplex_table->table), 0, i, -(get_vector_value(&(ppl->c), i)));
	}

	//TODO(Maybe): Refactor to copy_matrix_to_matrix
	for(uint i = 0; i < ppl->A.lines; i++) {
		for(uint j = 0; j < ppl->A.cols; j++) {
			set_matrix_value(&(simplex_table->table), i + 1, j, get_matrix_value(&(ppl->A), i, j));
		}
	}

	uint last_table_col_index = simplex_table->table.cols - 1;
	//TODO(Maybe): Refactor to copy_vector_to_matrix_col
	for(uint i = 0; i < ppl->b.size; i++) {
		set_matrix_value(&(simplex_table->table), i + 1, last_table_col_index, get_vector_value(&(ppl->b), i));
	}

	//README: Not sure if this value will be always zero
	set_matrix_value(&(simplex_table->table), 0, last_table_col_index, 0.0);
}

int main(int argc, char** argv){

	PPL ppl = {};
	Simplex_table simplex = {};

	get_PPL_from_file(stdin, &ppl);

	printf("\n\n\n");
	print_PPL(&ppl);

	//TODO: Will probably need to save tha PPL type before expanding it
	expand_PPL(&ppl);
	print_PPL(&ppl);

	fill_simplex_table(&simplex, &ppl);
	printf("\n\n\n");
	print_simplex_table(&simplex);

    return 0;
}
