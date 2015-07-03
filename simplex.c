#include <stdio.h>
#include <stdlib.h>
#include "simplex.h"

#define abs(n) (n < 0 ? -n : n)

//TODO: create macros to get PPL variables (and restrictions) number

int get_PPL_from_file(FILE* input, PPL* ppl) {
	uint n_variables;
	uint n_restrictions;

	printf("# of variables: ");
	fscanf(input, "%d", &n_variables);

	printf("\n# of restrictions: ");
	fscanf(input, "%d", &n_restrictions);

	ppl->restrictions_type = (Restriction*) malloc(n_restrictions * sizeof(Restriction));
	for(uint i = 0; i < n_restrictions; i++) {
		printf("\nRestriction type ( < | > | = ): ");
		char tmp_char = 0;
		do {
			fscanf(input, "%c", &tmp_char);
		} while(tmp_char != '<' && tmp_char != '>' && tmp_char != '=');

		ppl->restrictions_type[i] = (Restriction) tmp_char;
	}
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

		printf(" %c %.2lf\n", ppl->restrictions_type[i], get_vector_value(b, i));
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
		printf("============");
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
		printf("============");
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
void get_vector_from_matrix_col_with_offset(Matrix* matrix, Vector* vector, uint offset, uint col)

double inner_product(Vector* a, Vector* b) {
	double result = 0.0;
	for(uint i = 0; i < a->size; i++) {
		//TODO: Use pointer arithmetic here
		result += a->data[i] * b->data[i];
	}

	return result;
}

void expand_PPL(PPL* ppl) {
	// TODO: Rewrite this function

	// using namespace PPL
	Matrix* A = &(ppl->A);
	Vector* c = &(ppl->c);

	// For each restriction inequality restriction, we will have a slack variable
	for(uint i = 0; i < ppl->b.size; i++) {
		if(ppl->restrictions_type[i] != equal)
			ppl->slack_variables++;
	}

	int new_n_varibles = A->cols + ppl->slack_variables;

	Matrix new_matrix;
	new_matrix.cols = new_n_varibles;
	new_matrix.lines = A->lines;
	new_matrix.data = (double*) malloc(A->lines * new_n_varibles * sizeof(double));

	// Copying the old matrix to the new one
	for(uint i = 0; i < A->lines; i++) {
		for(uint j = 0; j < A->cols; j++) {
			set_matrix_value(&new_matrix, i, j, get_matrix_value(A, i, j));
		}
	}

	uint n_slack_variables_added = 0;
	for(uint i = 0; i < A->lines; i++) {
		for(uint j = 0; j < ppl->slack_variables; j++) {
			uint col = A->cols + j;
			double new_value = 0.0;
			if(ppl->restrictions_type[i] != equal && n_slack_variables_added == j) {
				switch(ppl->restrictions_type[i]) {
					case less_than:
						new_value = 1.0;
						ppl->identity_cols_indexes[i] = col;
						break;
					case greater_than:
						new_value = -1.0;
						break;
					default:
						break;
				}
				n_slack_variables_added++;
				ppl->restrictions_type[i] = equal;
			}
			
			set_matrix_value(&new_matrix, i, col, new_value);
		}
	}

	free(A->data);
	ppl->A = new_matrix;

	c->data = realloc(c->data, new_n_varibles * sizeof(double));
	for(uint i = c->size; i < new_n_varibles; i++) {
		set_vector_value(c, i, 0.0);
	}

	c->size = new_n_varibles;
}

// This function returns a zero value if the 2-phase method is needed, otherwise it returns a non-zero value.
// In the former case the 'simplex_table' is prepared to run the first phase, in the latter it is prepared to
// run the second phase directly.
int get_first_phase_table(Simplex_table* simplex_table, PPL* ppl) {
	uint n_artificial_variables = 0;

	simplex_table->n_variables_in_base = ppl->b.size;
	simplex_table->variables_in_base = (uint*) malloc(simplex_table->n_variables_in_base * sizeof(uint));

	// The identity_cols_indexes were filled by calling expand_PPL
	for(uint i = 0; i < simplex_table->n_variables_in_base; i++) {
		if(ppl->identity_cols_indexes[i] == -1) {
			// The artificial variables will occupy the last columns of the simplex_table and will be put in the base
			simplex_table->variables_in_base[i] = ppl->c.size + n_artificial_variables;
			n_artificial_variables++;
		}
		else {
			// A slack variable already has an indentity column, so this variable is put in the base
			simplex_table->variables_in_base[i] = ppl->identity_cols_indexes[i];
		}
	}

	simplex_table->table.cols = ppl->A.cols + n_artificial_variables + 1;
	simplex_table->table.lines = ppl->A.lines + 1;

	size_t table_size = simplex_table->table.cols * simplex_table->table.lines;
	simplex_table->table.data = (double*) calloc(table_size, sizeof(double));

	//TODO(Maybe): Refactor to copy_matrix_to_matrix_with_offset
	for(uint i = 0; i < ppl->A.lines; i++) {
		for(uint j = 0; j < ppl->A.cols; j++) {
			set_matrix_value(&(simplex_table->table), i + 1, j, get_matrix_value(&(ppl->A), i, j));
		}
	}

	// Sets the artificial variables columns to the correspondent indentity column
	// Since the table is filled with zeroes by calloc only the ones are filled
	for(uint i = 0; i < simplex_table->n_variables_in_base; i++) {
		uint variable_index = simplex_table->variables_in_base[i];
		if(variable_index >= ppl->c.size) // if this varible is an artificial variable
			set_matrix_value(&(simplex_table->table), i + 1, variable_index, 1.0);
	}

	uint last_table_col_index = simplex_table->table.cols - 1;
	//TODO(Maybe): Refactor to copy_vector_to_matrix_col
	for(uint i = 0; i < ppl->b.size; i++) {
		set_matrix_value(&(simplex_table->table), i + 1, last_table_col_index, get_vector_value(&(ppl->b), i));
	}

	Vector c_b;
	c_b.size = simplex_table->n_variables_in_base;
	c_b.data = (double*) malloc(c_b.size * sizeof(double));

	Vector tmp_vector;
	tmp_vector.size = ppl->b.size;
	tmp_vector.data = (double*) malloc(tmp_vector.size * sizeof(double));

	for(uint i = 0; i < simplex_table->n_variables_in_base; i++) {
		uint variable_index = simplex_table->variables_in_base[i];
		if(n_artificial_variables) {
			if(variable_index >= ppl->c.size) // It's an artificial variable
				c_b.data[i] = 1.0;
			else
				c_b.data[i] = 0.0;
		}
		else {
			c_b.data[i] = get_vector_value(&(ppl->c), variable_index);
		}
	}

	for(uint i = 0; i < simplex_table->table.cols - 1; i++) {
		//TODO: calculate c_B^t * b^(-1) * a_j - c_j
		double cost;
		if(n_artificial_variables) {
			if(i >= ppl->c.size)
				cost = 1.0;
			else
				cost = 0.0;
		}
		else
			cost = get_vector_value(&(ppl->c), i);

		//README: This code can probably be reused
		// B^-1 is always the identity matrix, therefore the reduced cost is calculated with
		// (c_b)' * a_j - c_j
		get_vector_from_matrix_col_with_offset(&(simplex_table->table), &tmp_vector, 1, i);
		double minus_reduced_cost = inner_product(&c_b, &tmp_vector) - cost;
		set_matrix_value(&(simplex_table->table), 0, i, minus_reduced_cost);
	}

	//README: Not sure if this value will be always zero
	set_matrix_value(&(simplex_table->table), 0, last_table_col_index, 0.0);

	// TODO(Maybe): try to reuse c_b so it doesn't need to be freed and reallocated
	free(c_b.data);

	//README: By returning n_artificial_variables we achieve the function requirement of returning zero if the
	//        2-phase method is not need.
	return n_artificial_variables;
}

int main(int argc, char** argv){

	PPL ppl = {};
	Simplex_table simplex = {};

	get_PPL_from_file(stdin, &ppl);

	printf("\n\n\n");
	print_PPL(&ppl);

	expand_PPL(&ppl);
	print_PPL(&ppl);

	get_first_phase_table(&simplex, &ppl);
	printf("\n\n\n");
	print_simplex_table(&simplex);

    return 0;
}
