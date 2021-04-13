//==============================================================================
// Recursive definition of determinate using expansion by minors.
//
// Notes: 1) arguments:
//             a (double **) pointer to a pointer of an arbitrary square matrix
//             n (int) dimension of the square matrix
//
//        2) Determinant is a recursive function, calling itself repeatedly
//           each time with a sub-matrix of the original till a terminal
//           2X2 matrix is achieved and a simple determinat can be computed.
//           As the recursion works backwards, cumulative determinants are
//           found till untimately, the final determinate is returned to the
//           initial function caller.
//
//        3) m is a matrix (4X4 in example)  and m13 is a minor of it.
//           A minor of m is a 3X3 in which a row and column of values
//           had been excluded.   Another minor of the submartix is also
//           possible etc.
//             m  a b c d   m13 . . . .
//                e f g h       e f . h     row 1 column 3 is elminated
//                i j k l       i j . l     creating a 3 X 3 sub martix
//                m n o p       m n . p
//
//        4) the following function finds the determinant of a matrix
//           by recursively minor-ing a row and column, each time reducing
//           the sub-matrix by one row/column.  When a 2X2 matrix is
//           obtained, the determinat is a simple calculation and the
//           process of unstacking previous recursive calls begins.
//
//                m n
//                o p  determinant = m*p - n*o
//
//        5) this function uses dynamic memory allocation on each call to
//           build a m X m matrix  this requires **  and * pointer variables
//           First memory allocation is ** and gets space for a list of other
//           pointers filled in by the second call to malloc.
//
//        6) C++ implements two dimensional arrays as an array of arrays
//           thus two dynamic malloc's are needed and have corresponsing
//           free() calles.
//
//        7) the final determinant value is the sum of sub determinants
//
//==============================================================================
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h> 

double Determinant(double** matrixA, int size);
void inv_matrix(double **matrix, int size, double det, double**inv);

int main(void)
{
	int size = 0;
	int row = 0;
	int column = 0;
	FILE* fa;
	double** matrixa;
	double det = 0;
	double** matrixa_inv;

	fa = fopen("matrix.txt", "r");
	fscanf(fa, "%d", &size);
	printf("size = %d x %d \n", size, size);

	matrixa = (double**)malloc(sizeof(double*) * size);
	matrixa_inv = (double**)malloc(sizeof(double*) * size);

	for (row = 0; row < size; row++)
	{
		*(matrixa + row) = (double*)malloc(sizeof(double) * size);
		*(matrixa_inv + row) = (double*)malloc(sizeof(double) * size);
	}

	for (row = 0; row < size; row++)
	{
		for (column = 0; column < size; column++)
		{
			fscanf(fa, "%lf", *(matrixa + row) + column);
			printf("%lf ", matrixa[row][column]);
		}
		printf("\n");
	}
	printf("\n");


	det = Determinant(matrixa, size);
	printf("\ndet = %.5g\n", det);

	row = 0; // row, column �ʱ�ȭ 
	column = 0;

	//<matrix inverse���ϱ�>
	printf("\n");
	printf("��� matrix inverse");
	printf("\n");

	for (row = 0; row < size; row++) //matrix_inv������� ��
	{
		for (column = 0; column < size; column++) //matrix_inv���࿭�� ��
		{
			matrixa_inv[row][column] = 0;//matrix_inv�࿭ �ʱ�ȭ
		}
	}

	inv_matrix(matrixa, size, det, matrixa_inv);//inv_matrix�Լ��� �̿��Ͽ� matrix�� ������� ���

	//matrix�� ����� ���
	for (row = 0; row < size; row++) {

		for (column = 0; column < size; column++) {
			printf("%lf ", matrixa_inv[row][column]);
		}
		printf("\n");
	}


	free(matrixa);

	return 0;
}


double Determinant(double** a, int n)
{
	int i, j, j1, j2;                    // general loop and matrix subscripts
	double det = 0;                   // init determinant
	double** m = NULL;                // pointer to pointers to implement 2d
							  // square array
	int sign = 1;

	if (n < 1) {}                                                    // error condition, should never get here

	else if (n == 1) {                                              // should not get here
		det = a[0][0];
	}

	else if (n == 2) {                                            // basic 2X2 sub-matrix determinate
														   // definition. When n==2, this ends the
		det = a[0][0] * a[1][1] - a[1][0] * a[0][1];              // the recursion series
	}


	// recursion continues, solve next sub-matrix
	else {                             // solve the next minor by building a
							   // sub matrix
		det = 0;                      // initialize determinant of sub-matrix

							   // for each column in sub-matrix
		for (j1 = 0; j1 < n; j1++) {
			// get space for the pointer list
			m = (double**)malloc((n - 1) * sizeof(double*));

			for (i = 0; i < n - 1; i++)
				m[i] = (double*)malloc((n - 1) * sizeof(double));

			for (i = 1; i < n; i++) {
				j2 = 0;               // start at first sum-matrix column position
								 // loop to copy source matrix less one column
				for (j = 0; j < n; j++) {
					if (j == j1) continue; // don't copy the minor column element

					m[i - 1][j2] = a[i][j];  // copy source element into new sub-matrix
									   // i-1 because new sub-matrix is one row
									   // (and column) smaller with excluded minors
					j2++;                  // move to next sub-matrix column position
				}
			}


			if (j1 % 2 == 1)
			{
				sign = -1;
			}
			else
			{
				sign = 1;
			} /* else */

			det += sign * a[0][j1] * Determinant(m, n - 1);
			// sum x raised to y power
			// recursively get determinant of next
			// sub-matrix which is now one
			// row & column smaller

			for (i = 0; i < n - 1; i++) free(m[i]);// free the storage allocated to
										  // to this minor's set of pointers
			free(m);                       // free the storage for the original
									// pointer to pointer
		}
	}
	return(det);
}

void inv_matrix(double **matrix, int size, double det, double**inv)
{
	int i, j, k, L, co_i, co_j, cofac_size = size - 1;
	double **cofactor;
	double **matrix_a;
	double cofac_det;
	

	//cofactor �࿭�� ũ�� �޸� �Ҵ�
	cofactor = (double**)malloc(sizeof(double*)*cofac_size);
	//matrix_a �࿭�� ũ�� �޸� �Ҵ�
	matrix_a = (double**)malloc(sizeof(double*)*size);

	for (co_i = 0; co_i < cofac_size; co_i++) // ���Ϸ��� matrix ���� 1 ���� size�� cofactor�� ���ϱ� ���� ��Ľ�
		*(cofactor + co_i) = (double*)malloc(sizeof(double)*cofac_size);

	for (i = 0; i < size; i++) 
		*(matrix_a + i) = (double*)malloc(sizeof(double)*size);

	i = 0;
	j = 0;

	//inv ����� ���ϱ� ���� cofactor(���μ����) ���ϱ�
	for (k = 0; k < size; k++)
	{
		for (L = 0; L < size; L++)
		{
			for (co_i = 0, co_j = 0; co_i < cofac_size; )
			{
				if ((i == k) || (j == L)) // ���õ� ��� ���� ������ submatrix ã��
				{
					j++;
					if (j >= size)
					{
						j = 0;
						i++;
					}
				}

				else {
					cofactor[co_i][co_j] = matrix[i][j]; //ã���� ����
					co_j++;
					j++;
					if (j >= size)
					{
						j = 0;
						i++;
					}
					if (co_j >= cofac_size)
					{
						co_j = 0;
						co_i++;

					}
				}
			}

			i = 0; // submatrix�� ã���� �ʱ�ȭ , if) a11�� ���� ��Ľ��� ã���� �ٽ� a12�� ���� ��Ľ��� ã�� ���� �ʱ�ȭ
			j = 0;
			cofac_det = Determinant(cofactor, cofac_size); // submatrix(���õ� ��� ���� ������ ��� ���� �� ��Ľ�)�� det 
			inv[k][L] = cofac_det / det; //cofactor�� ����

		}
	}
	for (i = 0; i < size; i++) //cofactor�� matrix_a ��Ŀ� ����
	{
		for (j = 0; j < size; j++)
		{
			matrix_a[i][j] = inv[i][j];
		}
	}

	for (i = 0; i < size; i++) // ������ �ٽ� �ٲ���,inv[i][j]=(-1)^(i+j)*matrix_a[i][j] ����, ���� inverse matrix�� ����
	{
		for (j = 0; j < size; j++)
		{
			inv[j][i] = matrix_a[i][j]; // cofactor�� adjoint(�������)�� ���� ���� �ڹٲ�
			if (((i + j) % 2) == 1 && inv[j][i] != 0) //i+j�� Ȧ���� ��� -, ���� inverse�� ���Ұ��� 0 �ش���� �ʰ� ����
			{
				inv[j][i] *= (double)(-1);
			}
		}
	}
}
