// fem.cpp : Defines the entry point for the console application.
//
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>

typedef struct
{
	float x;
	float y;
} ASSEMBLY;

typedef struct
{
	int first_top;
	int second_top;
	int third_top;
} TRIANGLE;

int num_assemblies; //number of the nodes
ASSEMBLY *assemblies; //array of the nodes coordinates
int num_triangles; //number of the triangles
TRIANGLE *triangles; //array of the triangles tops
int num_fixed_assemblies; //number of the fixed nodes
int *fixed_assemblies; //array of the fixed nodes
double *MElastic; //stiffness matrix
double *MMomentum; //inertia matrix
double *K1, *K2; //combinations of the stiffness and the inertia matrices
double *Right; //stress tensor
double *q, *diffq, *a, *b, *qp;
double Emodule = 1.0; //Young's modulus
double Poisson = 0.3; //Poisson's ratio
double Rho = 1.0; //density

#define Number(x, y) ((x)*2 * num_assemblies + (y))

//construction of the stiffness and the inertia matrices
void BuildMatrix(void)
{
	double x[3], y[3]; //nodes coordinates of the current triangle
	double A[3], B[3], C[3]; //coefficients of the three functions Ах+Ву+С
	double V; //triangle area
	int Point[3]; //tops of the current triangle
	int i, j, k, n, nemo1, nemo2;
	double DKoefLame, EKoefLame, XKoefLame;

	DKoefLame = Emodule * (1.0 - Poisson) / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
	EKoefLame = Emodule * Poisson / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
	XKoefLame = Emodule * 0.5 / (1.0 + Poisson);

	for (i = 0; i < 2 * num_assemblies * 2 * num_assemblies; i++)
	{
		MElastic[i] = 0.0;
		MMomentum[i] = 0.0;
	}

	for (k = 0; k < num_triangles; k++)
	{

		Point[0] = triangles[k].first_top;
		Point[1] = triangles[k].second_top;
		Point[2] = triangles[k].third_top;

		for (i = 0; i < 3; i++)
		{
			x[i] = assemblies[Point[i]].x;
			y[i] = assemblies[Point[i]].y;
		}

		V = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);

		if (fabs(V) < 1.0e-03)
		{
			printf("Error fabs(V) < 1.0e-03!\n");
			exit(1);
		}

		A[0] = (y[1] - y[2]) / V;
		A[1] = (y[2] - y[0]) / V;
		A[2] = (y[0] - y[1]) / V;

		B[0] = (x[2] - x[1]) / V;
		B[1] = (x[0] - x[2]) / V;
		B[2] = (x[1] - x[0]) / V;

		C[0] = (x[1] * y[2] - x[2] * y[1]) / V;
		C[1] = (x[2] * y[0] - x[0] * y[2]) / V;
		C[2] = (x[0] * y[1] - x[1] * y[0]) / V;

		V = fabs(0.5 * V);

		//diagonal blocks
		for (i = 0; i < 3; i++)
		{
			nemo1 = 2 * Point[i];
			MElastic[Number(nemo1 + 1, nemo1 + 1)] += V * (DKoefLame * A[i] * A[i] + XKoefLame * B[i] * B[i]);
			MElastic[Number(nemo1 + 0, nemo1 + 0)] += V * (DKoefLame * B[i] * B[i] + XKoefLame * A[i] * A[i]);
			MElastic[Number(nemo1 + 1, nemo1 + 0)] += V * (EKoefLame * A[i] * B[i] + XKoefLame * A[i] * B[i]);
			MElastic[Number(nemo1 + 0, nemo1 + 1)] += V * (EKoefLame * A[i] * B[i] + XKoefLame * A[i] * B[i]);
			MMomentum[Number(nemo1 + 0, nemo1 + 0)] += Rho * V / 6.0;
			MMomentum[Number(nemo1 + 1, nemo1 + 1)] += Rho * V / 6.0;
		}

		//subdiagonal blocks
		for (i = 0; i < 2; i++)
		{
			nemo1 = 2 * Point[i];

			for (j = i + 1; j < 3; j++)
			{
				nemo2 = 2 * Point[j];
				MElastic[Number(nemo1 + 1, nemo2 + 1)] += V * (DKoefLame * A[i] * A[j] + XKoefLame * B[i] * B[j]);
				MElastic[Number(nemo1 + 0, nemo2 + 0)] += V * (DKoefLame * B[i] * B[j] + XKoefLame * A[i] * A[j]);
				MElastic[Number(nemo1 + 1, nemo2 + 0)] += V * (EKoefLame * A[i] * B[j] + XKoefLame * A[j] * B[i]);
				MElastic[Number(nemo1 + 0, nemo2 + 1)] += V * (EKoefLame * A[j] * B[i] + XKoefLame * A[i] * B[j]);

				MElastic[Number(nemo2 + 1, nemo1 + 1)] += V * (DKoefLame * A[i] * A[j] + XKoefLame * B[i] * B[j]);
				MElastic[Number(nemo2 + 0, nemo1 + 0)] += V * (DKoefLame * B[i] * B[j] + XKoefLame * A[i] * A[j]);
				MElastic[Number(nemo2 + 0, nemo1 + 1)] += V * (EKoefLame * A[i] * B[j] + XKoefLame * A[j] * B[i]);
				MElastic[Number(nemo2 + 1, nemo1 + 0)] += V * (EKoefLame * A[j] * B[i] + XKoefLame * A[i] * B[j]);
			}

			MMomentum[Number(nemo1 + 0, nemo2 + 0)] += Rho * V / 12.0;
			MMomentum[Number(nemo1 + 1, nemo2 + 1)] += Rho * V / 12.0;
		}
	}
} //end of BuildMatrix

//consideration of the boundary conditions
void Fixup(void)
{
	int i, k, j, n = 2 * num_assemblies;

	for (i = 0; i < num_fixed_assemblies; i++)
	{
		k = 2 * fixed_assemblies[i] + 0;

		for (j = 0; j < n; j++)
		{
			K1[n * j + k] = 0.0;
			K1[n * k + j] = 0.0;
		}

		K1[n * k + k] = 1.0;
		q[k] = 0.0;
		k = 2 * fixed_assemblies[i] + 1;

		for (j = 0; j < n; j++)
		{
			K1[n * j + k] = 0.0;
			K1[n * k + j] = 0.0;
		}

		K1[n * k + k] = 1.0;
		q[k] = 0.0;
	}
}

//Gaussian elimination
void Gauss(void)
{
	int i, j, k, n = 2 * num_assemblies;
	double coef;

	for (i = 0; i < n; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			coef = K1[n * j + i] / K1[n * i + i];
			for (k = 0; k < n; k++)
			{
				K1[n * j + k] -= K1[n * i + k] * coef;
			}
			q[j] -= q[i] * coef;
		}
	}

	for (i = n - 1; i >= 0; i--)
	{
		for (j = i + 1; j < n; j++)
		{
			q[i] -= K1[n * i + j] * q[j];
		}
		q[i] /= K1[n * i + i];
	}
}

int main(int argc, char *argv[])
{
	FILE *fin, *fout;
	int i, j;
	const double T = 0.4; //time step

	if ((fin = fopen(argv[1], "rt")) == NULL)
	{
		printf("Input file can not be opened");
		exit(1);
	}

	//reading of the number of the nodes
	fscanf(fin, "%d", &num_assemblies);
	//array of the nodes coordinates
	assemblies = new ASSEMBLY[num_assemblies];
	//reading of the nodes coordinates
	for (i = 0; i < num_assemblies; i++)
	{
		fscanf(fin, "%f%f", &assemblies[i].x, &assemblies[i].y);
	}

	//reading of the number of the triangles
	fscanf(fin, "%d", &num_triangles);
	//array of the triangles tops
	triangles = new TRIANGLE[num_triangles];
	//reading of the triangles tops
	for (i = 0; i < num_triangles; i++)
	{
		fscanf(fin, "%d%d%d", &triangles[i].first_top, &triangles[i].second_top, &triangles[i].third_top);
	}

	//reading of the number of the fixed nodes
	fscanf(fin, "%d", &num_fixed_assemblies);
	//array of the fixed nodes
	fixed_assemblies = new int[num_fixed_assemblies];
	//reading of the fixed nodes
	for (i = 0; i < num_fixed_assemblies; i++)
	{
		fscanf(fin, "%d", &fixed_assemblies[i]);
	}

	MElastic = (double *)malloc(2 * num_assemblies * 2 * num_assemblies * sizeof(double));
	MMomentum = new double[2 * num_assemblies * 2 * num_assemblies];
	Right = new double[2 * num_assemblies];
	q = new double[2 * num_assemblies];
	qp = new double[2 * num_assemblies];
	diffq = new double[2 * num_assemblies];
	a = new double[2 * num_assemblies];
	b = new double[2 * num_assemblies];
	K1 = (double *)malloc(2 * num_assemblies * 2 * num_assemblies * sizeof(double));
	K2 = (double *)malloc(2 * num_assemblies * 2 * num_assemblies * sizeof(double));

	//reading of the stress tensor
	for (i = 0; i < num_assemblies; i++)
	{
		fscanf(fin, "%lf%lf", &Right[2 * i + 0], &Right[2 * i + 1]);
	}

	for (i = 0; i < 2 * num_assemblies; i++)
	{
		q[i] = 0.0;
		qp[i] = 0.0;
		diffq[i] = 0.0;
	}

	//openning of the output file
	if ((fout = fopen(argv[2], "wt")) == NULL)
	{
		printf("Output file can not be opened");
		exit(1);
	}

	//the stiffness and the inertia matrices are built
	BuildMatrix();
	fprintf(fout, "<html><body><table border=3>");

	for (double t = 0; t < 2; t += T)
	{

		//combinations of the stiffness and the inertia matrices
		for (i = 0; i < 2 * num_assemblies * 2 * num_assemblies; i++)
		{
			K1[i] = MMomentum[i] / T + MElastic[i] * T / 6.0;
			K2[i] = MMomentum[i] / T - MElastic[i] * T / 3.0;
		}

		for (i = 0; i < 2 * num_assemblies; i++)
		{
			a[i] = 0.0;
			b[i] = 0.0;
			qp[i] = q[i];
		}

		for (i = 0; i < 2 * num_assemblies; i++)
		{
			for (j = 0; j < 2 * num_assemblies; j++)
			{
				a[i] += K2[2 * num_assemblies * i + j] * q[j];
				b[i] += MMomentum[2 * num_assemblies * i + j] * diffq[j];
			}
			q[i] = a[i] + b[i] + Right[i] * T / 2.0;
		}

		Fixup();
		Gauss();

		for (i = 0; i < 2 * num_assemblies; i++)
		{
			diffq[i] = (q[i] - qp[i]) / T;
		}

		//Results output in html file
		fprintf(fout, "<tr><td colspan = 3>Time: %f - %f</td></tr>", t, t + T);

		///The bottom of the bar
		fprintf(fout, "<tr><td>Node</td><td>Bottom Shift X</td><td>Bottom Shift Y</td></tr>");

		fprintf(fout, "<tr><td> ");
		for (i = 0; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%d</br>", i);
		}
		fprintf(fout, "</td>");

		fprintf(fout, "<td>");
		for (i = 0; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%lf</br>", q[2 * i + 0]);
		}
		fprintf(fout, "</td>");

		fprintf(fout, "<td>");
		for (i = 0; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%lf</br>", q[2 * i + 1]);
		}
		fprintf(fout, "</td></tr>");

		///The middle of the bar
		fprintf(fout, "<tr><td>Node</td><td>Middle Shift X</td><td>Middle Shift Y</td></tr>");

		fprintf(fout, "<tr><td> ");
		for (i = 1; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%d</br>", i);
		}
		fprintf(fout, "</td>");

		fprintf(fout, "<td>");
		for (i = 1; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%lf</br>", q[2 * i + 0]);
		}
		fprintf(fout, "</td>");

		fprintf(fout, "<td>");
		for (i = 1; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%lf</br>", q[2 * i + 1]);
		}
		fprintf(fout, "</td></tr>");

		///The top of the bar
		fprintf(fout, "<tr><td>Node</td><td>Top Shift X</td><td>Top Shift Y</td></tr>");

		fprintf(fout, "<tr><td> ");
		for (i = 2; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%d</br>", i);
		}
		fprintf(fout, "</td>");

		fprintf(fout, "<td>");
		for (i = 2; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%lf</br>", q[2 * i + 0]);
		}
		fprintf(fout, "</td>");

		fprintf(fout, "<td>");
		for (i = 2; i < num_assemblies; i += 3)
		{
			fprintf(fout, "%lf</br>", q[2 * i + 1]);
		}
		fprintf(fout, "</td></tr>");

	} //time step iteration

	fprintf(fout, "</table></body></html>");

	delete[] assemblies;
	delete[] triangles;
	delete[] fixed_assemblies;
	delete[] Right;
	delete[] q;
	delete[] qp;
	delete[] diffq;
	delete[] a;
	delete[] b;
	free(MElastic);
	free(MMomentum);
	free(K1);
	free(K2);
	fclose(fin);
	fclose(fout);
	return 0;
}