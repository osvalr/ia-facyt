#include <stdio.h>
#include <stdlib.h>

# define MAX_UNIDADES_PROCESAMIENTO 1
# define MAX_ENTRENAMIENTO 148
# define MAX_COMPONENTES 5
# define MAX_SALIDAS 1

void set_fpu ( unsigned int mode )
{
	asm( "fldcw %0" : : "m" (*&mode) );
}

double serie_logistica( double x )
{
	return 4.0 * x  * ( 1.0 - x );
}

int main(int argc, char **argv)
{
	set_fpu(0x27F);
	double X0 = 0.0, x = 0.0, delta = 0.0;
	double X[200];
	double R[198][3];
	double vectores[200][6];
	int i, j;

	while ( ( X0 = ( rand() % 9999 + 1 ) / 10000.0 ) == 0.5 ) // Busqueda de un valor aleatorio distinto de 1/2
		;

	for ( i = 0; i < 500; ++i ) // 500 iteraciones para superar el transiente
		x = i? serie_logistica( x ) : serie_logistica( X0 );

	
	for ( i = 0, delta = 1 / 200.0; i < 200; ++i, delta += 1 / 200.0 ) // Generacion de 200 valores
		X[i] = x = serie_logistica( x );

	for ( i = 0, delta = 1 / 198.0; i < 198; ++i, delta += 1 / 198.0 ) // Pre-procesamiento polinomico
	{
		printf("%lf %lf\n", delta, X[i] );
		vectores[i][0] = X[i];
		vectores[i][1] = X[i + 1];
		vectores[i][2] = X[i] * X[i];
		vectores[i][3] = X[i] * X[i + 1];
		vectores[i][4] = X[i + 1] * X[i + 1];
		vectores[i][5] = X[i + 2];
	}

	printf("%d %d %d\n", MAX_ENTRENAMIENTO, MAX_SALIDAS, MAX_COMPONENTES);
	for ( i = 0; i < 198; ++i )
	{
		for ( j = 0; j < 6; ++j )
			printf("%lf ", vectores[i][j]);
		puts("");
	}
	return 0;
}
