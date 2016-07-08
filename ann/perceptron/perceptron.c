# include <stdio.h>
# include <unistd.h>
# include <stdint.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

enum {
	false,
	true
};

typedef uint8_t bool;

enum { // Tipos de entrenamiento
	NO_TYPE = 0,
	BATCH_TYPE,
	PATTERN_TYPE
};

enum { // Tipos de error
	NO_ERROR = 0,
	NO_VALOR,
	NO_ITERACIONES,
	NO_TIPO,
	NO_OPCION,
	ERROR_RUTA,
	ERROR_MALLOC
};

# define UMBRAL -1

int MAX_PATTERNS = 0;               // Cantidad de Patrones de Procesar
int MAX_UNIDADES_PROCESAMIENTO = 0; // Cantidad de Unidades de Procesamiento
int MAX_DENDRITAS = 0;              // Cantidad de Dentritas por neurona
bool modo_depuracion = false;

void set_fpu ( unsigned int mode )
{
	asm( "fldcw %0" : : "m" (*&mode) );
}
void mostrar_headers()
{
	int i;
	printf("| _______ | ________________ | __________ Salida _________ |\n|  Patron | ");
	for ( i = 1; i <= MAX_DENDRITAS; ++i )
		printf("   X%d%s  ", i , i == MAX_DENDRITAS - 1 ? " " : "  " );
	puts("|Deseada  |Calculada|Error Abs|");
}
void ayuda( char *arg )
{
	printf("\n Uso: %s [ OPCIONES ]\n", arg);
	puts("\nOpciones :");
	puts("\t-l [λ]: Parámetro de aprendizaje");
	puts("\t-i [It]: Numero de iteraciones");
	puts("\t-t <Patron/Lotes>: Tipo de entrenamiento: \"Patron\" ó \"Lotes\"");
	puts("\t-f [file.sample]: Archivo de entrada");
	puts("\t-p [P]: Numero de patrones a procesar");
	puts("\t-P [pasos]: Paso para mostrar avances ( valor mayor que 1 )");
	puts("\t-u [N]: Numero de unidades de procesamiento");
	puts("\t-d [M]: Numero de dendritas por neurona");
	puts("\t-v: Modo depuracion");
	puts("\t-h: Muestra esta ayuda\n");
}

void init_matriz_pesos( double **w )
{
	int i, j;
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
	{
		for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
		{
			w[i][j] = ( rand() % 6000 ) / 10000.0 - 0.3;
			if ( modo_depuracion )
				printf("%3.5lf ", w[i][j]);
		}
		if ( modo_depuracion )
			puts("");
	}
}

void init_matriz( double **dw )
{
	int i, j;
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
		for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
			dw[i][j] = 0.0;
}

int funcion_transferencia( double x ) // Heaviside
{
	if ( modo_depuracion )
		printf("x=%lf\n", x);
	return  x > 0.0 ? 1 : 0;

}

double actividad_neuronal( int **X, int fila_k, double **W, int col_r )
{
	int i;
	double r = 0;
	if ( modo_depuracion )
	{
		puts("Entrando Actividad Neuronal");
		printf("k= %d col_r= %d\n", fila_k, col_r);
		puts("    r     X[fila_k][i]     W[i][col_r]");
	}
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
	{
		r = r + X[fila_k][i] * W[i][col_r];
		if ( modo_depuracion )
			printf("%c%0.5lf%9d           %c%2.5lf\n", r < 0.0? 0 : ' ', r, X[fila_k][i], W[i][col_r] < 0.0? 0 : ' ', W[i][col_r]);
	}
	if ( modo_depuracion )
	puts("");
	return r;
}

void actualizar_pesos( double **W, double **dW )
{
	int i, j;
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
		for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
			W[i][j] = W[i][j] + dW[i][j];
}

int main(int argc, char **argv)
{
	int i, iter, j, k, r, opcion;
	double **W = NULL, **dW = NULL;
	int **X = NULL, **Y = NULL, **unidad_procesamiento = NULL;
	int error, pasos = -1;
	bool continuar;
	double lambda, dWri, error_acumulado;
	int8_t tipo_entrenamiento;
	int iteraciones;
	bool mostrar_ayuda = false;
	bool interrogar = false;
	int8_t error_flag = NO_ERROR;
	FILE *f = NULL;
	srand( getpid() ); // Implante de semilla de aletoriedad
	set_fpu (0x27F);  // use double-precision rounding

	while ( ( opcion = getopt (argc, argv, "l:i:t:p:P:f:u:d:vhI") ) != -1  &&  ( !error_flag || !mostrar_ayuda ) )
	{
		switch ( opcion )
		{
			case 'P': // el programa será interrogado después de el entrenamiento
				pasos = atoi( optarg );
				pasos = pasos <= 1? -1 : pasos;
				break;
			case 'I': // el programa será interrogado después de el entrenamiento
				interrogar = true;
				break;

			case 'l': // Parametro Lambda
				lambda = atof( optarg );
				break;

			case 'i': // Numero de iteraciones
				iteraciones = atoi( optarg );
				if ( iteraciones <= 0 )
					error_flag = NO_ITERACIONES;
				break;

			case 't': // Parametro que indica si la evaluacion es por patron o por lotes
				for ( i = 0; i < strlen( optarg ); ++i )
					optarg[i] = tolower( optarg[i] );
				if ( !( tipo_entrenamiento = ( !strcmp( "patron", optarg ) )? PATTERN_TYPE :  ( !strcmp( "lotes", optarg ) )? BATCH_TYPE : NO_TYPE ) )
					error_flag = NO_TIPO;
				break;

			case 'v': // Coloca al programa en modo verbose
				modo_depuracion = true;
				break;

			case 'h': // Muestra la ayuda
				mostrar_ayuda = true;
				break;

			case 'p':
			 	MAX_PATTERNS = atoi( optarg );
				break;

			case 'f':
				if ( ( f = fopen( optarg, "r" ) ) == NULL  )
				{
					error_flag = ERROR_RUTA;
					break;
				}
				fscanf(f, "%d %d %d", &MAX_PATTERNS, &MAX_UNIDADES_PROCESAMIENTO, &MAX_DENDRITAS );
				if ( modo_depuracion )
				{
					printf("P=%d N=%d M=%d\n", MAX_PATTERNS, MAX_UNIDADES_PROCESAMIENTO, MAX_DENDRITAS );
					printf("[INIT]Direcciones de memoria de matrices @ %d : %p   \t%p   \t%p   \t%p   \t%p   \n", __LINE__, unidad_procesamiento, W, dW, X, Y);
				}

				if ( ( unidad_procesamiento = malloc( sizeof( int * ) * MAX_PATTERNS ) ) == NULL )
					puts("\n\t - ERROR: No se pudo reservar memoria para las unidades de procesamiento");

				if ( ( W = malloc( sizeof( double * ) * ( MAX_DENDRITAS + 1 ) ) ) == NULL ) 
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz de Pesos [W]");

				if ( ( dW = malloc( sizeof( double * ) * ( MAX_DENDRITAS + 1 ) ) ) == NULL )
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz auxiliar de Entrenamiento por Lotes [dW]");

				if ( ( X = malloc( sizeof( int * ) * MAX_PATTERNS ) ) == NULL )
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz para las componentes de entrada de los patrones [X]");

				if ( ( Y = malloc( sizeof( int * ) * MAX_PATTERNS ) ) == NULL )
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz para las salidas deseadas [Y]");

				if ( !W || !dW || !X || !Y || !unidad_procesamiento )
				{
					error_flag = ERROR_MALLOC;
					break;
				}

				for ( i = 0; i <= MAX_DENDRITAS; ++i )
				{
					W[i] = malloc( sizeof( double ) * MAX_UNIDADES_PROCESAMIENTO );
					for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Inicializacion de matriz de pesos
						W[i][j] = i? 0.0 : UMBRAL;
				}
				
				for ( i = 0; i <= MAX_DENDRITAS; ++i )
				{
					dW[i] = malloc( sizeof( double ) * MAX_UNIDADES_PROCESAMIENTO );
					for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Inicializacion de la matriz Auxiliar para entrenamiento por lotes
						dW[i][j] = i? 0.0 : UMBRAL;
				}

				for ( i = 0; i < MAX_PATTERNS; ++i )
				{
					X[i] = malloc( sizeof( int ) *  ( MAX_DENDRITAS + 1 ) );
					for ( j = 0; j <= MAX_DENDRITAS; ++j ) // Inicializacion de las Componentes de Entrada X,
						X[i][j] = j? 0 : UMBRAL;
				}

				for ( i = 0; i < MAX_PATTERNS; ++i )
				{
					Y[i] = malloc( sizeof( int ) *  MAX_UNIDADES_PROCESAMIENTO );
					unidad_procesamiento[i] = malloc( sizeof( int ) *  MAX_UNIDADES_PROCESAMIENTO );
					for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
						unidad_procesamiento[i][j] = 0;
				}

				if ( modo_depuracion )
				{
					puts("Matrix W");
					for ( i = 0; i <= MAX_DENDRITAS; ++i )
					{
						for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
							printf("%c%2.2lf ", W[i][j] < 0.0? 0 : ' ', W[i][j]);
						puts("");
					}
					puts("Matrix dW");
					for ( i = 0; i <= MAX_DENDRITAS; ++i )
					{
						for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
							printf("%c%2.2lf ",dW[i][j] < 0.0? 0 : ' ', dW[i][j]);
						puts("");
					}
					puts("Matrix X");
					for ( i = 0; i < MAX_PATTERNS; ++i )
					{
						for ( j = 0; j <= MAX_DENDRITAS; ++j )
							printf("%d ", X[i][j]);
						puts("");
					}					
					puts("Matrix Y");
					for ( i = 0; i < MAX_PATTERNS; ++i )
					{
						for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
							printf("%d ", Y[i][j]);
						puts("");
					}
					
					printf("[RESV]Direcciones de memoria de matrices @ %d : %p   \t%p   \t%p   \t%p   \t%p   \n", __LINE__, unidad_procesamiento, W, dW, X, Y);
					for ( i = 0; i < MAX_DENDRITAS; ++i )
						printf("X%d ", i + 1 );
					for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
						printf("Y%d ", j + 1 );
					puts("");
				}
				
				for ( i = 0; i < MAX_PATTERNS; ++i )
				{
					for ( j = 1; j < MAX_DENDRITAS + 1; ++j )
					{
						fscanf( f, "%d", &X[i][j] );
						if ( modo_depuracion )
							printf( "%d  ", X[i][j] );
					}
					for ( k = 0; k < MAX_UNIDADES_PROCESAMIENTO; ++k )
					{
						fscanf( f, "%d", &Y[i][k] );
						if ( modo_depuracion )
							printf( "%d  ", Y[i][k] );
					}
					if ( modo_depuracion )
						puts("");
				}
				fclose ( f );
				break;

			case 'u':
				MAX_UNIDADES_PROCESAMIENTO = atoi( optarg );
				break;

			case 'd':
				MAX_DENDRITAS = atoi( optarg );
				break;

			case ':':
				error_flag = NO_VALOR;
				break;

			case '?':
				error_flag = NO_OPCION;
				break;
		}
		if ( error_flag )
			break;
	}

	if ( argc < 5 || mostrar_ayuda )
	{
		ayuda( argv[0] );
		return EXIT_SUCCESS;
	}

	if ( error_flag )
	{
		ayuda( argv[0] );
		puts("\n\t- ERROR: ");
		switch ( error_flag )
		{
			case NO_VALOR:
				printf("'-%c' requiere un valor\n", optopt );
				break;

			case NO_ITERACIONES:
				puts("Debe indicar la cantidad de iteraciones con un numero entero mayor a 0\n");
				break;

			case NO_TIPO:
				puts("Debe definir un tipo de entrenamiento\n");
				break;

			case NO_OPCION: // No se reconoce la opcion
				printf("No se reconoce la opcion '-%c'\n", optopt );
				break;
		}
		return EXIT_FAILURE;
	}

	if( f == NULL ) // si no se especificó algún archivo de entrada
		f = stdin;

	init_matriz_pesos( W );
	if ( modo_depuracion )
	{
		printf(">> Semilla = %5d\n", getpid() );
		printf(">> Iteraciones = %d\n", iteraciones);
		printf(">> Constante de aprendizaje = %lf\n", lambda);
		printf(">> Entrenamiento: %s\n", tipo_entrenamiento == PATTERN_TYPE? "Patron a patron" : "Por lotes" );
	}
	for ( iter = 0, continuar = true; continuar && iter < iteraciones ; ++iter )
	{
		continuar = false;
		if ( tipo_entrenamiento == BATCH_TYPE )
			init_matriz(dW);
		for ( k = 0; k < MAX_PATTERNS; ++k )
		{
			for ( r = 0; r < MAX_UNIDADES_PROCESAMIENTO; ++r )
			{
				unidad_procesamiento[k][r] = funcion_transferencia( actividad_neuronal( X, k, W, r ));
				error = Y[k][r] - unidad_procesamiento[k][r];
				if ( error != 0 )
				{
					continuar = true;
					if ( modo_depuracion )
						printf("continuar = false\n Y[%d][%d] - unidad_procesamiento[%d] = %d - %d \n", k, r, r, Y[k][r],  unidad_procesamiento[k][r] );
				}
				for ( i = 0; i <= MAX_DENDRITAS; ++i )
				{
					dWri = X[k][i] * error * lambda;
					if ( tipo_entrenamiento == PATTERN_TYPE ) // Patron a Patron
						W[i][r] = W[i][r] + dWri;
					else // Por lotes
						dW[i][r] = dW[i][r] + dWri;
				}
			}
		}
		// Mostrar tabla
		if ( !continuar || ( pasos != -1 && ( iter == 0 || !(( iter + 1 ) % pasos )) ) )
		{
			if ( !continuar || iter == iteraciones - 1 )
			{
				printf("** Ultima Iteracion ( Por : %s )**\n", continuar? "Nro de Iteraciones" : "Bandera" );
				printf("Nro. de Iteraciones: %5d\n", iteraciones );
				printf("Constante de Aprendizaje: %lf\n", lambda );
				printf("Tipo de Entrenamiento: %s\n", tipo_entrenamiento == PATTERN_TYPE? "Patron a patron" : "Por lotes" );
			}
			printf("Iteracion Actual:%5d\n", iter + 1);
			mostrar_headers();
			for ( k = 0; k < MAX_PATTERNS; ++k )
			{
				printf("| %7d |  ", k + 1);
				for ( j = 1; j <= MAX_DENDRITAS; ++j ) // Muestra la salida deseada para cada neurona
					printf("  %d     ", X[k][j]);
				putchar('|');
				for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra la salida deseada para cada neurona
					printf("    %d    ", Y[k][j]);
				putchar('|');
				for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra la salida obtenida para cada neurona
					printf("    %d    ", unidad_procesamiento[k][j]);
				putchar('|');
				for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra el valor absoluto de las salidas de las neuronas
					printf("    %d    ", abs( Y[k][j] - unidad_procesamiento[k][j] ) );
				puts("|");
			}
			puts("");
		}
		if ( tipo_entrenamiento == BATCH_TYPE )
			actualizar_pesos( W, dW );
		if ( modo_depuracion )
			printf("iter vale %d, continuar vale %s, iteraciones vale %d\n", iter + 1, continuar? "true": "false", iteraciones);
	}

	// Liberación de memoria
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
	{
		free( W[i] );
		free( dW[i] );
	}
	free( W );
	free( dW );

	for ( i = 0; i < MAX_PATTERNS; ++i )
	{
		free( X[i] );
		free( Y[i] );
	}
	free( X );
	free( Y );
	for ( i = 0; i < MAX_PATTERNS; ++i )
		free( unidad_procesamiento[i] );
	free( unidad_procesamiento );
	return EXIT_SUCCESS;
}
