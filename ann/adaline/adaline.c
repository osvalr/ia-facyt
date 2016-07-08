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

enum { // Tipo de entrenamiento
	NO_TYPE = 0,
	BATCH_TYPE,
	PATTERN_TYPE
};

enum { // Tipo de funcion de transferencia
	IDENTIDAD = 1,
	SIGMOIDE
};

enum { // ERRORES
	NO_ERROR = 0,
	NO_VALOR,
	NO_ITERACIONES,
	NO_TIPO,
	NO_FUNCION,
	NO_OPCION,
	ERROR_RUTA,
	ERROR_MALLOC,
};

# define UMBRAL -1.0
int MAX_PATTERNS = 0;               // Cantidad de Patrones de Procesar
int MAX_UNIDADES_PROCESAMIENTO = 0; // Cantidad de Unidades de Procesamiento
int MAX_DENDRITAS = 0;              // Cantidad de Dentritas por neurona
bool modo_depuracion = false;
bool mostrar_datos_grafica = false;
int8_t tipo_funcion;

void set_fpu ( unsigned int mode )
{
	asm( "fldcw %0" : : "m" (*&mode) );
}

void ayuda( char *arg )
{
	printf("\n Uso: %s [ OPCIONES ]\n", arg);
	puts("\n Opciones :");
	puts("\t-l [η]: Parámetro de aprendizaje");
	puts("\t-i [I]: Numero de iteraciones");
	puts("\t-t <Patron/Lotes>: Tipo de entrenamiento: \"Patron\" ó \"Lotes\"");
	puts("\t-f [file.sample]: Archivo de entrada");
	puts("\t-F <Identidad/Sigmoide>: Tipo de funcion de transferencia: \"Identidad\" ó \"Sigmoide\"");
	puts("\t-e [ϵ]: epsilon, por defecto 0.1");
	puts("\t-p [P]: Numero de patrones a procesar");
	puts("\t-u [N]: Numero de unidades de procesamiento");
	puts("\t-d [M]: Numero de dendritas por neurona");
	puts("\t-D: Muestra solo datos necesarios para la representacion de los datos en una grafica ( oculta la opcion -v )");
	puts("\t-I: Interrogar luego de entrenar, las entradas se tomaran desde el archivo de entrada especificado con -F");
	puts("\t-P [pasos]: Paso para mostrar avances ( valor mayor que 1 )");
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
			if ( modo_depuracion && !mostrar_datos_grafica )
				printf("%3.5lf ", w[i][j]);
		}
		if ( modo_depuracion && !mostrar_datos_grafica )
			puts("");
	}
}

void mostrar_headers()
{
	int i;
	printf("| _______ | ____________ | _______ Salida _______ |\n|  Patron | ");
	for ( i = 1; i <= MAX_DENDRITAS; ++i )
		printf("   X%d%s  ", i , i == MAX_DENDRITAS - 1 ? " " : "  " );
	puts("|Deseada  |Calculada|Error Abs|");
}

void init_matriz( double **dw )
{
	int i, j;
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
		for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
			dw[i][j] = 0.0;
}

double funcion_transferencia( double x ) // Heaviside
{
	if ( modo_depuracion  && !mostrar_datos_grafica )
		printf("Funcion de transferencia\n\tx=%lf\n", x);
	return tipo_funcion == SIGMOIDE? 1 / ( 1.0 + exp( -x ) ) : x; 
}

double actividad_neuronal( double *X, double **W, int col_r )
{
	int i;
	double r = 0;
	if ( modo_depuracion && !mostrar_datos_grafica )
	{
		puts("Entrando Actividad Neuronal");
		puts("    r     X[fila_k][i]     W[i][col_r]");
	}
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
	{
		r = r + X[i] * W[i][col_r];
		if ( modo_depuracion && !mostrar_datos_grafica )
			printf("%0.5lf     %2.5lf    %2.5lf\n", r, X[i], W[i][col_r]);
	}
	if ( modo_depuracion && !mostrar_datos_grafica )
		puts("\nSaliendo Actividad Neuronal");

	return r;
}

void actualizar_pesos( double **W, double **dW )
{
	int i, j;
	for ( i = 0; i <= MAX_DENDRITAS; ++i )
		for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
			W[i][j] = W[i][j] + dW[i][j];
}

double dF( double x )
{
	return tipo_funcion == IDENTIDAD? 1.0 : funcion_transferencia( -x ) * ( 1.0 - funcion_transferencia( -x ) );
}

int main(int argc, char **argv)
{
	int i, iter, j, k, r, opcion;
	double **W = NULL, **dW = NULL;
	double **X = NULL, **Y = NULL, **unidad_procesamiento = NULL;
	double error, error_acumulado;
	bool continuar;
	double eta, dWri, epsilon = 0.1;
	int8_t tipo_entrenamiento;
	int iteraciones, pasos = -1;
	bool mostrar_ayuda = false;
	bool interrogar = false;
	int8_t error_flag = NO_ERROR;
	FILE *f = NULL;

	srand( getpid() ); // Implante de semilla de aletoriedad
	set_fpu (0x27F);  // use double-precision rounding

	while ( ( opcion = getopt (argc, argv, "d:De:f:F:hi:Il:p:P:t:u:v") ) != -1  &&  ( !error_flag || !mostrar_ayuda ) )
	{
		switch ( opcion )
		{
			case 'P': // el programa será interrogado después de el entrenamiento
				pasos = atoi( optarg );
				pasos = pasos <= 1? -1 : pasos;
				break;

			case 'e': // epsilon
				epsilon = atof(optarg) == 0.0 ? epsilon : atof(optarg);
				break;

			case 'D': // Muestra solo la salida de los datos post-entrenamiento, post-interrogatorios, necesarios para representar una funcion
				mostrar_datos_grafica = true;
				break;

			case 'I': // el programa será interrogado después de el entrenamiento
				interrogar = true;
				break;

			case 'l': // Parametro eta
				eta = atof( optarg );
				break;

			case 'i': // Numero de iteraciones
				iteraciones = atoi( optarg );
				if ( iteraciones <= 0 )
					error_flag = NO_ITERACIONES;
				break;

			case 'F': // Parametro que indica la funcion de transferencia a usar
				for ( i = 0; i < strlen( optarg ); ++i )
					optarg[i] = tolower( optarg[i] );
				if ( !( tipo_funcion = ( !strcmp( "identidad", optarg ) )? IDENTIDAD :  ( !strcmp( "sigmoide", optarg ) )? SIGMOIDE : NO_TYPE ) )
					error_flag = NO_FUNCION;
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
				if ( modo_depuracion && !mostrar_datos_grafica )
				{
					printf(">> P=%d N=%d M=%d\n", MAX_PATTERNS, MAX_UNIDADES_PROCESAMIENTO, MAX_DENDRITAS );
					printf(">> [INIT]Direcciones de memoria de matrices @ %d : %p   \t%p   \t%p   \t%p   \t%p   \n", __LINE__, unidad_procesamiento, W, dW, X, Y);
				}

				if ( ( unidad_procesamiento = malloc( sizeof( double * ) *  MAX_PATTERNS ) ) == NULL )
				{
					puts("\n\t - ERROR: No se pudo reservar memoria para las unidades de procesamiento");
					return EXIT_FAILURE;
				}

				for( i = 0; i < MAX_PATTERNS; ++i)
				{
					unidad_procesamiento[i] = malloc( sizeof(double) * MAX_UNIDADES_PROCESAMIENTO );
					for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
						unidad_procesamiento[i][j] = 0.0;
				}
				if ( ( W = malloc( sizeof( double * ) * ( MAX_DENDRITAS + 1 ) ) ) == NULL ) 
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz de Pesos [W]");

				if ( ( dW = malloc( sizeof( double * ) * ( MAX_DENDRITAS + 1 ) ) ) == NULL )
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz auxiliar de Entrenamiento por Lotes [dW]");

				if ( ( X = malloc( sizeof( double * ) * MAX_PATTERNS ) ) == NULL )
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz para las componentes de entrada de los patrones [X]");

				if ( ( Y = malloc( sizeof( double * ) * MAX_PATTERNS ) ) == NULL )
					puts("\n\t - ERROR: No se pudo reservar memoria para la matriz para las salidas deseadas [Y]");

				if ( !W || !dW || !X || !Y )
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
					X[i] = malloc( sizeof( double ) *  ( MAX_DENDRITAS + 1 ) );
					for ( j = 0; j <= MAX_DENDRITAS; ++j ) // Inicializacion de las Componentes de Entrada X,
						X[i][j] = j? 0 : UMBRAL;
				}

				for ( i = 0; i < MAX_PATTERNS; ++i )
					Y[i] = malloc( sizeof( double ) *  MAX_UNIDADES_PROCESAMIENTO );

				if ( modo_depuracion && !mostrar_datos_grafica )
				{
					puts(">> Matrix W");
					for ( i = 0; i <= MAX_DENDRITAS; ++i )
					{
						for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
							printf("%2.2lf ", W[i][j]);
						puts("");
					}
					puts(">> Matrix dW");
					for ( i = 0; i <= MAX_DENDRITAS; ++i )
					{
						for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
							printf("%2.2lf ",dW[i][j]);
						puts("");
					}
					puts(">> Matrix X");
					for ( i = 0; i < MAX_PATTERNS; ++i )
					{
						for ( j = 0; j <= MAX_DENDRITAS; ++j )
							printf("%.2lf ", X[i][j]);
						puts("");
					}					
					puts(">> Matrix Y");
					for ( i = 0; i < MAX_PATTERNS; ++i )
					{
						for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
							printf("%.2lf ", Y[i][j]);
						puts("");
					}
					printf(">> [RESV]Direcciones de memoria de matrices @ %d : %p   \t%p   \t%p   \t%p   \t%p   \n", __LINE__, unidad_procesamiento, W, dW, X, Y);
					for ( i = 0; i < MAX_DENDRITAS; ++i )
						printf("X%d ", i + 1 );
					for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j )
						printf("Y%d ", j + 1 );
					puts("");
				}

				for ( i = 0; i < MAX_PATTERNS; ++i )
				{
					if ( modo_depuracion && !mostrar_datos_grafica )
						printf(">> ");
					for ( j = 1; j < MAX_DENDRITAS + 1; ++j )
					{
						fscanf( f, "%lf", &X[i][j] );
						if ( modo_depuracion && !mostrar_datos_grafica )
							printf( "%.2lf  ", X[i][j] );
					}
					if ( modo_depuracion && !mostrar_datos_grafica )
						printf("\n>> ");
					for ( k = 0; k < MAX_UNIDADES_PROCESAMIENTO; ++k )
					{
						fscanf( f, "%lf", &Y[i][k] );
						if ( modo_depuracion && !mostrar_datos_grafica )
							printf( "%.2lf  ", Y[i][k] );
					}
					if ( modo_depuracion && !mostrar_datos_grafica )
						puts("\n");
				}
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

			case NO_FUNCION: // No se reconoce la opcion
				printf("Debe definir una funcion de transferencia");
				break;
		}
		return EXIT_FAILURE;
	}

	if ( f == NULL ) // si no se especificó algún archivo de entrada
		f = stdin;

	init_matriz_pesos( W );
	if ( modo_depuracion && !mostrar_datos_grafica )
	{
		printf(">> Semilla = %5d\n", getpid() );
		printf(">> Iteraciones = %d\n", iteraciones);
		printf(">> Constante de aprendizaje = %lf\n", eta);
		printf(">> Entrenamiento: %s\n", tipo_entrenamiento == PATTERN_TYPE? "Patron a patron" : "Por lotes" );
		printf(">> Funcion de transferencia: %s\n", tipo_funcion == SIGMOIDE? "Sigmoide" : "Identidad" );
	}

	for ( iter = 0, continuar = true; continuar && iter < iteraciones ; ++iter )
	{
		continuar = false;
		if ( tipo_entrenamiento == BATCH_TYPE )
			init_matriz(dW);
		for ( k = 0; k < MAX_PATTERNS; ++k )
		{
			error_acumulado = 0.0;
			for ( r = 0; r < MAX_UNIDADES_PROCESAMIENTO; ++r )
			{
				unidad_procesamiento[k][r] = funcion_transferencia( actividad_neuronal( X[k], W, r ));
				error = Y[k][r] - unidad_procesamiento[k][r];
				error_acumulado += fabs(error);
				if ( fabs(error) > epsilon  )
				{
					continuar = true;
					if ( modo_depuracion && !mostrar_datos_grafica )
						printf(">> continuar = true\n>> Y[%d][%d] - unidad_procesamiento[%d][%d] = %lf - %lf \n", k, r, k, r, Y[k][r],  unidad_procesamiento[k][r] );
				}
				for ( i = 0; i <= MAX_DENDRITAS; ++i )
				{
					dWri = eta * error * dF( unidad_procesamiento[k][r] ) * X[k][i] ;
					if ( tipo_entrenamiento == PATTERN_TYPE ) // Patron a Patron
						W[i][r] = W[i][r] + dWri;
					else // Por lotes
						dW[i][r] = dW[i][r] + dWri;
				}
			}
			error_acumulado /= MAX_UNIDADES_PROCESAMIENTO;
		}
		// Mostrar tabla
		if ( !continuar || ( pasos != -1 && ( iter == 0 || !(( iter + 1 ) % pasos )) ) )
		{
			if ( !continuar || iter == iteraciones - 1 )
			{
				printf("** Ultima Iteracion ( Por : %s )**\n", continuar? "Nro de Iteraciones" : "Bandera" );
				printf("Nro. de Iteraciones: %5d\n", iteraciones );
				printf("Constante de Aprendizaje: %lf\n", eta );
				printf("Tipo de Entrenamiento: %s\n", tipo_entrenamiento == PATTERN_TYPE? "Patron a patron" : "Por lotes" );
				printf("Funcion de transferencia: %s\n", tipo_funcion == SIGMOIDE? "Sigmoide" : "Identidad" );
				printf("Epsilon: %lf\n", epsilon );
			}
			printf("Iteracion Actual:%5d\n", iter + 1);
			printf("Error promedio: %.8lf\n", error_acumulado);
			mostrar_headers();
			for ( k = 0; k < MAX_PATTERNS; ++k )
			{
				printf("| %7d |", k + 1);
				for ( j = 1; j <= MAX_DENDRITAS; ++j ) // Muestra la salida deseada para cada neurona
					printf("%lf ", X[k][j]);
				putchar('|');
				for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra la salida deseada para cada neurona
					printf("%lf ", Y[k][j]);
				putchar('|');
				for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra la salida obtenida para cada neurona
					printf("%lf ", unidad_procesamiento[k][j]);
				putchar('|');
				for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra el valor absoluto de las salidas de las neuronas
					printf("%lf ", fabs( Y[k][j] - unidad_procesamiento[k][j] ) );
				puts("|");
			}
			puts("");
		}
		if ( tipo_entrenamiento == BATCH_TYPE )
			actualizar_pesos( W, dW );
		if ( modo_depuracion && !mostrar_datos_grafica )
			printf(">> iter = %d\n>> continuar vale %s\n>> iteraciones vale %d\n\n", iter + 1, continuar? "true": "false", iteraciones);
		
	}
	if ( mostrar_datos_grafica )
	{
		printf(">> Patrones Entrenados: %5d\n", MAX_PATTERNS);
		epsilon = 1.0 / ( MAX_PATTERNS + 50 );
		printf(">> Delta: %lf\n", epsilon);
		for ( j = 0, eta = epsilon; j < MAX_PATTERNS + 50; ++j, eta += epsilon )
			printf("%lf %lf\n", eta, j < MAX_PATTERNS ? unidad_procesamiento[j][0] : 0.0 ); // Solo se considera la salida de la primera neurona (tos) CAIMAN DEL ORINOCO
	}
	// se interroga la red neuronal si así se especificó en la entrada
	if( interrogar )
	{
		mostrar_headers();
		for( i = 1; !feof( f ); ++i )
		{
			for( j = 1; j <= MAX_DENDRITAS; ++j ) // Se carga el patrón a interrogar
				fscanf(f, "%lf", &X[0][j] );

			for( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Se carga la salidas deseadas para el patron
				fscanf(f, "%lf ", &Y[0][j] );

			for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Se interroga a la red neuronal por el patron
				unidad_procesamiento[i][j] = funcion_transferencia ( actividad_neuronal( X[0], W, j ) );

			printf("| %7d |", i);
			for ( j = 1; j <= MAX_DENDRITAS; ++j ) // Muestra las componentes de entrada de un patron
				printf("%lf ", X[0][j]);
			putchar('|');
			for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra la salida deseada para cada neurona
				printf("%lf ", Y[0][j]);
			putchar('|');
			for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra la salida obtenida para cada neurona
				printf("%lf ", unidad_procesamiento[i][j]);
			putchar('|');
			for ( j = 0; j < MAX_UNIDADES_PROCESAMIENTO; ++j ) // Muestra el valor absoluto de las salidas de las neuronas
				printf("%lf ", fabs( Y[0][j] - unidad_procesamiento[i][j] ) );
			puts("|");
		}
		if ( mostrar_datos_grafica )
		{
			printf(">> Patrones interrogados: %5d\n", --i);
			epsilon = 1.0 / ( MAX_PATTERNS + i );
			printf(">> Delta: %lf\n", epsilon);
			for ( j = 0, eta = epsilon; j < MAX_PATTERNS + i; ++j, eta += epsilon )
				printf("%lf %lf\n", eta, j < MAX_PATTERNS ? 0.0 : unidad_procesamiento[ j - MAX_PATTERNS ][0]); // Solo se considera la salida de la primera neurona (tos) CAIMAN DEL ORINOCO pt. 1.2
		}
	}

	if ( f != stdin )
		fclose( f );

	for ( i = 0; i <= MAX_DENDRITAS; ++i ) // Liberación de memoria
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
