#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <signal.h>

# define MAX_GENES 31
# define MAX_GENOTIPO 0x7FFFFFFF
# define putbit( bitarray , pos ) putchar( bitarray>>( MAX_GENES - 1 - i ) & 1 ? '1' : '0' )

typedef uint8_t bool;
enum {
	false,
	true
};

typedef uint32_t cromosoma_t;
typedef cromosoma_t *poblacion_t;
typedef double fenotipo_t;

typedef struct
{
	double probabilidad;
	int indice;
} adaptacion;

poblacion_t p                = NULL;
adaptacion *a                = NULL;
bool *flags                  = NULL;
int MAX_ITERACIONES          = 0;
int MAX_INDIVIDUOS           = 0;
int GAP_GENERACIONAL         = 0;
double PROBABILIDAD_MUTACION = 0.0;
double PROBABILIDAD_CRUCE    = 0.0;

#define MAX_FLAGS 9 /* Cantidad de banderas que depende de la cantidad de opciones en el menu */

enum {
	NO_ITER,
	NO_POBL,
	NO_PM,
	NO_PC,
	NO_GAP,
	NO_VALOR,
	NO_OPCION,
	DEBUG,
	AYUDA
};

void set_fpu ( unsigned int mode )
{
	asm( "fldcw %0" : : "m" (*&mode) );
}

double f( double x ) // función a minimizar
{
	return 0.45 * ( 2.0 - x ) + ( 1.0 - fabs( sin( 32.0 * x ) ) );
}

int cromosoma_punto_cruce()
{
	return rand() % (MAX_GENES - 1) + 1;
}

bool es_mutable()
{
	return ( PROBABILIDAD_MUTACION >= rand() / (double) RAND_MAX );
}

bool existe_cruce()
{
	return ( PROBABILIDAD_CRUCE >= rand() / (double) RAND_MAX );
}

fenotipo_t decodificar ( cromosoma_t c )
{
	return ( c * 1.3 ) / MAX_GENOTIPO;
}

cromosoma_t codificar ( fenotipo_t f )
{
	return ( f * MAX_GENOTIPO ) / 1.3;
}

double funcion_adaptacion( double y )
{
	return 1.0 / ( 1.0 + y );
}

void poblacion_liberar( poblacion_t p )
{
	free( p );
	p = NULL;
}

void signal_handler( int s )
{
	puts("[!] Preparando para salir..");
	printf("[!] Liberando memoria..");
	poblacion_liberar ( p );
	free( a );
	printf("OK\n[!] Restaurando terminal..");
	system("stty cooked echo");
	puts("OK\n[!] Saliendo..\n");
	exit( s );
}

void getkey ( char *arg )
{ 
	puts ( arg );
	system("stty cbreak -echo");
	fflush(stdout);
	getchar();
	system("stty cooked echo");
}

poblacion_t poblacion_crear ()
{
	int i, j;
	poblacion_t p = NULL;
	cromosoma_t individuo = 0;

	if( flags[DEBUG] )
		printf("Creando población: %p > ", p );

	if ( !( p = malloc( MAX_INDIVIDUOS * sizeof(cromosoma_t) ) ) )
		return;

	if( flags[DEBUG] )
		printf("%p\n", p );

	for ( i = 0; i < MAX_INDIVIDUOS; ++i )
		do
		{
			individuo = rand() % MAX_GENOTIPO;
			if ( flags[DEBUG] )
				printf("p[%d]: %10u - %lf\n", i, individuo , decodificar( individuo ) );
			p[i] = individuo;
			for ( j = 0; j < i && p[j] != p[i]; ++j )
				;
		}while( j != i );

	if( flags[DEBUG] )
		puts("\n");

	return p; // se devuelve la población creada
}

void cromosoma_mostrar( cromosoma_t c )
{
	int i;
	for ( i = 0; i < MAX_GENES; ++i )
		putbit( c, i );
}

int comparar_adaptacion(const void *a, const void *b)
{
	if ( ((adaptacion *)a)->probabilidad < ((adaptacion *)b)->probabilidad)
		return 1;
	if ( ((adaptacion *)a)->probabilidad > ((adaptacion *)b)->probabilidad)
		return -1;
	return 0;
}

void cromosoma_mutar( cromosoma_t *c )
{
	int i;
	for ( i = 0; i < MAX_GENES; ++i )
		if ( es_mutable() )
			*c ^= 1<<i;
}

bool cromosoma_existe ( poblacion_t p, int n, cromosoma_t c )
{
	int i;
	for ( i = 0; i < n && p[i] != c; ++i )
		;
	return i != n; 
}

void poblacion_cruzar( poblacion_t padres )
{
	cromosoma_t hijos[MAX_INDIVIDUOS], pi=0, pj=0, hijo1=0, hijo2=0, mascara=0;
	int i, j, hijos_almacenados, punto_cruce;
	adaptacion c[MAX_INDIVIDUOS];
	double sumatoria = 0, ruleta = 0, n, ruleta_aux;

	if( flags[DEBUG] )
	{
		puts("(1/)Presion Ambiental");
		puts("Indice\tProbabilidad\tSumatoria");
	}

	for( i = 0; i < MAX_INDIVIDUOS; ++i ) // 
	{
		hijos[i] = 0;
		c[i].indice = i;
		c[i].probabilidad = funcion_adaptacion( f( decodificar( padres[i] ) ) ); /* Calculo de aptitud de cada cromosoma */
		sumatoria += c[i].probabilidad; /* Calculo de aptitud total */
		if( flags[DEBUG] )
			printf("%4d\t%.9lf\t%.9lf\n", c[i].indice + 1, c[i].probabilidad, sumatoria );
	}

	if( flags[DEBUG] )
		puts("\n(2/)Probabilidad Real");

	for( i = 0; i < MAX_INDIVIDUOS; ++i )
	{
		c[i].probabilidad /= sumatoria; /* Calculando probabilidad de seleccion */
		if( flags[DEBUG] )
			printf("%.9lf\n", c[i].probabilidad);
	}

	if( flags[DEBUG] )
	{
		puts("\nIndice\tProbabilidad");
		for ( i = 0; i < MAX_INDIVIDUOS; ++i )
			printf("%4d\t%.9lf\n", i + 1, c[i].probabilidad );
		puts("\n(3/)Sumas Acumuladas");
	}

	for( i = 1; i < MAX_INDIVIDUOS; ++i ) /* Calculo de probabilidad acumulada */
		c[i].probabilidad += c[i - 1].probabilidad;

	if( flags[DEBUG] )
		for ( i = 0; i < MAX_INDIVIDUOS; ++i )
			printf("%4d\t%.9lf\n", i + 1, c[i].probabilidad);

	hijos_almacenados = ceil( MAX_INDIVIDUOS * ( 100 - GAP_GENERACIONAL ) / 100.0 );

	if( flags[DEBUG] )
		printf("\n(4/)Cantidad de individuos de la poblacion a mantener: %d\nCopiando Padres a la siguiente generacion...", hijos_almacenados);

	for ( i = 0; i < hijos_almacenados; ++i )
		hijos[i] = padres[i];

	if( flags[DEBUG] )
	{
		printf("Listo!\n\n(5/) Entrando a calcular los nuevos %d individuos de la poblacion\n", MAX_INDIVIDUOS - hijos_almacenados);
		getkey("Presione cualquier tecla para continuar...[[Ctrl + C para salir]]");
	}

	n = -1;
	do {

		ruleta_aux = ruleta = rand() / (double)( RAND_MAX - 1 ); /* Rueda la ruleta  P1 */
		for( pi = 0; pi < MAX_INDIVIDUOS - 1; ++pi ) // se busca a cual cromosoma corresponde el valor obtenido
			if ( c[pi].probabilidad <= ruleta && ruleta < c[pi + 1].probabilidad )
				break;

		do {
			ruleta = rand() / (double)( RAND_MAX - 1 ); /* Rueda la ruleta P2 */
			for ( pj = 0; pj < MAX_INDIVIDUOS - 1 ; ++pj)
				if ( c[pj].probabilidad <= ruleta && ruleta < c[pj + 1].probabilidad )
					break;
		} while ( pi == pj );

		if( flags[DEBUG] )
		{
			printf("\n\n\n\nRuleta P1   %.9lf\n", ruleta_aux);
			printf("Ruleta P2   %.9lf\nP1 %3d ----> ", ruleta, pi + 1 );
			cromosoma_mostrar( padres[c[pi].indice] );
			printf(" [%.9lf] \nP2 %3d ----> ", decodificar( padres[pi] ) , pj + 1 );
			cromosoma_mostrar( padres[c[pj].indice] );
			printf(" [%.9lf]\n", decodificar( padres[pj] ) );
		}

		if( flags[DEBUG] )
			printf("Existe cruce? ");

		if ( existe_cruce() )
		{
			if( flags[DEBUG] )
				puts("SI");

			punto_cruce = cromosoma_punto_cruce(); // Se obtiene el punto de cruce

			if( flags[DEBUG] )
				printf("Punto de Cruce: %2d\n", punto_cruce );

			for ( mascara = 0, j = 0; j < punto_cruce; ++j ) // Generación de la mascara
				mascara |= 1<<j;

			if( flags[DEBUG] )
			{
				printf("Mascara:     ");
				cromosoma_mostrar(mascara);
				puts("\n");
			}

			// Sucede el cruce
			if( flags[DEBUG] )
				puts("Cruce:");

			hijo1 = (padres[c[pi].indice] &  mascara) | (padres[c[pj].indice] & ~mascara);
			hijo2 = (padres[c[pi].indice] & ~mascara) | (padres[c[pj].indice] &  mascara);
		}
		else // si no se cruzan, pasan directamente a la nueva generación
		{
			if( flags[DEBUG] )
				puts("NO <<Los padres pasan igual>>");

			hijo1 = padres[c[pi].indice];
			hijo2 = padres[c[pj].indice];
		}

		if( flags[DEBUG] )
		{
			printf("\tH1 : ");
			cromosoma_mostrar(hijo1);
			printf("\n\tH2 : ");
			cromosoma_mostrar(hijo2);
			putchar('\n');
		}

		// Sucede la mutación
		cromosoma_mutar( &hijo1 );
		cromosoma_mutar( &hijo2 );

		if( flags[DEBUG] )
		{
			printf("\nLuego de la mutacion:\n\tH1 : ");
			cromosoma_mostrar(hijo1);
			printf("\n\tH2 : ");
			cromosoma_mostrar(hijo2);
			puts("");
		}

		if ( cromosoma_existe ( hijos, i, hijo1 ) || cromosoma_existe ( hijos, i, hijo2 ) ) /* Se comprueba que ninguno de los hijos existan en la siguiente generación */
		{
			if ( flags[DEBUG] )
				puts("################################################################################\nAlgunos de los dos individuos generados, ya existen en la generacion siguiente\n################################################################################");
			continue;
		}
		if ( i + 1 == MAX_INDIVIDUOS )
		{
			if ( !( n = rand() % 2 ) ) // Hijo1 es quien pasa a la siguiente generacion;
				hijos[i] = hijo1;
			else // Hijo2
				hijos[i] = hijo2;
			i++;
		}
		else
		{
			hijos[i] = hijo1;
			hijos[i + 1] = hijo2;
			i += 2;
		}

		if( flags[DEBUG] )
		{
			printf("\n  - Numero de individuos procesados: %d de %d\n\n", i , MAX_INDIVIDUOS);
			if ( i == MAX_INDIVIDUOS )
			{
				if ( n >= 0 )
					printf("    NOTA: Se debe tener en cuenta que en el ultimo cruce se asignó el Hijo%c\n", n? '2' : '1' );
				puts("\nSe completó el número de individuos para esta generacion\n");
			}
			else
				getkey("Presione cualquier tecla continuar con la siguiente seleccion de padres [[Ctrl + C para salir]]\n");	
			puts("\n\n");
		}
	} while ( i < MAX_INDIVIDUOS );

	for ( i = 0; i < MAX_INDIVIDUOS; ++i ) // se sustituye la generación actual con la nueva
		padres[i] = hijos[i];
}


void mostrar_ayuda ( char *arg )
{
	puts("\n ------------------------------------------------------------");
	puts("| IA - Algoritmos Geneticos, Febrero 2012                    |");
	puts(" ------------------------------------------------------------");
	puts("| Desarrollado por Luis G. <luisg123v@gmail.com>             |");
	puts("|                  Osval R. <oreyes2@uc.edu.ve>              |");
	puts(" ------------------------------------------------------------");
	printf("\n Uso: %s [ OPCIONES ]\n", arg);
	puts("\nOpciones :");
	puts("\t-i <valor>: Cantidad de iteraciones");
	puts("\t-n <valor>: Numero de individuos en la población");
	puts("\t-m <valor>: Probabilidad de mutación");
	puts("\t-c <valor>: Probabilidad de cruce");
	puts("\t-g <valor>: Gap generacional");
	puts("\t-v: Modo depuracion");
	puts("\t-h: Muestra esta ayuda\n");
}

bool *flags_inicializar () /* Reservar espacio e inicializa arreglo de flags */
{
	int i;
	bool *f;
	if ( !( f = malloc( MAX_FLAGS ) ) )
		return NULL;
	for ( i = 0; i < MAX_FLAGS; ++i )
		f[i] = false;

	return f; // se devuelve el arreglo de banderas creado
}

bool flags_verificar ( bool *f ) /* Verifica si al menos un flag está activado, retorna true si no hay ningun problema,  false en caso contrario */
{
	int i;
	for ( i = 0; i < MAX_FLAGS && f && !f[i]; ++i )
		;
	if ( i < MAX_FLAGS - 4 )
		switch ( i )
		{
			case NO_ITER:
				printf("ERROR %d : Debe indicar el numero de iteraciones con el cual se va a ejecutar\n", __LINE__ );
				break;
			case NO_POBL:
				printf("ERROR %d : Debe indicar la cantidad de individuos para la poblacion\n", __LINE__ );
				break;
			case NO_PM:
				printf("ERROR %d : Debe indicar la probabilidad de mutacion\n", __LINE__ );
				break;
			case NO_PC:
				printf("ERROR %d : Debe indicar la probabilidad de cruce\n", __LINE__ );
				break;
			case NO_GAP:
				printf("ERROR %d : GAP Generacional debe ubicarse entre 0-100%%\n", __LINE__ );
				break;
		}
	return i == MAX_FLAGS;
}

double poblacion_adaptacion_promedio( poblacion_t p )
{
	int i;
	double promedio = funcion_adaptacion( f( decodificar( p[0] ) ) );
	for ( i = 1; i < MAX_INDIVIDUOS; ++i )
		promedio += funcion_adaptacion( f( decodificar( p[i] ) ) );
	return promedio / (double)MAX_INDIVIDUOS;
}

int main( int argc, char **argv )
{
	int i, j;
	char opcion;
	double tmp;

	signal( SIGINT, signal_handler ); /* Manejo de señales */
	srand( getpid() );                /* Ajuste de la semilla de aleatoriedad */
	set_fpu (0x27F);                  /* Ajuste de la aritmetica de la FPU con redondeo de doble precision */
	flags = flags_inicializar();      /* Inicializacion de banderas */

	while ( ( opcion = getopt( argc, argv, "i:n:m:c:g:vh") ) != -1 && flags_verificar(flags) )
	{
		switch ( opcion )
		{
			case '?':
				flags[NO_OPCION] = true;
				break;

			case 'h':
				flags[AYUDA] = true;
				break;

			case 'v':
				flags[DEBUG] = true;
				break;

			case ':':
				flags[NO_VALOR] = true;
				break;

			case 'i': // Numero de iteraciones
				if ( ( MAX_ITERACIONES = atoi( optarg ) ) <= 0 )
					flags[NO_ITER] = true;
				break;

			case 'n': // Numeros de individuos en la población
				if ( ( MAX_INDIVIDUOS = atoi( optarg ) ) <= 1 ) // La cantidad de individuos deberá ser mayor a 1 
					flags[NO_POBL] = true;						// para ejecutar al menos 1 cruce entre dos individuos de la poblacion
				break;

			case 'm':
				if ( ( PROBABILIDAD_MUTACION = atof( optarg ) ) == 0.0 )
					flags[NO_PM] = true;
				break;

			case 'c':
				if ( ( PROBABILIDAD_CRUCE = atof( optarg ) ) == 0.0 )
					flags[NO_PC] = true;
				break;

			case 'g':
				if ( ( GAP_GENERACIONAL = atoi( optarg ) ) <= 0 || GAP_GENERACIONAL > 100 )
					flags[NO_GAP] = true;
				break;
		}
	}

	if ( flags[DEBUG] )
	{
		printf("Parámetros usados\n");
		printf("Numero de Individuos: %d\n", MAX_INDIVIDUOS);
		printf("Numero de Genes por Individuos: %d\n", MAX_GENES);
		printf("Probabilidad de Cruce: %lf\n", PROBABILIDAD_CRUCE);
		printf("Probabilidad de Mutacion: %lf\n", PROBABILIDAD_MUTACION);
		printf("Gap generaacional: %d\n\n", GAP_GENERACIONAL);
	}
	else
		if ( !flags_verificar(flags) || argc < 9 || flags[AYUDA] )
		{
			if ( argc != 1 && argc < 9 && !flags[NO_OPCION] )
				puts(" - Insuficientes argumentos para ejecutarse");
			mostrar_ayuda( argv[0] );
			return EXIT_FAILURE;
		}

	if ( !( p = poblacion_crear() ) )
		return EXIT_FAILURE;

	a = malloc( sizeof(adaptacion) * MAX_INDIVIDUOS );

	puts("\n\n  GENERACION INICIAL");
	for ( tmp = 0.0, j = 0; j < MAX_INDIVIDUOS; ++j )
	{
		a[j].probabilidad = funcion_adaptacion ( f ( decodificar ( p[j] ) ) );
		a[j].indice = j;
		tmp += a[j].probabilidad;
	}

	qsort( a, MAX_INDIVIDUOS, sizeof(adaptacion), comparar_adaptacion);

	tmp /= (double) MAX_INDIVIDUOS;
	printf("  Mejor adaptación: %.9lf,    Adaptación promedio: %.9lf\n\n", a[0].probabilidad, tmp );
	puts("  No  Cromosoma                         Valor Real    Evaluación    Adaptación");
	for ( j = 0; j < MAX_INDIVIDUOS; ++j )
	{

		if ( MAX_INDIVIDUOS == 100 && j + 1 == 11 )
			printf("      ...\n");

		if ( MAX_INDIVIDUOS == 10 || ( j < 10 || MAX_INDIVIDUOS - j < 9 ) )
		{
			printf("%4d  ", j + 1 );
			cromosoma_mostrar( p[a[j].indice] );
			printf("   %.9lf   %.9lf   %.9lf\n", decodificar ( p[a[j].indice] ), f(decodificar ( p[a[j].indice] )), a[j].probabilidad );
		}
	}

	for( i = 0; i < MAX_ITERACIONES; ++i )
	{
		poblacion_cruzar(p);
		if ( i == MAX_ITERACIONES - 1 || i == MAX_ITERACIONES / 2 - 1 )
		{
			printf("\n\n  GENERACION ACTUAL: %5d\n", i + 1 );
			for ( tmp = 0.0, j = 0; j < MAX_INDIVIDUOS; ++j )
			{
				a[j].probabilidad = funcion_adaptacion ( f ( decodificar ( p[j] ) ) );
				a[j].indice = j;
				tmp += a[j].probabilidad;
			}

			qsort( a, MAX_INDIVIDUOS, sizeof(adaptacion), comparar_adaptacion);

			tmp /= (double) MAX_INDIVIDUOS;
			printf("  Mejor adaptación: %.9lf,    Adaptación promedio: %.9lf\n\n", a[0].probabilidad, tmp );
			puts("  No  Cromosoma                         Valor Real    Evaluación    Adaptación");
			for ( j = 0; j < MAX_INDIVIDUOS; ++j )
			{
				if ( MAX_INDIVIDUOS == 100 && j + 1 == 11 )
					printf("      ...\n");

				if ( MAX_INDIVIDUOS == 10 || ( j < 10 || MAX_INDIVIDUOS - j < 9 ) )
				{
					printf("%4d  ", j + 1);
					cromosoma_mostrar(p[a[j].indice]);
					printf("   %.9lf   %.9lf   %.9lf\n", decodificar ( p[a[j].indice] ), f(decodificar ( p[a[j].indice] )), a[j].probabilidad);
				}
			}

			if ( flags[DEBUG] )
				getkey( i + 1 == MAX_ITERACIONES ? "Presione cualquier tecla para finalizar..[[Ctrl + C para salir]]" : "\nFinalizó la generación.. Presione para continuar[[Ctrl + C para salir]]");
			printf("\n       Mínimo: f(Xopt) = %.9lf  | Alcanzado en Xopt = %.9lf\n\n", f(decodificar ( p[a[0].indice] )), decodificar ( p[a[0].indice] ));
		}
	}

	poblacion_liberar ( p ); /* liberación de la memoria reservada */
	free(a);
	return EXIT_SUCCESS;
}
