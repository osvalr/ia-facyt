
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <signal.h>
# include <getopt.h>
# include <ctype.h>
# include <stdint.h>
# include <math.h>

#ifdef DIM
#  define PARQUE_DIMENSION DIM                /* Dimension del parque definido en tiempo de compilacion*/
#else
#  define PARQUE_DIMENSION 32                 /* Dimension predeterminada */
#endif

enum {
	SIN_TIPO,
	NEUMANN,
	MOORE
};

# define NUM_FLAGS 5
enum {
	VECINDAD,
	ITER,
	PASOS,
	PPI,
	ARGS,
	DEBUG,
	AYUDA
};

typedef int8_t bool;
enum {
	false,
	true
};

typedef uint8_t estado_t;
typedef uint8_t tiempo_t;

typedef struct {
	estado_t e;
	tiempo_t t;
} area_t;

typedef area_t **parque_t;

enum
{
	FUEGO,
	TERRENO_VACIO,
	PASTO,
	ARBUSTOS,
	BOSQUE
};

# define AREA_INICIAR_INCENDIO { FUEGO, 0 } /* Inicializador estatico para un area en incendio */ 

const area_t AREA_ESTADOS_INICIALES[] =     /* Diferentes tipos de estados iniciales a tomar por area */
{
	{ FUEGO,         0 },
	{ TERRENO_VACIO, 1 },
	{ PASTO,         3 },
	{ ARBUSTOS,      7 },
	{ BOSQUE,       52 }
}; 

/* Lineas y Columnas de la terminal */
static int lineas                = 0;
static int columnas              = 0;

/* Definicion de colores para cada uno de los estados */
const char COLOR_TERRENO_VACIO[] = "\e[37;1m"; /* Gris */
const char COLOR_FUEGO[]         = "\e[31;1m"; /* Rojo */
const char COLOR_ARBUSTOS[]      = "\e[32;1m"; /* Verde */
const char COLOR_BOSQUE[]        = "\e[36;1m"; /* Cyan */
const char COLOR_PASTO[]         = "\e[33;1m"; /* Marron/Amarillo */
const char COLOR_NORMAL[]        = "\e[0m";    /* Normal */

static bool *flags               = NULL;     /* Flags array para los parametros del automata celular */
static int num_iteraciones       = 0;        /* Numero de iteraciones */
static int intervalo_muestreo    = 1;        /* Intervalo de muestreo del estado del array bidimensional */
static double ppi                = 0.0;      /* Probabilidad de propagación del incendio en un area en el vecindario */
static int tipo_vecindad         = SIN_TIPO; /* Tipo de vecindario usado: Moore o Von Neumann */
static int dimension_vecindad    = 0;        /* Tamaño del vecindario, este valor depende directamente del tipo de vecindad */

void mostrar_ayuda ( char *arg )
{
	puts("\n ------------------------------------------------------------");
	puts("| IA - Autómatas Celulares, Marzo 2012                       |");
	puts(" ------------------------------------------------------------");
	puts("| Desarrollado por Luis G. <luisg123v@gmail.com>             |");
	puts("|                  Osval R. <oreyes2@uc.edu.ve>              |");
	puts(" ------------------------------------------------------------");
	printf("\n Uso: %s [ OPCIONES ]\n", arg);
	puts("\nOpciones :");
	puts("\t-t <MOORE|NEUMANN>: Tipo de vecindad");
	puts("\t-p <valor>: PPI: Probabilidad de Propagación del Incendio");
	puts("\t-i <valor>: Cantidad de iteraciones");
	puts("\t-m <valor>: Mostrar cierta cantidad de iteraciones");
	puts("\t-v: Modo depuracion");
	puts("\t-h: Muestra esta ayuda\n");
}

inline bool vecindad_existe_uno_en_llamas( area_t *v, int n )
{
	int i;
	for ( i = 0; i < n && v[i].e != FUEGO; ++i )
		;
	return i != n;
}

void vecindad_liberar( area_t *v, int dimension_vecindad )
{
	if ( v )
	{
		printf("[!] Liberando vecindad..");
		free( v );
		v = NULL;
		puts("OK!");
	}
}

inline bool probabilidad_incendio()
{
	return ppi >= rand() % 100;
}

bool * flags_inicializar ()
{
	int i;
	flags = malloc( NUM_FLAGS );
	for ( i = 0; i < NUM_FLAGS; ++i )
		flags[i] = false;
	return flags;
}

bool flags_verificar()
{
	int i;
	for ( i = 0; i < NUM_FLAGS && !flags[i]; ++i )
		;
	return i == NUM_FLAGS;
}

inline void flags_liberar ( bool *f )
{
	if ( f )
	{
		printf("[!] Liberando banderas..");
		free( f );
		f = NULL;
		puts("OK!");
	}
}

inline area_t area_inicializar()
{
	return AREA_ESTADOS_INICIALES[ rand() % BOSQUE + TERRENO_VACIO ];
}

inline bool area_es_vegetacion( area_t *a )
{
	return a->e == PASTO || a->e == ARBUSTOS || a->e == BOSQUE; 
}

area_t area_transicion( area_t *a, area_t *v, int dimension_vecindad )
{
	area_t t = *a;
	if ( area_es_vegetacion(a) 
		&& vecindad_existe_uno_en_llamas( v, dimension_vecindad ) 
		&& probabilidad_incendio() )
	{
		t = AREA_ESTADOS_INICIALES[FUEGO];
	}
	else
	{
		if ( t.e == FUEGO && t.t == 0 )
		{
			t = AREA_ESTADOS_INICIALES[TERRENO_VACIO];
			return t;
		}

		if ( t.e == TERRENO_VACIO && t.t == 1 )
		{
			t.t++;
			return t ;
		}

		if ( t.e == TERRENO_VACIO && t.t == 2 )
		{
			t = AREA_ESTADOS_INICIALES[PASTO];
			return t;
		}

		if ( t.e == PASTO && t.t == 6 )
		{
			t = AREA_ESTADOS_INICIALES[ARBUSTOS];
			return t;
		}

		if ( t.e == ARBUSTOS && t.t == 51 )
		{
			t = AREA_ESTADOS_INICIALES[BOSQUE];
			return t;
		}
		t.t++;
	}
	return t;
}

parque_t parque_crear()
{
	int i;
	parque_t p = NULL;
	p = malloc ( PARQUE_DIMENSION * sizeof( area_t * ) );
	for ( i = 0; i < PARQUE_DIMENSION; ++i )
		p[i] = malloc( PARQUE_DIMENSION * sizeof( area_t ) );
	return p;
}

void parque_inicializar ( parque_t p )
{
	int i, j;
	for ( i = 0; i < PARQUE_DIMENSION; ++i )
		for ( j = 0; j < PARQUE_DIMENSION; ++j )
			p[i][j] = area_inicializar();
}

void parque_mostrar ( parque_t p )
{
	int i, j;
	if ( flags[DEBUG] )
	{
		system("clear");
		puts("\n");
	}
	for ( i = 0; i < PARQUE_DIMENSION; ++i )
	{
		if ( flags[DEBUG] )
			printf("\t");
		else
			printf("\e[%d;%df", lineas + i, columnas); /* Ubica la matriz en el centro de la terminal */
		for ( j = 0; j < PARQUE_DIMENSION; ++j )
		{
			switch ( p[i][j].e )
			{
				case FUEGO:
					printf(COLOR_FUEGO);
					break;
				case TERRENO_VACIO:
					printf(COLOR_TERRENO_VACIO);
					break;
				case PASTO:
					printf(COLOR_PASTO);
					break;
				case ARBUSTOS:
					printf(COLOR_ARBUSTOS);
					break;
				case BOSQUE:
					printf(COLOR_BOSQUE);
					break;
			}
			printf("█");
			putchar( ( j + 1 ) < PARQUE_DIMENSION? ' ' : '\n');
		}
		printf(COLOR_NORMAL);
	}
}

void parque_incendiar_areas( parque_t p, int *numero_incendios )
{
	int i, j, k, celda;
	area_t incendio = AREA_INICIAR_INCENDIO;
	*numero_incendios = ceil( fmod( rand(), ppi ) ) + 1;
	if ( flags[DEBUG] )
		printf("Incendios: { ");
	
	for ( k = 0; k < *numero_incendios; ++k )
	{
		celda = rand() % (PARQUE_DIMENSION * PARQUE_DIMENSION);
		i = floor( (double) celda / PARQUE_DIMENSION ); 
		j = celda % PARQUE_DIMENSION - 1;
		if ( flags[DEBUG] )
			printf("[%2d][%2d]%s", i + 1, j + 1, (k + 1 == *numero_incendios? " }\n\n" : ", "));
		p[i][j] = incendio;
	}
}

void parque_liberar ( parque_t p )
{
	int i;
	if ( p )
	{
		printf("[!] Liberando parque..");
		for ( i = 0; i < PARQUE_DIMENSION; ++i )
			free( p[i] );
		p = NULL;
		puts("OK!");
	}
}

void signal_handler( int s )
{
	printf("[!] Restaurando terminal..");
	system("stty cooked echo");
	puts("OK\n[!] Saliendo..\n");
	exit( s );
}

void dimension_terminal()
{
	int f, c;
	char tmp[20], *termtype = getenv("TERM");
	if ( tgetent(tmp, termtype) >= 0 )
	{
		f = tgetnum("li");
		c = tgetnum("co");
	}

	if ( f < PARQUE_DIMENSION  )
	{
		sprintf(tmp, "resize -s %d %d > /dev/null", f = PARQUE_DIMENSION, c);
		system(tmp);
	}
	if ( c < PARQUE_DIMENSION  )
	{
		sprintf(tmp, "resize -s %d %d > /dev/null", f, c = PARQUE_DIMENSION);
		system(tmp);
	}
	lineas = f / 2 - PARQUE_DIMENSION / 2;
	columnas = c / 2 - PARQUE_DIMENSION;
}

void getkey ( void )
{
	system("stty cbreak -echo");
	fflush(stdout);
	getchar();
	system("stty cooked echo");
}

int main( int argc, char **argv )
{
	static parque_t p                = NULL;
	static parque_t p_sig            = NULL;
	static area_t *vecindad          = NULL;
	int i, j, t;
	char opcion;
	int numero_incendios;
	signal( SIGINT, signal_handler ); /* Manejo de señales */
	flags = flags_inicializar();
	srand( getpid() );

	if ( argc < 8 )
	{
		puts(" - Insuficientes argumentos");
		mostrar_ayuda(argv[0]);
		return EXIT_FAILURE;
	}

	while ( ( opcion = getopt( argc, argv, "t:i:p:m:vh") ) != -1 && flags_verificar() )
	{
		switch ( opcion )
		{
			case 't':
				for ( i = 0; i < strlen( optarg ) ; ++i )
					optarg[i] = tolower( optarg[i] );
				tipo_vecindad = !strcmp("moore", optarg)? MOORE : ( !strcmp("neumann", optarg)? NEUMANN : SIN_TIPO ); 
				dimension_vecindad = tipo_vecindad == MOORE ? 8 : 4;
				if ( !tipo_vecindad )
					flags[VECINDAD] = true;
				break;

			case 'i':
				if ( ( num_iteraciones = atoi( optarg ) ) <= 0 )
					flags[ITER] = true; 
				break;

			case 'p':
				if ( ( ppi = atof( optarg ) ) <= 0.0 )
					flags[PPI] = true;
				break;

			case 'm':
				if ( ( intervalo_muestreo = atoi( optarg ) ) <= 0 )
					flags[PASOS] = true;
				break;

			case '?':
				flags[ARGS] = true;
				break;

			case 'h':
				flags[AYUDA] = true;
				break;

			case 'v':
				flags[DEBUG] = true;
				break;

			case ':':
				flags[ARGS] = true;
				break;
		}
	}

	if ( !flags_verificar() || flags[AYUDA] )
	{
		mostrar_ayuda(argv[0]);
		if ( !flags[AYUDA] )
		{
			if ( flags[VECINDAD] )
				puts("El tipo de vecindad introducido no es válido!");
			if ( flags[ITER] )
				puts("Revise el numero de iteraciones");
			if ( flags[PASOS] )
				puts("Revise el numero de pasos");
			if ( flags[PPI] )
				puts("Revise la Probabilidad de Propagación del Incendio ( a.k.a PPI )");
			if ( flags[ARGS] )
				puts("Revise los argumentos de las opciones usadas");
		}

		return EXIT_FAILURE;
	}

	dimension_terminal();

	p = parque_crear();
	p_sig = parque_crear();

	vecindad = malloc( dimension_vecindad * sizeof( area_t ) );
	parque_inicializar( p );
	parque_incendiar_areas( p, &numero_incendios );

	if ( flags[DEBUG] )
	{
		puts("Parametros de ejecucion");
		printf(" - Tipo de vecindad      : %s\n", ( tipo_vecindad == MOORE? "Moore" : "Neumann" ) );
		printf(" - Numero de iteraciones : %d\n", num_iteraciones );
		printf(" - Numeros de pasos      : %d\n", intervalo_muestreo );
		printf(" - Numeros de Incendios  : %d\n\n", numero_incendios);
	}

	printf("\033[2J \e[%d;0f", lineas);
	parque_mostrar( p ); /* Muestra estado inicial */
	if ( flags[DEBUG] )
		getkey();
	else 
		usleep(150000);
	for ( t = 0; t < num_iteraciones; ++t )
	{
		for ( i = 0; i < PARQUE_DIMENSION; ++i )
		{
			for ( j = 0; j < PARQUE_DIMENSION; ++j )
			{
				if ( tipo_vecindad == NEUMANN )
				{
					vecindad[0] = p[(PARQUE_DIMENSION + (i - 1)) % PARQUE_DIMENSION][j];
					vecindad[1] = p[i][(PARQUE_DIMENSION + (j - 1)) % PARQUE_DIMENSION];
					vecindad[2] = p[i][(j + 1) % PARQUE_DIMENSION];
					vecindad[3] = p[(i + 1) % PARQUE_DIMENSION][j];
				}
				else /* Vecindad de Moore */
				{
					vecindad[0] = p[(PARQUE_DIMENSION + (i - 1)) % PARQUE_DIMENSION][(PARQUE_DIMENSION + (j + 1)) % PARQUE_DIMENSION];
					vecindad[1] = p[(PARQUE_DIMENSION + (i - 1)) % PARQUE_DIMENSION][j];
					vecindad[2] = p[(PARQUE_DIMENSION + (i - 1)) % PARQUE_DIMENSION][(j + 1) % PARQUE_DIMENSION];
					vecindad[3] = p[i][(PARQUE_DIMENSION + (j - 1)) % PARQUE_DIMENSION];
					vecindad[4] = p[i][(j + 1) % PARQUE_DIMENSION];
					vecindad[5] = p[(i + 1) % PARQUE_DIMENSION][(PARQUE_DIMENSION  + (j - 1)) % PARQUE_DIMENSION];
					vecindad[6] = p[(i + 1) % PARQUE_DIMENSION][j];
					vecindad[7] = p[(i + 1) % PARQUE_DIMENSION][(j + 1) % PARQUE_DIMENSION];
				}
				p_sig[i][j] = area_transicion( &p[i][j], vecindad, dimension_vecindad);
			}
		}

		for ( i = 0; i < PARQUE_DIMENSION; ++i )
			for ( j = 0; j < PARQUE_DIMENSION; ++j )
				p[i][j] = p_sig[i][j];

		if ( (t + 1) % intervalo_muestreo == 0 )
		{
			parque_mostrar( p );
			if ( flags[DEBUG] )
				getkey();
			else
				usleep(150000);
		}
	}

	vecindad_liberar ( vecindad, dimension_vecindad );
	parque_liberar( p );
	parque_liberar( p_sig );
	flags_liberar( flags );
	return EXIT_SUCCESS;
}
