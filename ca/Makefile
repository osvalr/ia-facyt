#!/usr/bin/make
GCC=gcc
CFLAGS=-g -ltermcap -lm
POBJETO=bin/main
FUENTE=source/main.c
ifdef DIM
	CFLAGS+= -D DIM=$(DIM) 
endif
ARGS=$(FUENTE) $(CFLAGS) -o $(POBJETO)
make :
	@$(GCC) $(ARGS)
	@bash -c  \
	' echo -n "Compilando..."; \
	if [ -f bin/main ]; then \
		echo " Listo! "; \
	fi; \
	echo -n "Generando enlace simbólico..."; \
	if [ ! -f main ]; then \
		ln -s bin/main main; \
		echo " Listo! "; \
	else \
		echo " Ya existe! "; \
	fi '
clean :
	@bash -c  \
	' echo -n "Eliminando enlace simbólico..."; \
	if [ -f main ]; then \
		unlink main; \
		echo " Listo! "; \
	else \
		echo " No se hizo nada "; \
	fi; \
	echo -n "Eliminando binarios..."; \
	if [ -f bin/main ]; then \
		echo " Listo! "; \
		rm bin/main; \
	else \
		echo " No se hizo nada "; \
	fi '

