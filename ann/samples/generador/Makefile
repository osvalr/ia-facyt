make :
	gcc generador.c -o generador

datos :
	./generador > generador.out
	head -198 generador.out >  serie_logistica.orig
	tail -199 generador.out > ../serie_logistica.sample
	rm generador.out
graphic :
	gnuplot < plotall.plt
