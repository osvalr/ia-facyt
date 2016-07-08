set xrange [0:1]
set yrange [0:1]
set grid
set key off
set terminal png font ubuntu 14 size 800, 600
set title "Serie Logística Original"
set output "serie_logistica_original.png"
plot "serie_logistica.orig" with lines
set title "Serie Logística Interrogatorio"
set output "serie_logistica_interrogatorio.png"
plot "serie_logistica.interrogatorio" with lines
set title "Serie Logística Entrenamiento"
set output "serie_logistica_entrenamiento.png"
plot "serie_logistica.entrenamiento" with lines
set title "Original/Entrenamiento/Interrogatorio"
set output "serie_logistica_mix_lines.png"
plot "serie_logistica.interrogatorio" with lines, "serie_logistica.orig" with lines, "serie_logistica.entrenamiento" with lines
set output "serie_logistica_mix_dld.png"
plot "serie_logistica.interrogatorio" with points, "serie_logistica.orig" with lines, "serie_logistica.entrenamiento" with points
