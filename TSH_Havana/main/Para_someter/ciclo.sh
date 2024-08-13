#!/bin/bash


#***************************************************************|
#******************Ernesto Garcia Alfonso***********************|
#***        Facultad de Fisica, Universidad de La Habana    ****|
#                                                               |
# Es importante que la forma que esta abajo se mantenga, en caso|
# que se tenga que usar otro formato de script se copia o se    |
# modifique este de manera que funcione de igual manera         |
#                                                               |
#***************************************************************|




path=$1 # es la direccion, por ejemplo ../../main/v=0/esque1-16.int 
nv=$2 # es el nivel vibracional de partida, por ejemplo 16,....,29
pathresult=$3 # es la direccion de salida donde guardo mis resultados
#echo "path-> "$path, "nivel de partida-> "$nv



#pone los trabajos en background

# el numero 5000 es por la cantidad de condiciones iniciales
for ((i=0;i<5000;i++)) 
do

./compila$nv.out $path $nv $i >salida$nv 

done

cp Tiempo_disociacion.txt Disangular Erotacional Canal_disociacion.txt $pathresult
 
