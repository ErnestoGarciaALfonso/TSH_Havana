#!/bin/bash

#***************************************************************|
#******************Ernesto Garcia Alfonso***********************|
#***        Facultad de Fisica, Universidad de La Habana    ****|
#                                                               |
#    Programa hecho para automatizar todos mis calculos         |
#        sistema de estudio NeBr2 con TSH                       |
#                       2020                                    |
#                                                               |
#***************************************************************|


FCOMPILER=/home/ernesto/opt/intel/bin/icpc # direccion del compilador 

$FCOMPILER rescalado.cpp # compilas el programa principal que se llama rescalado.cpp 


for ((i=20;i<21;i++)) 
do 
echo "corrida $i"
./a.out ../CI_automaticas/21/esque1-21.int 21 $i >salida 
done 
#pongo 27 porque en Liouville.h tengo para que me imprima las poblacionde de 12 con perdida de un cuanto
