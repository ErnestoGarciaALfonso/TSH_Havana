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


source /share/apps/intel/2015.1/bin/compilervars.sh intel64
#lo de arriba para que no me de problema con el compilador nuevo

FCOMPILER=icpc # direccion del compilador 

$FCOMPILER rescalado.cpp # compilas el programa principal que se llama rescalado.cpp 


continuar=$1


if(($continuar==1))
then

############################################## Step 1 Making folders
cd ../
mkdir v=0 # se crea la carpeta v=0
cd v=0  # estamos dentro de la carpeta v=0

 # creamos las carpetas de v=16,.....,29
for ((i=16;i<26;i++)) 
do 
mkdir v$i
done


##################################################### Step 2  Copies and links
cd ../ # estoy en el directorio principal

ls
cd v=0 #entramos a v=0

for ((i=16;i<26;i++)) 
do 
cd v$i
cp ../../main/fonda* .  #copio todas las funciones de onda
cp ../../main/final .
cp ../../main/Para_someter/qsubmit . #script donde ejecuto el programa principal
cp ../../main/final .
cp ../../main/parametros_iniciales .
ln -s -f ../../main/a.out compila$i.out 
#./compila$i.out ../../main/v=0/esque1-$i.int $i $j  >salida$i # ejecuta 
cd ../
done



################################################# Step 3 Making results folders

cd ../ # estoy en el directorio principal 
mkdir Resultadosv0
cd Resultadosv0

for((i=16;i<26;i++))
do
mkdir resultadosv$i
done

################################################ Step 4 qsub qsubmit

cd ../

cd v=0 #estos son reset coef con metodo de Truhlar
# pone en background todos los programas
for((i=16;i<26;i++)) 
do
cd v$i
qsub -N ArBr2nv$i qsubmit ../../main/v=0/esque1-$i.int $i ../../Resultadosv0/resultadosv$i 
cd ../

done


else


if(($continuar==2))
then

cd ../ # estoy en el directorio principal 

cd v=0
# pone en background todos los programas
for((i=16;i<26;i++)) 
do

cd v$i
qsub -N ArBr2nv$i qsubmit ../../main/v=0/esque1-$i.int $i ../../Resultadosv0/resultadosv$i 
cd ../

done

else

echo "******************************************************************"
echo "* ./ejecutar.sh 1 es para iniciar los programas                  *"
echo "* ./ejecutar.sh 2 es para continuar los programas                *" 
echo "******************************************************************"

fi

fi



######################################################### Step 6 The end

echo "Todos los programas estan en background......"
echo "Gracias por ejecutarme!!!!!!, Have a nice day "




 

