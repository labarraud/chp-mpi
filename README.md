# Projet de Calcul Haute Performance

By [Clément DELIBES](https://github.com/Cle-md) and  [Laurent Barraud](https://github.com/labarraud)  

## Description
Projet de parallélisation de schéma numérique pour la résolution de l'équation de la chaleur en 2D au moyen de la bibliothèque open MPI.  Projet réalisé au cours de l'option CHP (Calcul Hautes Performance) sous la supervision de Héloïse Beaugendre en 2eme Annee de l'ENSEIRB-MATMECA.

## Compilation
Avant de compiler les sources, indiquer où se trouve le chemin du compilateur `mpif90` dans la première ligne du [Makefile] (https://github.com/labarraud/chp-mpi/blob/master/Makefile), en initialisant la variable `FC`

```makefile
FC=mpif90
```
Ensuite, compiler les sources en exécutant la commande `make` en étant positionné à la racine du répertoire contenant le projet. 
```console
$> make
```

Cette commande compilera tous les fichiers sources contenus dans le dossier [Src] (https://github.com/labarraud/chp-mpi/blob/master/Src), créra les fichiers objets et les modules dans leurs dossiers respectifs [Object] (https://github.com/labarraud/chp-mpi/blob/master/Object) et [Module] (https://github.com/labarraud/chp-mpi/blob/master/Module), et l'exécutable `chppar`

## Fichier de configuration

Dans le fichier [param.txt] (https://github.com/labarraud/chp-mpi/blob/master/param.txt), indiquer selon la norme décrite ci-dessous les valeurs des variables globales du problème, le chemin relatif du dossier où seront exportés les fichiers des solutions calculés par chaque processeur et le cas test choisi. 

Par exemple, si l'on souhaite initialiser le programme à `epsilon=0.00001` (condition d'arrêt du gradiant conjugué), `Lx=1.0`, `Ly=1.0`, `Nx=10`, `Ny=10`, `D=1.0`, `dt=1.0`, `Nt=10`, pour le cas test 2 et que l'on souhaite que les fichiers solutions soient exportés dans le dossier `Result`, le fichier de configuration est le suivant :  

```txt
0.00001		!epsilon
1.0	    	!Lx
1.0	    	!Ly
10    		!Nx
10    		!Ny
1.    		!D
1.     		!dt
10    		!Nt
Result	  !folder result
2       	!cas
```

## Exécution

Une fois le fichier [param.txt] (https://github.com/labarraud/chp-mpi/blob/master/param.txt) comfiguré, on exécute le programme au moyen de la commande `mpirun`. Par exemple, si l'on souhaite exécuter le programme sur 4 processeurs on tapera la commande suivante :

```console
$> mpirun -n 4 ./chpmpi
```

## Visualisation des résultats

Les fichiers résultats `.dat`  sont créés dans le dossier `Result` pour visualiser le résultat sous la forme d'un graphique. Il suffit alors d'excécuter le script gnuplot `plotresult.gnu` créé lors de l'exécution du programme.

```console
$> gnuplot plotresult.gnu
```

Cela créé une image eps. A titre d'exemple, et en gardant les paramètres explicités plus haut, le résultat est le suivant :

![myimage-alt-tag](http://labarraud.vvv.enseirb-matmeca.fr/chp/plotresult.png)

## Calcul de la charge, du speedup et de l'éfficacité

Avant d'excuter le script shell `computestat.sh`, il faut paramétrer le nombre maximum de processeurs pour lesquels on veut calculer le speed up et la charge, la variable `max` ainsi que le chemin de l'exécutable `mpirun`. De plus, on peut modifier le nombre de fois que sera exécuté le programme pour un nombre de processeurs donné avec la variable `nbmoy`. Ces executions répétées sont faites pour se soustraire au xdifférentes interférences causées par les processus externes qui ont tendance à interrompre l'exécution des processus du programme et augmenter le temps d'exécution. En prenant le minimum de la charge maximun pour un nombre de processeurs donné on s'approche au plus près de la charge maximum réelle.

Pour un `max = 20`, `nbmax=10` et un chemin de l'executatble mpirun `mpirun` on aura comme début de script :

```shell
#!/bin/bash

max=20
nbmoy=10
mpirun="mpirun"

(...)
```

Ensuite on donnera au script l'autorisation de s'exécuter par la commande :

```console
$> chmod +x computestat.sh
```

Enfin on executera le script :

```console
$> ./computestat.sh
```

Les résultats de la charge (temps d'exécution), du speedup et de l'éfficacité se trouveront respectivement dans les fichiers créés `time.dat`, `speedup.dat` et `efficacite.dat`
