# Projet de Calcul Haute Performance

By [Clément DELIBES](https://github.com/Cle-md) and  [Laurent Barraud](https://github.com/labarraud)  

## Description
Projet de parallélisation de schéma numérique pour la résolution de l'équation de la chaleur en 2D au moyen de la bibliothèque open MPI.  Projet réalisé au cours de l'option CHP (Calcul Hautes Performance) sous la supervision de Héloïse Beaugendre en 2eme Annee de l'ENSEIRB-MATMECA.

## Compilation
Avant de compiler les sources, tout d'abord indiqué ou se trouve le chemin du compilateur `mpif90` dans la première ligne du [Makefile] (https://github.com/labarraud/chp-mpi/blob/master/Makefile), en initialisant la variable `FC`

```makefile
FC=mpif90
```
Ensuite, compiler les sources en executant la commande `make` en étant positionée au sein de la racine du répertoire de projet. 
```console
$> make
```

## Fichier de configuration

Au sein du fichier [param.txt] (https://github.com/labarraud/chp-mpi/blob/master/param.txt), indiquer selon la norme suivante les valeur des variables globales du problème, le chemain relatif du dossier où seront exporté les fichiers des solutions calculer par chaque processuset  le cas test choisi. 

Par exemple si l'on souhaite inisialisé le programme à `epsilon=0.00001` (condition d'arrêt du grandiant conjugué), `Lx=1.0`, `Ly=1.0`, `Nx=10`, `Ny=10`, `D=1.0`, `dt=1.0`, `Nt=10`, pour le cas test 2 et que l'on souhaite que les fichiers solutions soit exporté dans le dossier `Result`, le fichier de configuration est le suivant :  

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

## Usage

Un fois que le fichier [param.txt] (https://github.com/labarraud/chp-mpi/blob/master/param.txt)

```erb
<%= your_code_goes @here do |f| %>
  <%= f.input :example %>
  <%= f.input :example %>
  <%= f.button :example %>
<% end %>
```


## Configuration

This block of text should explain how to configure your application:

`rails generate my_example_gem:install`


## Information

Screenshots of your application below:

![Screenshot 1](http://placekitten.com/400/300)

![Screenshot 2](http://placekitten.com/400/300)


### Known Issues

If you discover any bugs, feel free to create an issue on GitHub fork and
send us a pull request.

[Issues List](Github Issues List URL goes here).

## Authors

* Your Name (Your Github URL goes here)
* Additional Author's name (Their Github URL goes here)


