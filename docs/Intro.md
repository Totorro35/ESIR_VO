# Introduction

Les techniques d'asservissement se résument à l'élaboration d'une loi de commande
sur le mouvement d'une caméra. L'objectif est d'utilisé les différentes informations visuelles
à notre disposition pour les intégrer dans une boucle de commande. Ainsi, on peut estimer les différents degrées de liberté de notre robot à partir de la différence est les informations visuelles reçues et désirées. On note __x__ les informations issu de la caméra et __x<sub>d</sub>__ les valeurs désirées. L'objectif est de calculer __v__, la vitesse de la caméra de manière à minimiser __x__ - __x<sub>d</sub>__. __v__ est un vecteur de taille 6 correspond aux 3 paramètres de translation et aux 3 paramètres de rotations. On peut calculer __v__ de la manière suivante :

<img src="http://latex.codecogs.com/gif.latex?v=-\lambda L_x^+(x-x^*)" border="0"/>

avec __L<sub>x</sub>__, la matrice d'intéraction. Elle définit la relation en une information visuelle et le déplacement du robot. Tout l'objectif de l'asservissement visuel est d'exprimer correctement __v__. Elle dépend énormement des informations que l'on souhaite utiliser en entrée. En effet, on peut prendre en compte des points 2D ou 3D, des informations de géométries ou même d'autres informations issu de la sémantique de la scène.

## Changement de repère

Une fonction importante de tout système 3D est de se positionner dans le bon repère. Il est souvent nécessaire ce changer de repère. La première fonctionnalité de ce TP est donc d'établire cette fonction de changement de repère. En Vision par Ordinateur, le changement de repère le plus fréquents est le passage des coordonnées monde au repère de la caméra.

<img src="http://latex.codecogs.com/gif.latex?~^cX=~^cT_w~^wX" border="0"/>

La matrice de transformation  <b><sup>c</sup>T<sub>w</sub></b>  correspond à l'expression du repère monde selon le repère Caméra. Cette Matrice de transformation est une matrice <b>4x4</b> et les vecteurs <b><sup>c</sup>X</b> et <b><sup>w</sup>X</b> correspondent aux vecteurs homogènes de position.  La matrice de transformation s'écrit avec <b><sup>c</sup>R<sub>w</sub></b> et <b><sup>c</sup>t<sub>w</sub></b> la rotation et la translation entre les deux repères.


Dans le cadre de ce TP, il faut donc créer une fonction qui prend un vecteur à 3 dimensions et une matrice 4x4 pour retourner un vecteur à 3 dimensions.  Il est donc nécessaire de transformer le vecteur d'entrée en version homogène puis de le normaliser pour le retourner.

## Projection

Un deuxième point important en vision par ordinateur est d'être capable de passer du monde 3D à 2D. On appelle cela la projection. La projection basique et la plus simple que nous allons utilisé est la projection perspective.

<img src="pinhole1.png" alt="drawing" width="30%"/>

De manière calculatoire, on peut définir la projection avec les équations suivantes.