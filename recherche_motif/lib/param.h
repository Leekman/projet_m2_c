#include "includes.h"
#ifndef PARAM_H
#define PARAM_H

/******************/
/* Affiche l'aide */
/******************/
void notice();


/*********************/
/* liste des options */
/*********************/
extern struct option long_options[];


///////////////////////////////////////////////////
/*Fonction permettant de récupérer les paramètres*/
///////////////////////////////////////////////////

void getParam (char **chemin, int* longueurMotif, int* nbFenetre, int argc, char *argv[]);
#endif