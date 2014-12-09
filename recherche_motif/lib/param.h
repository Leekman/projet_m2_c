#include "includes.h"
#ifndef PARAM_H
#define PARAM_H

//////////////////
/*AFFICHE L'AIDE*/
//////////////////
void notice();


/////////////////////
/*LISTE DES OPTIONS*/
/////////////////////
extern struct option long_options[];


///////////////////////////////////////////////////
/*FONCTION PERMETTANT DE RECUPERER LES PARAMETRES*/
///////////////////////////////////////////////////

void getParam (char **cheminEntree, char **output,int *nbIterations, int *nbErreurMax, int* longueurMotif, int* nbFenetre, int argc, char *argv[]);

#endif