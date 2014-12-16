#ifndef PARAM_H
#define PARAM_H
#include "includes.h"

/////////////////////
/*AFFICHE LA NOTICE*/
/////////////////////

void notice();

/////////////////////
/*LISTE DES OPTIONS*/
/////////////////////
extern struct option long_options[];

/////////////////////////////////////////////////////////////
/*REMPLISSAGE DES PARAMETRES AVEC LES VALEURS DES ARGUMENTS*/
/////////////////////////////////////////////////////////////

void getParam (int *nbErreur, int *nbSeq, int *tailleSeq, char **motif, char **cheminFasta, char** cheminInfo, int *variation, int argc, char *argv[]);

#endif