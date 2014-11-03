#ifndef PARAM_H
#define PARAM_H
#include "includes.h"


/******************/
/* Affiche l'aide */
/******************/
void notice();


/*********************/
/* liste des options */
/*********************/
extern struct option long_options[];

/*************************************************************/
/* remplissage des paramètres avec les valeurs des arguments */
/*************************************************************/

void getParam (int* nbErreur, int* nbSeq, int* tailleSeq, char* motif, int argc, char *argv[]);

 #endif