#include "includes.h"
#ifndef OUTPUT_H
#define OUTPUT_H




void enTeteSortieTerm (int l, int k, int i, int *masque);

void enTeteSortieFichier (FILE *sortie, int l, int k, int i, int *masque);

void sortieTerm (int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l);

void sortieFichier (FILE *sortie, int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l);



#endif
