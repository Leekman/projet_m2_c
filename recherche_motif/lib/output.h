#include "includes.h"
#ifndef OUTPUT_H
#define OUTPUT_H

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define RESET "\033[0m"

void enTeteSortieTerm (int l, int k, int i, int *masque);

void enTeteSortieFichier (FILE *sortie, int l, int k, int i, int *masque);

void sortieTerm (int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l, char **ensembleT);

void sortieFichier (FILE *sortie, int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l, char **ensembleT);

void afficheMotifConsensus(double **motifConsensusPSSM, char *motifConsensus);

#endif
