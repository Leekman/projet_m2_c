#ifndef FONCTIONS_UTILES_H
#define FONCTIONS_UTILES_H
#include "includes.h"

char **fasta_to_2Dtable(FILE *fichierSequences, int nombreSequences, int longueurSequencesMax);

void conversionMinMaj(char **tableauSequences, int nombreSequences);

double *calculerMotifDeFond(char **tableauSequences, int nombreSequences);

void copieProfondePSSM(double ***p_pssmVide, double **pssmACopier, int dim1, int dim2);

#endif