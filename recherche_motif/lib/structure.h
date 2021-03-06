#ifndef STRUCTURE_H
#define STRUCTURE_H
#include "includes.h"

//structures
typedef struct occurence{
	int position;
	struct occurence *nextOccurence;
} occurence;

typedef struct sequence{
	int numSequence;
	struct occurence *firstOccurence;
	struct sequence *nextSequence;
} sequence;

typedef struct k_mer{
	char *k_mer;
	struct sequence *firstSequence;
	struct k_mer *nextK_mer;
} k_mer;

typedef struct dictionnaire{
	struct k_mer *firstK_mer;
} dictionnaire;

typedef struct resultat{
	int l, k, i;
	int *masque;
	int scoreMasque;
	int **infoEnsembleT;
	int nbSequenceDuMotifConsensus;
	char *motifConsensus;
	double **motifConsensusPSSM;
	char **ensembleT;
	struct resultat *nextRes;
} resultat;

#endif