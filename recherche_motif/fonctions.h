#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//fonctions_utiles

char **fasta_to_2Dtable(FILE *fichierSequences, int nombreSequences, int longueurSequencesMax);

double *calculerMotifDeFond(char **tableauSequences, int nombreSequences);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//masque

int compare(void const *a, void const *b);

bool estDoublon(int n, int *doublon, int count);

int *generateurMasque(int l, int k);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//recherche_motif

void recherche_motif (int *masque, int l, int k, double **pssm, char **tableauSequences, int nombreSequences, dictionnaire **p_p_dictionnaire, double *p_score, double *motifDeFond); 

double **construirePSSM(k_mer *k_merCandidat, char **tableauSequences, int nombreSequences, int k);

double calculDuScore(k_mer *k_merCandidat, char **tableauSequences, int nombreSequences, int k, double **pssm, double *motifDeFond);

double calculProbaPSSM(k_mer *k_merCandidat, char **tableauSequences, int k, double **pssm, double *motifDeFond);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//gestion_dictionnaire

void ajoutAuDictionnaire(dictionnaire **p_p_dictionnaire, char *k_mer, int num_sequence, int position_occurence);

void ajoutK_mer(k_mer **p_p_parcoureurK_mer_precedent, char *k_mer_courant, int numSequence, int positionOccurence);

void initialiserDictionnaire(dictionnaire **p_p_dictionnaire, char *k_merCourant, int numSequence, int positionOccurence);

void initialiserK_mer(k_mer **p_p_k_mer, char *k_merCourant, int numSequence, int positionOccurence);

void initialiserSequence(sequence **p_p_sequence, int numSequence, int positionOccurence);

void initialiserOccurence(occurence **p_p_occurence, int positionOccurence);

void updateK_mer(k_mer **p_p_parcoureurK_mer, int numSequence, int positionOccurence);

void updateSequence(sequence **p_p_parcoureurK_sequence, int positionOccurence);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

