#ifndef MOTIF_H
#define MOTIF_H
#include "includes.h"
#include "structure.h"
#include "gestion_dictionnaire.h"
#include "fonctions_utiles.h"

//recherche_motif

void recherche_motif (int *masque, int l, int k, double **pssm, int ***infoPssmCourante, char **tableauSequences, int nombreSequences, dictionnaire **p_p_dictionnaire, double *p_score, double *motifDeFond, char ***p_ensembleT, double ***p_motifConsensusPSSM, char **p_motifConsensus, int *nbSequenceDuMotifCandidat, int *p_scoreMasque);

double **construirePSSM(k_mer *k_merCandidat, char **tableauSequences, int nombreSequences, int k);

double calculDuScore(k_mer *k_merCandidat, char **tableauSequences, int nombreSequences, int k, double **pssm, double *motifDeFond);

double calculScoreK_mer(k_mer *k_merCandidat, char **tableauSequences, int k, double **pssm, double *motifDeFond);

void ameliorerMotif(int ***infoPssmCourante, double **pssmCourante, double *p_scoreCourant, char **tableauSequences,int nombreSequences, int nombreOccurence, int k, int l, k_mer *p_k_merCandidat, double *motifDeFond);

void affinerMotif(char ***p_ensembleT, int ***infoPssmCourante, char **tableauSequences, int nombreOccurence, int l, double **pssm, double *motifDeFond, int *nbSequenceDuMotifCandidat);

void creerMotifConsensus(double ***p_motifConsensusPSSM, char **p_motifConsensus, char **ensembleT, int l, int *nbSequenceDuMotifCandidat);

int calculScoreMasque(int k, int l, char *motifConsensus, char **ensembleT, int *nbSequencesDuMotifCandidat);
#endif
