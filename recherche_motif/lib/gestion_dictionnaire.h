#ifndef GESTION_DICTIONNAIRE_H
#define GESTION_DICTIONNAIRE_H
#include "includes.h"
#include "structure.h"


//gestion_dictionnaire

void ajoutAuDictionnaire(dictionnaire **p_p_dictionnaire, char *k_mer, int num_sequence, int position_occurence);

void ajoutK_mer(k_mer **p_p_parcoureurK_mer_precedent, char *k_mer_courant, int numSequence, int positionOccurence);

void initialiserDictionnaire(dictionnaire **p_p_dictionnaire, char *k_merCourant, int numSequence, int positionOccurence);

void initialiserK_mer(k_mer **p_p_k_mer, char *k_merCourant, int numSequence, int positionOccurence);

void initialiserSequence(sequence **p_p_sequence, int numSequence, int positionOccurence);

void initialiserOccurence(occurence **p_p_occurence, int positionOccurence);

void updateK_mer(k_mer **p_p_parcoureurK_mer, int numSequence, int positionOccurence);

void updateSequence(sequence **p_p_parcoureurK_sequence, int positionOccurence);

// void liberationDictionnaire(dictionnaire **p_p_dictionnaire, k_mer *pk, sequence *ps, occurence *po, k_mer *nextPk, sequence *nextPs, occurence *nextPo);

void liberationDictionnaire(dictionnaire *p_dictionnaire);

// int recupNombreOccurence(k_mer *p_k_merCandidat);
#endif