#include "fonctions.h"

void ajoutAuDictionnaire(dictionnaire **p_p_dictionnaire, char *k_merCourant, int numSequence, int positionOccurence){
	
	k_mer *p_parcoureurK_mer = NULL;
	k_mer *p_parcoureurK_merPrecedent = NULL;

	//printf("%s\n", k_merCourant);

	if (*p_p_dictionnaire==NULL)
	{
		initialiserDictionnaire(p_p_dictionnaire, k_merCourant, numSequence, positionOccurence);	
	}

	p_parcoureurK_mer=(*p_p_dictionnaire)->firstK_mer;
	p_parcoureurK_merPrecedent=(*p_p_dictionnaire)->firstK_mer;


	while (p_parcoureurK_mer != NULL && strcmp(p_parcoureurK_mer->k_mer,k_merCourant)!=0)
	{
		p_parcoureurK_merPrecedent=p_parcoureurK_mer;
		p_parcoureurK_mer=p_parcoureurK_mer->nextK_mer;
	}

	if (p_parcoureurK_mer == NULL)
	{
		//On est allé jusq'au bout de la liste de k_mer, et on a pas trouvé le bon k_mer
		ajoutK_mer(&p_parcoureurK_merPrecedent, k_merCourant, numSequence, positionOccurence);
		return ;
	}
	if (p_parcoureurK_mer != NULL)
	{
		//On est pas allé jusqu'au bout de la liste, donc le k_mer existe deja
		updateK_mer(&p_parcoureurK_mer, numSequence, positionOccurence);
		return ;
	}
	//ajoutQueueDictionnaire(p_p_dictionnaire, k_mer, numSequence, positionOccurence);
}

void ajoutK_mer(k_mer **p_p_parcoureurK_merPrecedent, char *k_merCourant, int numSequence, int positionOccurence){

	k_mer *p_k_mer = NULL;

	initialiserK_mer(&p_k_mer, k_merCourant, numSequence, positionOccurence);
	(*p_p_parcoureurK_merPrecedent)->nextK_mer=p_k_mer;
	
}

void initialiserDictionnaire(dictionnaire **p_p_dictionnaire, char *k_merCourant, int numSequence, int positionOccurence){

	k_mer *p_k_mer = NULL;

	*p_p_dictionnaire=(dictionnaire*)malloc(sizeof(dictionnaire));
	initialiserK_mer(&p_k_mer, k_merCourant, numSequence, positionOccurence);
	(*p_p_dictionnaire)->firstK_mer = p_k_mer;
}

void initialiserK_mer(k_mer **p_p_k_mer, char *k_merCourant, int numSequence, int positionOccurence){

	sequence *p_sequence = NULL;

	*p_p_k_mer=(k_mer*)malloc(sizeof(k_mer));
	(*p_p_k_mer)->k_mer=(char*)malloc(sizeof(char)*strlen(k_merCourant));
	strcpy((*p_p_k_mer)->k_mer,k_merCourant);
	initialiserSequence(&p_sequence, numSequence, positionOccurence);
	(*p_p_k_mer)->firstSequence = p_sequence;
	(*p_p_k_mer)->nextK_mer = NULL;
}

void initialiserSequence(sequence **p_p_sequence, int numSequence, int positionOccurence){

	occurence *p_occurence = NULL;

	*p_p_sequence=(sequence*)malloc(sizeof(sequence));
	(*p_p_sequence)->numSequence = numSequence;
	initialiserOccurence(&p_occurence, positionOccurence);
	(*p_p_sequence)->firstOccurence = p_occurence;
	(*p_p_sequence)->nextSequence = NULL;
}

void initialiserOccurence(occurence **p_p_occurence, int positionOccurence){

	*p_p_occurence=(occurence*)malloc(sizeof(occurence));
	(*p_p_occurence)->position = positionOccurence;
	(*p_p_occurence)->nextOccurence = NULL;
}

void updateK_mer(k_mer **p_p_parcoureurK_mer, int numSequence, int positionOccurence){

	sequence *p_parcoureurSequence = NULL;
	sequence *p_parcoureurSequencePrecedente = NULL;
	sequence *p_sequence = NULL;

	p_parcoureurSequence = (*p_p_parcoureurK_mer)->firstSequence;
	p_parcoureurSequencePrecedente = (*p_p_parcoureurK_mer)->firstSequence;

	while (p_parcoureurSequence != NULL && p_parcoureurSequence->numSequence != numSequence)
	{
		p_parcoureurSequencePrecedente = p_parcoureurSequence;
		p_parcoureurSequence = p_parcoureurSequence->nextSequence;
	}

	if (p_parcoureurSequence == NULL)
	{
		initialiserSequence(&p_sequence, numSequence, positionOccurence);
		p_parcoureurSequencePrecedente->nextSequence = p_sequence;
	}
	else
	{
		updateSequence(&p_parcoureurSequence, positionOccurence);
	}
}

void updateSequence(sequence **p_p_parcoureurK_sequence, int positionOccurence){

	occurence *p_parcoureurOccurence = NULL;
	occurence *p_occurence = NULL;

	p_parcoureurOccurence = (*p_p_parcoureurK_sequence)->firstOccurence;

	while (p_parcoureurOccurence->nextOccurence != NULL)
	{
		p_parcoureurOccurence = p_parcoureurOccurence->nextOccurence;
	}

	initialiserOccurence(&p_occurence, positionOccurence);
	p_parcoureurOccurence->nextOccurence = p_occurence;
}