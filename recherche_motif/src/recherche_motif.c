#include "../lib/masque.h"
#include "../lib/motif.h"
#include "../lib/param.h"
#include "../lib/output.h"
int main(int argc, char *argv[]){


	////////////////////////////////////////////////
    /*Initialisation du rand au time du processeur*/
    ////////////////////////////////////////////////
	srand(time(NULL));


    ///////////////////////////////////////////////
    /*DECLARATION ET INITIALISATION DES VARIABLES*/
    ///////////////////////////////////////////////
	int i;
	int scoreMasque;
	int nombreSequences, nbSequenceDuMotifConsensus;
	int l,k;
	int nbErreurMax = 0;
	int nbIterations;
	double score;
	int *masque = NULL;
	double *motifDeFond = NULL;
	char *cheminEntree = NULL;
	char *cheminSortie = NULL;
	char *motifConsensus = NULL;
	dictionnaire *p_dictionnaire = NULL;
	FILE *fichierSequences = NULL;
	FILE *sortie = NULL;
	int **infoEnsembleT = NULL;
	double **pssm = NULL;
	double **motifConsensusPSSM = NULL;
	char **tableauSequences = NULL;
	char **ensembleT = NULL;



	///////////////////////////////
	/*RECUPERATION DES PARAMETRES*/
	///////////////////////////////
	getParam(&cheminEntree, &cheminSortie, &nbIterations, &nbErreurMax, &l, &k, argc, argv);
	
	/////////////////////////////////////////////////////////////
	/*RECUPERATION DU NOMBRE DE SEQUENCES DANS LE FICHIER FASTA*/
	/////////////////////////////////////////////////////////////
 	fichierSequences = fopen(cheminEntree, "r"); // ouverture du fichier en mode lecture uniquement

 	nombreSequences = recupNbSeq (fichierSequences);
	
	tableauSequences=fasta_to_2Dtable(fichierSequences, nombreSequences);

	/////////////////////////////////////////
	/*MODIFICATION MINUSCULES EN MAJUSCULES*/
	/////////////////////////////////////////

	conversionMinMaj(tableauSequences, nombreSequences);
	
 	///////////////////////////
 	/*CALCUL DU MOTIF DE FOND*/
 	///////////////////////////

	motifDeFond=calculerMotifDeFond(tableauSequences, nombreSequences);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	for (i = 0 ; i < nbIterations; i++)
	{
		sortie = fopen(cheminSortie, "a+");
		score = 0;

		masque=generateurMasque(l,k);

		enTeteSortieTerm(l, k, i, masque);
		enTeteSortieFichier(sortie, l, k , i, masque);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		recherche_motif(masque, l, k, pssm, &infoEnsembleT, tableauSequences, nombreSequences, &p_dictionnaire, &score, motifDeFond, &ensembleT, &motifConsensusPSSM, &motifConsensus, &nbSequenceDuMotifConsensus, &scoreMasque, nbErreurMax);		

		/////////////////////////////////////////
		/*SORTIE DANS LE TERMINAL ET LE FICHIER*/
		/////////////////////////////////////////

		if (motifConsensus != NULL)
		{
			sortieTerm(scoreMasque, infoEnsembleT, nbSequenceDuMotifConsensus, motifConsensus, motifConsensusPSSM, l, ensembleT);
			sortieFichier(sortie, scoreMasque, infoEnsembleT, nbSequenceDuMotifConsensus, motifConsensus, motifConsensusPSSM, l, ensembleT);
		}
		else
		{
			fprintf(sortie, "Aucun motif commun trouvé avec ce masque\n");
			fclose(sortie);	
		}

		free(motifConsensus);
		free (masque);
	}

	for (i = 0; i < nombreSequences; i++)
	{
		free (tableauSequences[i]);
	}
	free(tableauSequences);
	free(cheminEntree); 
	free(cheminSortie);
	free(motifDeFond);

	return 0;
}

