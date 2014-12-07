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
	char *chemin = NULL;
	FILE *fichierSequences = NULL;
	FILE *sortie = NULL;
	char **tableauSequences = NULL;
	int nombreSequences;
	double *motifDeFond;
	int scoreMasque;
	int *masque = NULL;	
	int l,k;
	int *p_l = NULL;
	int *p_k = NULL;
	p_l = &l;
	p_k = &k;
	dictionnaire *p_dictionnaire = NULL;
	double **pssm = NULL;
	double score;
	int **infoEnsembleT = NULL;
	char **ensembleT = NULL;
	double **motifConsensusPSSM;
	char *motifConsensus;
	int nbSequenceDuMotifConsensus = 0;


	///////////////////////////////
	/*RECUPERATION DES PARAMETRES*/
	///////////////////////////////
	getParam(&chemin, p_l, p_k, argc, argv);
	
	/////////////////////////////////////////////////////////////
	/*RECUPERATION DU NOMBRE DE SEQUENCES DANS LE FICHIER FASTA*/
	/////////////////////////////////////////////////////////////
	fichierSequences = fopen(chemin, "r"); // ouverture du fichier en mode lecture uniquement

	nombreSequences = recupNbSeq (fichierSequences);
	rewind (fichierSequences); //on remet le curseur de lecture du fichier au début du fichier
	//longueurSequencesMax = 300;
	tableauSequences=fasta_to_2Dtable(fichierSequences, nombreSequences);

	/////////////////////////////////////////
	/*MODIFICATION MINUSCULES EN MAJUSCULES*/
	/////////////////////////////////////////

	conversionMinMaj(tableauSequences, nombreSequences);

	motifDeFond=calculerMotifDeFond(tableauSequences, nombreSequences);

	// printf("Motif de fond\n");
 //    for (i=0; i<4; i++)
 //    {
 //        printf("%f\n", motifDeFond[i]);
 //    }
 //    printf("\n\n");

	/*for (i=0; i<nombreSequences; i++)
		printf("%s", tableauSequences[i]);*/
    /*i=0;
    while (tableauSequences[i]){
    	test++;
    	i++;
    }
    printf("%d\n", test);
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i=0; i<2; i++)
	{
		sortie = fopen("output/resultats.info", "a");
		score = 0;

		masque=generateurMasque(l,k);

		printf("\n\n/////////////////////////////Masque %d/////////////////////////////\n\n", i+1);
		fprintf(sortie,"\n\n/////////////////////////////Masque %d/////////////////////////////\n\n", i+1);
		/*for (j=0;j<(l-k);j++)
			printf("masque[%d] = %d\n",j,masque[j]);*/


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		recherche_motif(masque, l, k, pssm, &infoEnsembleT, tableauSequences, nombreSequences, &p_dictionnaire, &score, motifDeFond, &ensembleT, &motifConsensusPSSM, &motifConsensus, &nbSequenceDuMotifConsensus, &scoreMasque);		

	///////////////////////////
	/*SORTIE DANS LE TERMINAL*/
	///////////////////////////
		sortieTerm(scoreMasque, infoEnsembleT, nbSequenceDuMotifConsensus, motifConsensus, motifConsensusPSSM, l);

	//////////////////
	/*SORTIE FICHIER*/
	//////////////////
		sortieFichier(sortie, scoreMasque, infoEnsembleT, nbSequenceDuMotifConsensus, motifConsensus, motifConsensusPSSM, l);

	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	free (masque);
	free (motifDeFond);
	return 0;
}

