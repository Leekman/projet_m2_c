#include "../lib/masque.h"
#include "../lib/motif.h"
#include "../lib/param.h"
#include "../lib/output.h"
int main(int argc, char *argv[]){

	srand(time(NULL));

	int i;

	char *chemin = NULL;
	FILE *fichierSequences = NULL;

	char **tableauSequences = NULL;
	int nombreSequences, longueurSequencesMax;
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



	printf("\n");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	chemin = getParam(p_l, p_k, argc, argv);
	printf("%s\n", chemin);
	fichierSequences = fopen(chemin, "r+");

	nombreSequences = 100;
	longueurSequencesMax = 150;

	tableauSequences=fasta_to_2Dtable(fichierSequences, nombreSequences, longueurSequencesMax);
	
	/*for (i = 0; tableauSequences[i]; i++);
	test=i;
	printf("%d\n", test);*/

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

		score = 0;

		masque=generateurMasque(l,k);

		printf("\n\n/////////////////////////////Masque %d/////////////////////////////\n\n", i+1);
		/*for (j=0;j<(l-k);j++)
			printf("masque[%d] = %d\n",j,masque[j]);*/


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		recherche_motif(masque, l, k, pssm, &infoEnsembleT, tableauSequences, nombreSequences, &p_dictionnaire, &score, motifDeFond, &ensembleT, &motifConsensusPSSM, &motifConsensus, &scoreMasque);		

	///////////////////////////
	/*SORTIE DANS LE TERMINAL*/
	///////////////////////////

		//sortieTerm(scoreMasque, infoEnsembleT, nombreSequences, motifConsensus, motifConsensusPSSM, l);

	/*SORTIE FICHIER*/
		//sortieFichier(scoreMasque, infoEnsembleT, l, motifConsensus, motifConsensusPSSM);

	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	free (masque);
	free (motifDeFond);
	return 0;
}

