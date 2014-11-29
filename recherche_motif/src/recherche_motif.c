#include "../lib/masque.h"
#include "../lib/motif.h"

int main(){

	srand(time(NULL));

	int i;
	//int test=0;

	FILE *fichierSequences = NULL;
	char **tableauSequences = NULL;
	int nombreSequences, longueurSequencesMax;
	double *motifDeFond;

	int *masque = NULL;	
	int l,k;
	
	dictionnaire *p_dictionnaire = NULL;
	double **pssm = NULL;
	double score;

	char **ensembleT = NULL;
	double **motifConsensusPSSM;
	char *motifConsensus;

	printf("\n");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fichierSequences = fopen("input/sequences.fasta", "r+");

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

	printf("Motif de fond\n");
    for (i=0; i<4; i++)
    {
        printf("%f\n", motifDeFond[i]);
    }
    printf("\n\n");

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

	for (i=0; i<1; i++)
	{
		l=4;
		k=2;
		score = 0;

		masque=generateurMasque(l,k);

		printf("\n\n/////////////////////////////Masque %d///////////////////////////\n\n", i+1);
		/*for (j=0;j<(l-k);j++)
			printf("masque[%d] = %d\n",j,masque[j]);*/


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		recherche_motif(masque, l, k, pssm, tableauSequences, nombreSequences, &p_dictionnaire, &score, motifDeFond, &ensembleT, &motifConsensusPSSM, &motifConsensus);		
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	free (masque);
	free (motifDeFond);
	return 0;
}

