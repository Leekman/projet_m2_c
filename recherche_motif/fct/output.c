#include "../lib/output.h"
#include "../lib/fonctions_utiles.h"





void enTeteSortieTerm (int l, int k, int i, int *masque){

	int j;
	printf("\n\n/////////////////////////////Masque %d/////////////////////////////\n", i+1);
	printf("\nInformations sur le masque :\n");
	for (j=0;j<(l-k);j++)
		{
			printf("\nMasque[%d] : case %d masquée\n",j+1, masque[j]+1); //on ajoute +1 pour une lecture plus compréhensible par l'utilisateur
		}
}

void enTeteSortieFichier (FILE *sortie, int l, int k, int i, int *masque){
	
	int j;
	fprintf(sortie, "\n\n/////////////////////////////Masque %d/////////////////////////////\n", i+1);
	fprintf(sortie, "\nInformations sur le masque :\n");
	for (j=0;j<(l-k);j++)
		{
			fprintf(sortie, "\nMasque[%d] : case %d masquée\n",j+1, masque[j]+1);
		}	
}



void sortieTerm (int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l, char **ensembleT){


	int i,j,k;
	char nucleotide[] = "ATCG";
	printf("\nScore du masque %d\n", scoreMasque);
	printf("\nNumero de sequence | Position du motif | Motif \n");
	for (i = 0; i < nbSequenceDuMotifConsensus; i++)
	{
		printf("%d | %d | %s\n", infoEnsembleT[i][0] +1, infoEnsembleT[i][1], ensembleT[i]);
	}
	printf("\nMotif consensus obtenu : %s\n", motifConsensus);
	printf("PSSM de ce motif consensus : \n\n");
    for (j=0;j<4;j++)
    {
    	printf("%c\t", nucleotide[j]);
    	for(k = 0; k < l; k++)
        {
            printf("%f ", motifConsensusPSSM[j][k]);
        }
        printf("\n");
    } 
}

void sortieFichier (FILE* sortie, int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l, char **ensembleT){

	int i,j,k;
	char nucleotide[] = "ATCG";

	if (sortie != NULL)
	{
		fprintf(sortie, "\nScore du masque %d\n", scoreMasque);
		fprintf(sortie, "\nNumero de sequence | Position du motif | Motif\n");
		for (i = 0; i < nbSequenceDuMotifConsensus; i++)
		{
			fprintf(sortie, "%d | %d | %s\n", infoEnsembleT[i][0] +1, infoEnsembleT[i][1], ensembleT[i]);
		}
		fprintf(sortie, "\nMotif consensus obtenu : %s\n", motifConsensus);
		fprintf(sortie, "PSSM de ce motif consensus : \n\n");
	    for (j=0;j<4;j++)
	    {
	    	fprintf(sortie, "%c\t", nucleotide[j]);
	    	for(k = 0; k < l; k++)
	        {
	            fprintf(sortie, "%f ", motifConsensusPSSM[j][k]);
	        }
	        fprintf(sortie, "\n");
   		} 

	}

	for (i = 0; i < nbSequenceDuMotifConsensus; i++)
	{	
		free(infoEnsembleT[i]);
		free(ensembleT[i]);
	}
	free(infoEnsembleT);
	free(ensembleT);
	liberationMemoirePSSM(motifConsensusPSSM);
	fclose(sortie);
}