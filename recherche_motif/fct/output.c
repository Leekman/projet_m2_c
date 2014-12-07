#include "../lib/output.h"


void sortieTerm (int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l){


	int i,j,k;
	char nucleotide[] = "ATCG";
	printf("Score du masque %d\n", scoreMasque);
	printf("Numero de sequence | Position du motif\n");
	for (i = 0; i < nbSequenceDuMotifConsensus; i++)
	{
		printf("%d | %d\n", infoEnsembleT[i][0] +1, infoEnsembleT[i][1] );
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

void sortieFichier (FILE* sortie, int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l){

	int i,j,k;
	char nucleotide[] = "ATCG";

	if (sortie != NULL)
	{
		fprintf(sortie, "Score du masque %d\n", scoreMasque);
		fprintf(sortie, "Numero de sequence | Position du motif\n");
		for (i = 0; i < nbSequenceDuMotifConsensus; i++)
		{
			fprintf(sortie, "%d | %d\n", infoEnsembleT[i][0] +1, infoEnsembleT[i][1] );
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
	fclose(sortie);
}