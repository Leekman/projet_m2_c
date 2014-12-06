#include "../lib/output.h"


void sortieTerm (int scoreMasque, int **infoEnsembleT, int nombreSequences, char *motifConsensus, double **motifConsensusPSSM, int l){


	int i,j,k;

	printf("Score du masque %d\n", scoreMasque);
	printf("Numero de sequence\tPosition du motif\n");
	for (i = 0; i < 99; i++)
	{
		printf("%d\t%d\n", infoEnsembleT[i][0], infoEnsembleT[i][1] );
	}
	//printf("Motif consensus obtenu : %s\n", motifConsensus);
	/*printf("\nPSSM de ce motif consensus : \n");
    for (j=0;j<4;j++)
    {
        for(k = 0; k < l; k++)
        {
            printf("%f ", motifConsensusPSSM[j][k]);
        }
        printf("\n");
    } */
}