#include "../lib/output.h"
#include "../lib/fonctions_utiles.h"





void enTeteSortieTerm (int l, int k, int i, int *masque){

	int j;
	printf("\n\n/////////////////////////////Masque %d/////////////////////////////\n", i+1);
	printf("\nInformations sur le masque :\n");
	printf("Fenetre(s) ");
	for (j=0;j<(l-k);j++)
	{
		printf("%d ", masque[j]+1); //on ajoute +1 pour une lecture plus compréhensible par l'utilisateur
	}
	printf(" masquee(s)\n");
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


	int i,j;
	char nucleotide[] = "ATCG";
	printf("\nScore du masque %d\n", scoreMasque);
	/*printf("\nNumero de sequence | Position du motif | Motif \n");
	for (i = 0; i < nbSequenceDuMotifConsensus; i++)
	{
		printf("%d | %d | %s\n", infoEnsembleT[i][0] +1, infoEnsembleT[i][1], ensembleT[i]);
	}*/
	afficheMotifConsensus(motifConsensusPSSM, motifConsensus);
	//printf("\nMotif consensus obtenu : %s\n", motifConsensus);
	printf("PSSM de ce motif consensus : \n\n");
    for (i=0;i<4;i++)
    {
    	printf("%c\t", nucleotide[i]);
    	for(j = 0; j < l; j++)
        {
            printf("%f ", motifConsensusPSSM[i][j]);
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

void afficheMotifConsensus(double **motifConsensusPSSM, char *motifConsensus){
	int i;
	int tailleMotif;
	tailleMotif = strlen(motifConsensus);
	printf("\nMotif consensus obtenu : ");
	for (i = 0; i < tailleMotif; i++)
	{
		switch(motifConsensus[i])
		{
			case 'A' : 
			{
				if (motifConsensusPSSM[0][i] < 0.50) 
				{
					printf(KRED "A" RESET);
				}
				else 
				{
					if (motifConsensusPSSM[0][i] < 0.75)
					{
						printf(KYEL "A" RESET);
					}
					else
					{
						printf(KGRN "A" RESET);
					}
				}
				break;
			}
			case 'T' : 
			{
				if (motifConsensusPSSM[1][i] < 0.50) 
				{
					printf(KRED "T" RESET);
				}
				else 
				{
					if (motifConsensusPSSM[1][i] < 0.75)
					{
						printf(KYEL "T" RESET);
					}
					else
					{
						printf(KGRN "T" RESET);
					}
				}				
				break;
			}
			case 'C' : 
			{
				if (motifConsensusPSSM[2][i] < 0.50) 
				{
					printf(KRED "C" RESET);
				}
				else 
				{
					if (motifConsensusPSSM[2][i] < 0.75)
					{
						printf(KYEL "C" RESET);
					}
					else
					{
						printf(KGRN "C" RESET);
					}
				}				
				break;
			}
			case 'G' :
			{
				if (motifConsensusPSSM[3][i] < 0.50) 
				{
					printf(KRED "G" RESET);
				}
				else 
				{
					if (motifConsensusPSSM[3][i] < 0.75)
					{
						printf(KYEL "G" RESET);
					}
					else
					{
						printf(KGRN "G" RESET);
					}
				}			
				break;
			}
			default :
			{
				printf("X");
			}
		}
	}
	printf("\nCode couleur : ");
	printf(KGRN "\tTres peu variable (> 0,75), " RESET);
	printf(KYEL "peu variable (> 0,5), " RESET);
	printf(KRED "variable (< 0,5)." RESET);
	printf("\n\n");
}