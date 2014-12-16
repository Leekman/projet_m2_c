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
	printf(" ouverte(s)\n");
}

void enTeteSortieFichier (FILE *sortie, int l, int k, int i, int *masque){
	
	int j;
	fprintf(sortie, "\n\n/////////////////////////////Masque %d/////////////////////////////\n", i+1);
	fprintf(sortie, "\nInformations sur le masque :\n");
	for (j=0;j<(l-k);j++)
		{
			fprintf(sortie, "\nMasque[%d] : case %d ouverte\n",j+1, masque[j]+1);
		}	
}



void sortieTerm (int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l, char **ensembleT, int nombreDeSequences){


	//int i, j;
	//char nucleotide[] = "ATCG";
	printf("\nScore du masque : %d\n", scoreMasque);
	/*printf("\nNumero de sequence | Position du motif | Motif \n");
	for (i = 0; i < nbSequenceDuMotifConsensus; i++)
	{
		printf("%d | %d | %s\n", infoEnsembleT[i][0] +1, infoEnsembleT[i][1], ensembleT[i]);
	}*/
	afficheMotifConsensus(motifConsensusPSSM, motifConsensus);
	printf("Quorum du motif : %d%%\n", (int)(((double)nbSequenceDuMotifConsensus/(double)nombreDeSequences)*100));
	//printf("\nMotif consensus obtenu : %s\n", motifConsensus);
	/*printf("PSSM de ce motif consensus : \n\n");
    for (i=0;i<4;i++)
    {
    	printf("%c\t", nucleotide[i]);
    	for(j = 0; j < l; j++)
        {
            printf("%f ", motifConsensusPSSM[i][j]);
        }
        printf("\n");
    } */
}

void sortieFichier (FILE* sortie, int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, int l, char **ensembleT, int nombreDeSequences){

	int i,j,k;
	char nucleotide[] = "ATCG";

	if (sortie != NULL)
	{
		fprintf(sortie, "\nScre du masque : %d\n", scoreMasque);
		fprintf(sortie, "\nNumero de sequence | Position du motif | Motif\n");
		for (i = 0; i < nbSequenceDuMotifConsensus; i++)
		{
			fprintf(sortie, "%d | %d | %s\n", infoEnsembleT[i][0] +1, infoEnsembleT[i][1], ensembleT[i]);
		}
		fprintf(sortie, "\nMotif consensus obtenu : %s\n", motifConsensus);
		fprintf(sortie, "Quorum du motif : %f\n", (double)nbSequenceDuMotifConsensus/(double)nombreDeSequences);
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

	/*for (i = 0; i < nbSequenceDuMotifConsensus; i++)
	{	
		free(infoEnsembleT[i]);
		free(ensembleT[i]);
	}
	free(infoEnsembleT);
	free(ensembleT);
	liberationMemoirePSSM(motifConsensusPSSM);*/
	//fclose(sortie);
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
	printf("\n\n");
}

void ajouterResultat(resultat **p_listeResultats,resultat* p_resultat, int l, int k, int i, int *masque, int scoreMasque, int **infoEnsembleT, int nbSequenceDuMotifConsensus, char *motifConsensus, double **motifConsensusPSSM, char **ensembleT){

	int m;
	resultat *parcourResultat = NULL;
	resultat *parcourResultatPrec = NULL;

	p_resultat->l = l;
	p_resultat->k = k;
	p_resultat->i = i;
	p_resultat->masque = (int*)calloc(sizeof(int),l-k);
	for (m = 0; m < (l-k); m++)
	{
		//printf("%d : %d\n", m, masque[m]);
		p_resultat->masque[m] = masque[m];
	}
	p_resultat->scoreMasque = scoreMasque;
	p_resultat->infoEnsembleT = NULL;
	copieProfondeInt2D(&(p_resultat->infoEnsembleT), infoEnsembleT, nbSequenceDuMotifConsensus, 2);
	p_resultat->nbSequenceDuMotifConsensus = nbSequenceDuMotifConsensus;
	p_resultat->motifConsensus = (char*)calloc(sizeof(char),(l+1));
	strcpy(p_resultat->motifConsensus, motifConsensus);
	p_resultat->motifConsensusPSSM = NULL;
	copieProfondePSSM(&(p_resultat->motifConsensusPSSM), motifConsensusPSSM, 4, l);
	p_resultat->ensembleT = NULL;
	copieProfondeTabString(&(p_resultat->ensembleT), ensembleT, nbSequenceDuMotifConsensus, l+1);
	p_resultat->nextRes = NULL;
	//free(ensembleT);
	parcourResultat = (*p_listeResultats);

	if (parcourResultat == NULL)
	{
		(*p_listeResultats) = p_resultat;
		p_resultat = NULL;
		//printf("First try. \n");
		return ;
	}

	while(parcourResultat != NULL && parcourResultat->scoreMasque <= scoreMasque)
	{
		parcourResultatPrec = parcourResultat;
		parcourResultat = parcourResultat->nextRes;
	}

	if (parcourResultat == NULL && parcourResultatPrec != NULL)
	{
		//printf("Ajout en queue.\n");		
		parcourResultatPrec->nextRes = p_resultat;
		p_resultat = NULL;

	}
	else
	{
		if (parcourResultatPrec == NULL)
		{
			//printf("Ajout en tete.\n");
			p_resultat->nextRes = parcourResultat;
			(*p_listeResultats) = p_resultat;
			p_resultat = NULL;

		}
		else
		{
			//printf("Ajout en milieu.\n");
			p_resultat->nextRes = parcourResultat;
			parcourResultatPrec->nextRes = p_resultat;
			p_resultat = NULL;
			
		}
	}	
	free(p_resultat);
}

void afficherSortie(FILE *sortie, resultat *listeResultats, int nombreDeSequences){
	if (listeResultats != NULL)
	{
		//ECRIT
		
		enTeteSortieFichier(sortie, listeResultats->l, listeResultats->k , listeResultats->i, listeResultats->masque);
		sortieFichier(sortie, listeResultats->scoreMasque, listeResultats->infoEnsembleT, listeResultats->nbSequenceDuMotifConsensus, listeResultats->motifConsensus, listeResultats->motifConsensusPSSM, listeResultats->l, listeResultats->ensembleT, nombreDeSequences);
		
		//APPEL RECURSIF
		afficherSortie(sortie, listeResultats->nextRes, nombreDeSequences);

		//AFFICHE
		enTeteSortieTerm(listeResultats->l, listeResultats->k, listeResultats->i, listeResultats->masque);
		sortieTerm(listeResultats->scoreMasque, listeResultats->infoEnsembleT, listeResultats->nbSequenceDuMotifConsensus, listeResultats->motifConsensus, listeResultats->motifConsensusPSSM, listeResultats->l, listeResultats->ensembleT, nombreDeSequences);

		//FREE
		//printf("%s\n", listeResultats->motifConsensus);
		libererMemoireResultat(listeResultats);
		free(listeResultats);
	}
}

void libererMemoireResultat(resultat *p_resultat){

	int i;

	for (i = 0; i < p_resultat->nbSequenceDuMotifConsensus; i++)
	{
		free((p_resultat->infoEnsembleT)[i]);
	}
	free(p_resultat->infoEnsembleT);

	free(p_resultat->motifConsensus);

	free(p_resultat->masque);

	liberationMemoirePSSM(p_resultat->motifConsensusPSSM);

	for (i = 0; i < p_resultat->nbSequenceDuMotifConsensus; i++)
	{
		free((p_resultat->ensembleT)[i]);
	}
	free(p_resultat->ensembleT);

}