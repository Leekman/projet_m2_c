#include "fonctions.h"

char **fasta_to_2Dtable(FILE *fichierSequences, int nombreSequences, int longueurSequencesMax){

	char **tableauSequences;
	char poubelle[100];

	int i;

	if (fichierSequences != NULL)
    {
    	tableauSequences=(char**)malloc(sizeof(char*)*nombreSequences);
    	for (i=0;i<nombreSequences;i++)
    	{    		
    		tableauSequences[i]=(char*)malloc(sizeof(char)*longueurSequencesMax+1);
    		fgets(poubelle, longueurSequencesMax, fichierSequences);
    		fgets(tableauSequences[i], longueurSequencesMax, fichierSequences);
    	}
        
    }
    else
    {
        printf("Impossible d'ouvrir le fichier test.txt");
    }

    fclose(fichierSequences);

    return tableauSequences;
}

double *calculerMotifDeFond(char **tableauSequences, int nombreSequences){

    double *motifDeFond;
    int nombreDeBaseTotal;

    int i,j;

    nombreDeBaseTotal = 0;
    motifDeFond=(double*)malloc(sizeof(double)*4);

    for (i=0; i<nombreSequences; i++)
    {
        nombreDeBaseTotal += (strlen(tableauSequences[i])-1);
    }

    //printf("Nombre de base : %d\n", nombreDeBaseTotal);

    for (i=0; i<nombreSequences; i++)
    {
        for (j=0; j<(strlen(tableauSequences[i])-1); j++)
        {
            switch (tableauSequences[i][j])
            {
                case 'A' : motifDeFond[0] += 1; break;
                case 'T' : motifDeFond[1] += 1; break;
                case 'C' : motifDeFond[2] += 1; break;
                case 'G' : motifDeFond[3] += 1; break;
                default : printf("Motif de fond : format de base incorect.\n"); break;
            }
        }
    }

    for (i=0; i<4; i++)
    {
        motifDeFond[i] /= nombreDeBaseTotal;
    }

    return motifDeFond;    
}