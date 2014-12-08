#include "../lib/fonctions_utiles.h"


int recupNbSeq(FILE *fichierSequences){

    int c;
    int nbSeq = 0;
 

    if (fichierSequences != NULL)
    {
        while((c=fgetc(fichierSequences)) != EOF)
        {
            if(c=='\n')
                nbSeq++;
        }
            nbSeq++;
            rewind (fichierSequences);
    }
    return nbSeq/2;
}



char **fasta_to_2Dtable(FILE *fichierSequences, int nombreSequences){

	char **tableauSequences;
	//char *poubelle;

	int i;
    int c;
    int taillePoubelle = 0;
    int tailleSequence = 0;
    int positionCurseurLectureSequence = 0;
	if (fichierSequences != NULL)
    {
    	tableauSequences=(char**)malloc(sizeof(char*)*nombreSequences);
    	for (i = 0 ; i < nombreSequences; i++)
    	{   
            taillePoubelle = 0;
            tailleSequence = 0;
            positionCurseurLectureSequence = 0;
            while ((c=fgetc(fichierSequences)) != '\n')
            {
                taillePoubelle++;
            }
            //printf("taille poubelle :%d\n", taillePoubelle);
            //printf("Position du curseur après seq1> : %d\n", ftell (fichierSequences));
            do 
            {
                c=fgetc(fichierSequences); 
                tailleSequence++;
                positionCurseurLectureSequence++;
                if (c == EOF)
                {
                    positionCurseurLectureSequence--;
                }
            } while (c != '\n' && c != EOF);
            //printf("taille sequence :%d\n", tailleSequence);
            tableauSequences[i] = (char*) calloc (tailleSequence, sizeof(char));
            fseek (fichierSequences, -positionCurseurLectureSequence, SEEK_CUR);
            //printf("%d\n", tailleSequence);
            //printf("Position du curseur  après fseek: %d\n", ftell (fichierSequences));
    		fgets(tableauSequences[i], tailleSequence, fichierSequences);
            fseek (fichierSequences, 1, SEEK_CUR);
            //printf("sequence : %s\n", tableauSequences[i]);
    	}
        
    }
    else
    {
        printf("Impossible d'ouvrir le fichier");
        exit(1);
    }

    fclose(fichierSequences);

    return tableauSequences;
}


void conversionMinMaj(char **tableauSequences, int nombreSequences){

	int i,j;

	for (i = 0; i < nombreSequences; i++)
		{
			for (j = 0; j < strlen(tableauSequences[i]); j++)
			{
					tableauSequences[i][j] = toupper(tableauSequences[i][j]);
			}
		}
}



double *calculerMotifDeFond(char **tableauSequences, int nombreSequences){

    double *motifDeFond;
    int nombreDeBaseTotal;

    int i,j;

    nombreDeBaseTotal = 0;
    motifDeFond=(double*)calloc(4,sizeof(double));

    for (i=0; i<nombreSequences; i++)
    {
        nombreDeBaseTotal += (strlen(tableauSequences[i]));
    }

    //printf("Nombre de base : %d\n", nombreDeBaseTotal);

    for (i=0; i<nombreSequences; i++)
    {
        for (j=0; j<(strlen(tableauSequences[i])); j++)
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

void copieProfondePSSM(double ***p_pssmVide, double **pssmACopier, int dim1, int dim2){

    int i, j;


    if (*p_pssmVide != NULL)
    {
        free(*p_pssmVide);
        *p_pssmVide = NULL;
    }
    (*p_pssmVide)=(double**)malloc(sizeof(double*)*dim1);
    for (i=0; i<dim1; i++)
    {
        (*p_pssmVide)[i]=(double*)malloc(sizeof(double)*dim2);
        for (j=0; j<dim2; j++)
        {
            (*p_pssmVide)[i][j]=pssmACopier[i][j];
        }

    }
}