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
            tableauSequences[i] = (char*) calloc (tailleSequence, sizeof(char));
            fseek (fichierSequences, -positionCurseurLectureSequence, SEEK_CUR);
    		fgets(tableauSequences[i], tailleSequence, fichierSequences);
            fseek (fichierSequences, 1, SEEK_CUR);
            
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
    motifDeFond=(double*) calloc(4, sizeof(double));

    for (i=0; i<nombreSequences; i++)
    {
        nombreDeBaseTotal += (strlen(tableauSequences[i]));
    }


    for (i=0; i<nombreSequences; i++)
    {
        for (j=0; j < (strlen(tableauSequences[i])); j++)
        {
            switch (tableauSequences[i][j])
            {
                case 'A' : motifDeFond[0] += 1; break;
                case 'T' : motifDeFond[1] += 1; break;
                case 'C' : motifDeFond[2] += 1; break;
                case 'G' : motifDeFond[3] += 1; break;
                default : printf("Motif de fond : format de base incorrect.\n"); break;
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
        liberationMemoirePSSM(*p_pssmVide);
    }
    /*ALLOC DOUBLE PUIS REBOUCLE OU JE LAISSE COMME CA ?????????????????*/
    (*p_pssmVide)=(double**)calloc(sizeof(double*),dim1);
    for (i=0; i<dim1; i++)
    {
        (*p_pssmVide)[i]=(double*)calloc(sizeof(double),dim2);
        for (j=0; j<dim2; j++)
        {
            (*p_pssmVide)[i][j]=pssmACopier[i][j];
        }
    }
}

void liberationMemoirePSSM(double **pssm){

    int i;

    for (i = 0; i < 4; i++)
    {
        free(pssm[i]);
    }
    free(pssm);
}

double ** allocDoubleDeuxDim (double **variable, int dim1, int dim2){

    int i;

    variable = (double**) calloc (dim1, sizeof(double*));
    for (i = 0; i < dim1; i++)
    {
        variable[i] = (double*) calloc(dim2, sizeof(double));
    }
    return variable;
}

int ** allocIntDeuxDim (int **variable, int dim1, int dim2){

    int i;

    variable = (int**) calloc (dim1, sizeof(int*));
    for (i = 0;  i < dim1; i++)
    {
        variable[i] = (int*) calloc(dim2, sizeof(int));
    }
    return variable;
}

void copieProfondeInt2D(int ***p_tabIntVide, int **tabIntACopier, int dim1, int dim2){

    int i, j;

    if (*p_tabIntVide != NULL)
    {
        for (i = 0; i < dim1; i++)
        {
            free(p_tabIntVide[0][i]);
        }
    }
    /*ALLOC DOUBLE PUIS REBOUCLE OU JE LAISSE COMME CA ?????????????????*/
    (*p_tabIntVide)=(int**)calloc(sizeof(int*),dim1);
    for (i=0; i<dim1; i++)
    {
        (*p_tabIntVide)[i]=(int*)calloc(sizeof(int),dim2);
        for (j=0; j<dim2; j++)
        {
            (*p_tabIntVide)[i][j]=tabIntACopier[i][j];
        }
    }
}

void copieProfondeTabString(char ***p_tabStringVide, char **tabStringACopier, int dim1, int dim2){

    int i, j;

    if (*p_tabStringVide != NULL)
    {
        for (i = 0; i < dim1; i++)
        {
            free(p_tabStringVide[0][i]);
        }
    }
    /*ALLOC DOUBLE PUIS REBOUCLE OU JE LAISSE COMME CA ?????????????????*/
    (*p_tabStringVide)=(char**)calloc(sizeof(char*),dim1);
    for (i=0; i<dim1; i++)
    {
        (*p_tabStringVide)[i]=(char*)calloc(sizeof(char),dim2);
        for (j=0; j<dim2; j++)
        {
            (*p_tabStringVide)[i][j]=tabStringACopier[i][j];
        }
    }
}