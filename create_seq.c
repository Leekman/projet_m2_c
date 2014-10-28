#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

bool estDoublon(int n, int *doublon, int count){

    int i;

    for (i=0; i<count; i++) 
    {
        if (n==doublon[i])
        {
            return true;
        }
    }
    return false;
}


char genererNuc(){
	char nucleotide[5]={"ATCG"};
	return nucleotide[rand()%4];
}


int modifMotif(int nbErreurMax,char *motif, int tailleMotif){

    int nbErreur, positionChangement;
    int i;
    int *tabErreur;
    char *nuclActuel;
    nbErreur=rand()%(++nbErreurMax); //permet de faire le rand entre 0 et le nombre d'erreur maximum (car la borne supérieure est exclue)$
    tabErreur=(int*)malloc(nbErreur*sizeof(int));
    nuclActuel=(char*)malloc(nbErreur*sizeof(char));
    for (i = 0; i < nbErreur; i++)
    {
        do
        {
            positionChangement=rand()%tailleMotif;
            tabErreur[i]=positionChangement;
            nuclActuel[i]=motif[positionChangement]; //stock le nucleotide correspondant pour être sûr de ne pas remplacer ce nucl par le meme
        } while (estDoublon(tabErreur[i], tabErreur, i));
        do  
        {
            motif[tabErreur[i]]=genererNuc();     //assure nucl non identique
        } while (nuclActuel[i]==motif[tabErreur[i]]);  
    }
    free(tabErreur);
    free(nuclActuel);
    return nbErreur;
}

char *createSeq(int tailleSeq){
    int i;
    char* sequence=NULL;//initialisation du premier tableau de caractère
    sequence=(char*) malloc((tailleSeq+1)*sizeof(char));
    for (i=0;i<tailleSeq;i++)
    {
        sequence[i]=genererNuc();
    }
    sequence[i]='\0';
    return sequence;
}

void insertionMotif(char *motif, int positionMotif, int tailleMotif,char *sequence){
	int i,j;
	j=0;
	for (i = positionMotif; i <(positionMotif+tailleMotif) ; i++)
	{
		sequence[i]=motif[j];
		j++;
	}
}


int main(int argc, char *argv[]) {
    
    srand(time(NULL));
    int positionMotif, tailleMotif, tailleMotifMax, positionMotifMax, nbErreurMax;
    int tailleSeq, nbSeq, tailleTemp, proportionTailleSeq, intervalle;
    int i;
    char *motif;
    char *motifFixe;
    char **tabSeq;
    int *tabNbErreur;
    int *tabPosition;
    FILE* fasta=NULL;
    tailleSeq=500;
    nbSeq=100;
    tailleMotifMax=(int)(tailleSeq*10)/100;
    motif=(char*)malloc((tailleMotifMax+1)*sizeof(char));// perte de mémoire si motif plus petit mais c'est ça ou alloc statique avec grosse perte. 
    strcpy(motif,"MMMMMMMMMMMMMMM");
    tailleMotif=strlen(motif);
    motifFixe=(char*)malloc((tailleMotif+1)*sizeof(char));
    tabNbErreur=(int*)malloc(nbSeq*sizeof(int));
    tabPosition=(int*)malloc(nbSeq*sizeof(int));
    nbErreurMax=3;
    proportionTailleSeq=20;
    intervalle=(int)(tailleSeq*proportionTailleSeq)/100;
    //positionMotifMax=(--tailleSeq-tailleMotif);
    tabSeq =(char**) malloc(nbSeq * sizeof(char*));//tableau a deux dimensions de séquences



    for(i=0 ; i < nbSeq ; i++)
    {
    	strcpy(motifFixe,motif);
    	tailleTemp=tailleSeq+rand()%((intervalle*2)+1)-intervalle;
    	positionMotifMax=(--tailleTemp-tailleMotif);
    	positionMotif=rand()%positionMotifMax;
    	tabPosition[i]=positionMotif;
    	tabNbErreur[i]=modifMotif(nbErreurMax,motifFixe,tailleMotif);
    	printf("%d\n", tabNbErreur[i]);
    	tabSeq[i]=createSeq(tailleTemp);
    	insertionMotif(motifFixe,positionMotif,tailleMotif,tabSeq[i]);
    
    }

    fasta=fopen("output/sequences.fasta","w");

    if (fasta != NULL)
    {
    	i=0;
   		for (i=0; i<nbSeq ; i++)
    	{
        	fprintf(fasta,">seq%d\n%s\n",i+1,tabSeq[i]); //affichage tableau
    	}
    	fclose(fasta);    }
   


	return 0;
}
