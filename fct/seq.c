#include "../lib/includes.h"
#include "../lib/seq.h"
#include "../lib/doublon.h"
#define PROPORTION 20


//////////////////////////////////////////////////////////////
/*Fonction permettant de générer aléatoirement un nucléotide*/
//////////////////////////////////////////////////////////////
char genererNuc(){
    char nucleotide[5]={"ATCG"};
    return nucleotide[rand()%4];
}

/////////////////////////////////////////////////////////////
/*Fonction permettant de créer une séquences de nucléotides*/
/////////////////////////////////////////////////////////////

char *sequenceSeule(int tailleSeq){
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
/////////////////////////////////////////////////////////////////////////////////////////////
/*Fonction permettant de modifier le motif en fonction du nombre d'erreur maximum autorisée*/
/////////////////////////////////////////////////////////////////////////////////////////////
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


///////////////////////////////////////////////////////////////////////////////////
/*Procédure permettant l'insertion du motif dans la séquence généré préalablement*/
///////////////////////////////////////////////////////////////////////////////////
void insertionMotif(char *motif, int positionMotif, int tailleMotif,char *sequence){
    int i,j;
    j=0;
    for (i = positionMotif; i <(positionMotif+tailleMotif) ; i++)
    {
        sequence[i]=motif[j];
        j++;
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Procédure permettant la création des séquences avec l'insertion des motifs avec taille randomisée et position du motif randomisée*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void creationSeq(int nbErreurMax, int nbSeq, int tailleSeq, char *motif, char **tabSeq, int *tabPosition, int *tabNbErreur){

    int i;
    int positionMotif, tailleMotif, positionMotifMax;
    int tailleTemp, intervalle;
    char *motifFixe;

    tailleMotif=strlen(motif);
    motifFixe=(char*)malloc((tailleMotif+1)*sizeof(char));
    intervalle=(int)(tailleSeq*PROPORTION)/100;

    for(i=0 ; i < nbSeq ; i++)
        {
         strcpy(motifFixe,motif);
         tailleTemp=tailleSeq+rand()%((intervalle*2)+1)-intervalle;
         positionMotifMax=(--tailleTemp-tailleMotif);
         positionMotif=rand()%positionMotifMax;
         tabPosition[i]=positionMotif;
         tabNbErreur[i]=modifMotif(nbErreurMax,motifFixe,tailleMotif);
         tabSeq[i]=sequenceSeule(tailleTemp);
         insertionMotif(motifFixe,positionMotif,tailleMotif,tabSeq[i]);
        }

}
///////////////////////////////////////////////////////////////////////////////////
/*Procédure permettant de créer un fichier fasta contenant les séquences générées*/
///////////////////////////////////////////////////////////////////////////////////
void insertFasta(char **tabSeq, int nbSeq){

    FILE* fasta=NULL;
    int i;
    fasta=fopen("../output/sequences.fasta","w");

    if (fasta != NULL)
    {
        i=0;
        for (i=0; i<nbSeq ; i++)
        {
            fprintf(fasta,">seq%d\n%s\n",i+1,tabSeq[i]); //affichage tableau
        }
        fclose(fasta);    
    }
}