#include "lib/includes.h"
#include "lib/doublon.h"
#include "lib/seq.h"
#include "lib/param.h"
int main(int argc, char *argv[]) {
    
    
    ////////////////////////////////////////////////
    /*Initialisation du rand au time du processeur*/
    ////////////////////////////////////////////////
    srand(time(NULL));

    ////////////////////////////////
    /*DECLARATION ET INITIALISATION DES VARIABLES*/
    ////////////////////////////////
    int tailleSeq, nbSeq, nbErreurMax;
    int tailleMotif;
    tailleSeq = 0;
    nbSeq = 0;
    nbErreurMax = 0;
    int *tabNbErreur=NULL;
    int *tabPosition=NULL;
	int *p_tailleSeq = NULL;
	int *p_nbSeq = NULL;
    int *p_nbErreurMax = NULL;
    char motif[1000];//allocation dynamique serait mieux mais je ne sais pas comment faire
    char **tabSeq = NULL;
    double **PSSM = NULL;
    p_nbSeq = &nbSeq;
    p_tailleSeq = &tailleSeq;
    p_nbErreurMax = &nbErreurMax;
    //FILE* donnees=NULL;
    tailleMotif=strlen(motif);

    ///////////////////////////////
    /*RECUPERATION DES PARAMETRES*/
    ///////////////////////////////

    getParam(p_nbErreurMax, p_nbSeq, p_tailleSeq, motif, argc, argv);
    //////////////////////////////////////
    /*ALLOCATION DES DIFFERENTS TABLEAUX*/
    //////////////////////////////////////
    tabNbErreur=(int*)malloc(nbSeq*sizeof(int));
    tabPosition=(int*)malloc(nbSeq*sizeof(int));
    tabSeq =(char**) malloc(nbSeq * sizeof(char*));

    /////////////////////////////////////
    /*CREATION DES SEQUENCES AVEC MOTIF*/
    /////////////////////////////////////

    creationSeq(nbErreurMax, nbSeq, tailleSeq, motif, tabSeq, tabPosition, tabNbErreur);


    ///////////////////////
    /*CREATION DE LA PSSM*/
    ///////////////////////

    PSSM=construirePSSM(tailleMotif, tabSeq, nbSeq,tabPosition);

    /////////////////////////////////////////////////
    /*INSERTION DES SEQUENCES DANS UN FICHIER FASTA*/
    /////////////////////////////////////////////////

    creationFasta(tabSeq, nbSeq);


    /////////////////////////////////////
    /*CREATION DU DEUXIEME FICHIER INFO*/
    /////////////////////////////////////

    creationInfo(PSSM, motif, tabPosition, tabNbErreur, tailleSeq, tailleMotif, nbSeq, nbErreurMax);

    ////////////////////////////
    /*LIBERATION DE LA MEMOIRE*/
    ////////////////////////////

    free(tabNbErreur);
    free(tabPosition);
    free(PSSM);
    free(tabSeq);

	return 0;
}
