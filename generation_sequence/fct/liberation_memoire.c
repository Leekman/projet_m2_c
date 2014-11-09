#include "../lib/includes.h"

//////////////////////////////////////////////////////////
/*PROCEDURE PERMETTANT DE LIBERER LE TABLEAU DE SEQUENCE*/
//////////////////////////////////////////////////////////

void freeTabSeq(char **tabSeq, int nbSeq){

	int i;

	for (i = 0; i < nbSeq; i++)
	{
		free(tabSeq[i]);
	}
	free(tabSeq);
}

///////////////////////////////////////////
/*PROCEDURE PERMETTANT DE LIBERER LA PSSM*/
///////////////////////////////////////////

void freePSSM(double **PSSM){

	int i;

	for (i = 0; i < 4; i++)
	{
		free(PSSM[i]);
	}
	free(PSSM);

}
