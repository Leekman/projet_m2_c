#include "fonctions.h"

 void recherche_motif (int *masque, int l, int k, double **pssm, char **tableauSequences, int nombreSequences, dictionnaire **p_p_dictionnaire) {
 
	int i,j,m;
	int longueurSequencesCourante;
	char *k_merCourant = NULL;
	double quorum;
	k_mer *k_merCandidat;

	k_mer *pk;
	sequence *ps;
	occurence *po;
	k_mer *nextPk;
	sequence *nextPs;
	occurence *nextPo;

	quorum = 0;
	
	//construit le dictionnaire avec le masque courant
	k_merCourant=(char*)malloc(sizeof(char)*((l-k)+1));

	for (i=0; i<nombreSequences; i++) //Sequence 
	{
		longueurSequencesCourante=(strlen(tableauSequences[i])-1);

		for (j=0; j<((longueurSequencesCourante-l)); j++) //Position dans la séquence
		{
			for (m=0; m<(l-k); m++) //Position dans le k_mer
			{
				k_merCourant[m]=tableauSequences[i][j+masque[m]];								
			}
			k_merCourant[m]='\0';

			ajoutAuDictionnaire(p_p_dictionnaire, k_merCourant, i, j);			
		}
		//printf("%d\n", strlen(tableau_sequences[i]));
	}

	//Verification dictionnaire

	/*pk=(*p_p_dictionnaire)->firstK_mer;
	while (pk != NULL)
	{								
		if (strcmp(pk->k_mer,"AAAAAAAA")==0)
		{
			ps = pk->firstSequence;
			printf("K_mer : %s\n", pk->k_mer);
			while (ps != NULL)
			{
				po = ps->firstOccurence;
				printf("seqence : %d position : ", ps->numSequence+1);
				while (po != NULL)
				{
					printf("%d ", po->position+1);
					po = po->nextOccurence;
				}
				ps = ps->nextSequence;
				printf("\n");
			}
		}
		pk = pk->nextK_mer;
	}*/

	//Selection candidat et calcul quorum
	j=0;

	pk=(*p_p_dictionnaire)->firstK_mer;
	while (pk != NULL)
	{				
		j=0;				
		ps = pk->firstSequence;
		while (ps != NULL)
		{
			ps = ps->nextSequence;
			j++;
			//printf("%d", j);
		}
		if ((double)j/(double)nombreSequences > quorum)
		{
			quorum = (double)j/(double)nombreSequences;
			k_merCandidat = pk;			
		}
		pk = pk->nextK_mer;
	}
	
	//construit la PSSM du motif le plus courant du dictionnaire
	printf("k_merCandidat : %s, quorum : %f\n", k_merCandidat->k_mer, quorum);
	pssm=construirePSSM(k_merCandidat, tableauSequences, nombreSequences, k);

	printf("\nPSSM\n");
	for (i=0;i<4;i++)
	{
		for(j=0;j<(strlen(k_merCandidat->k_mer)+k);j++)
		{
			printf("%f ", pssm[i][j]);
		}
		printf("\n");
	}

	//free du dictionnaire
	pk=(*p_p_dictionnaire)->firstK_mer;
	while (pk != NULL)
	{								
		ps = pk->firstSequence;
		while (ps != NULL)
		{
			po = ps->firstOccurence;
			while (po != NULL)
			{
				nextPo = po->nextOccurence;
				free(po);
				po = nextPo;
			}
			nextPs = ps->nextSequence;
			free(ps);
			ps = nextPs;
		}
		nextPk = pk->nextK_mer;
		free(pk);
		pk = nextPk;
	}
	free(*p_p_dictionnaire);
	(*p_p_dictionnaire)=NULL;

} 

double **construirePSSM(k_mer *k_merCandidat, char **tableauSequences, int nombreSequences, int k){

	double **pssm = NULL;
	sequence *parcoureurSequence;
	occurence *parcoureurOccurence;
	int tailleMotif;

	int nombreDeMotif;

	int i, j, m, n;

	nombreDeMotif = 0;

	pssm = (double**)malloc(sizeof(double*)*4);

	tailleMotif = (strlen(k_merCandidat->k_mer)+k);

	for (i=0; i<4; i++)
	{
		pssm[i]=(double*)malloc(sizeof(double)*tailleMotif);		
	}

	parcoureurSequence = k_merCandidat->firstSequence;

	while (parcoureurSequence != NULL)
		{
			parcoureurOccurence = parcoureurSequence->firstOccurence;
			m=0;
			while (parcoureurOccurence != NULL)
			{
				for (n=0; n<tailleMotif; n++)
				{
					//printf("PLOP %d %d\n", j, (parcoureurOccurence->position)+n);
					//printf("%c", tableauSequences[j][(parcoureurOccurence->position)+n]);
					switch (tableauSequences[parcoureurSequence->numSequence][(parcoureurOccurence->position)+n])
					{
						case 'A' : pssm[0][n] += 1; break;
						case 'T' : pssm[1][n] += 1; break;
						case 'C' : pssm[2][n] += 1; break;
						case 'G' : pssm[3][n] += 1; break;
						default : printf("PSSM : format de base incorect. %c\n", tableauSequences[parcoureurSequence->numSequence][(parcoureurOccurence->position)+n]); break;
					}					
				}
				//printf("\n");				
				parcoureurOccurence = parcoureurOccurence->nextOccurence;
				m++;
			}
			//printf("\nNombre d'occurence dans sequence %d : %d ",parcoureurSequence->numSequence+1, m);
			nombreDeMotif += m;
			parcoureurSequence = parcoureurSequence->nextSequence;
			j++;
		}
		printf("Nombre de sequence : %d. Nombre de motifs : %d\n", j, nombreDeMotif);

	for (i=0;i<4;i++)
	{
		for(j=0;j<(strlen(k_merCandidat->k_mer)+k);j++)
		{
			pssm[i][j] /= nombreDeMotif;
		}
	}

	return pssm;
}