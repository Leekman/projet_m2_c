#include "fonctions.h"

 void recherche_motif (int *masque, int l, int k, double **pssm, char **tableauSequences, int nombreSequences, dictionnaire **p_p_dictionnaire, double *p_score, double *motifDeFond) {
 
	int i,j,m;
	int nbSequenceDuMotif;
	int nbSequenceDuMotifCandidat;
	int longueurSequencesCourante;
	char *k_merCourant = NULL;
	double quorum;
	double scoreCourant;
	double **pssmCourante = NULL;
	k_mer *p_k_merCandidat = NULL;
	int **infoPssmCourante = NULL;
	int nombreOccurence;

	k_mer *pk;
	sequence *parcoureurSequence;
	occurence *p_parcoureurOccurence;
	sequence *ps;
	occurence *po;
	k_mer *nextPk;
	sequence *nextPs;
	occurence *nextPo;

	quorum = 0;
	*p_score = 0;
	
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

	//Calcul du quorum et du score des motif commun, selection du k_mer avec le score le plus haut
	nbSequenceDuMotif=0;

	pk=(*p_p_dictionnaire)->firstK_mer;
	while (pk != NULL)
	{				
		nbSequenceDuMotif=0;				
		parcoureurSequence = pk->firstSequence;
		while (parcoureurSequence != NULL)
		{
			p_parcoureurOccurence = parcoureurSequence->firstOccurence;
			parcoureurSequence = parcoureurSequence->nextSequence;
			nbSequenceDuMotif++;
		}
		if (nbSequenceDuMotif > 30) //Si le motif est présent dans plus de deux séquences
		{
			pssmCourante=construirePSSM(pk, tableauSequences, nombreSequences, k);
			scoreCourant=calculDuScore(pk, tableauSequences, nombreSequences, k, pssmCourante, motifDeFond);
			if (scoreCourant > *p_score)
			{
				quorum = (double)nbSequenceDuMotif/(double)nombreSequences;
				*p_score = scoreCourant;
				copieProfondePSSM(&pssm, pssmCourante, 4, l);
				nbSequenceDuMotifCandidat = nbSequenceDuMotif;
				p_k_merCandidat = pk;
			}
			free(pssmCourante);			
		}
		pk = pk->nextK_mer;
	}

	if (p_k_merCandidat == NULL)
	{
		printf("Aucun motif commun trouvé avec ce masque.\n");

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

		return ;
	}
	//Affichage des infos du k_mer candidat
	printf("p_k_merCandidat : %s, quorum : %f, score : %e\n", p_k_merCandidat->k_mer, quorum, *p_score);
	printf("\nPSSM\n");
	for (i=0;i<4;i++)
	{
		switch (i)
		{
			case 0 : printf("A : "); break;
			case 1 : printf("T : "); break;
			case 2 : printf("C : "); break;
			case 3 : printf("G : "); break;
			default : break;
		}
		for(j=0;j<(strlen(p_k_merCandidat->k_mer)+k);j++)
		{
			printf("%f ", pssm[i][j]);
		}
		printf("\n");
	}

	//Recupération des info du k_mer_candidat dans un tableau

	nombreOccurence = 0;
	pk=p_k_merCandidat;							
	ps = pk->firstSequence;
	while (ps != NULL)
	{
		po = ps->firstOccurence;
		while (po != NULL)
		{
			nextPo = po->nextOccurence;
			po = nextPo;
			nombreOccurence++;
		}
		nextPs = ps->nextSequence;
		ps = nextPs;
	}
	printf("Nombre d'occurence = %d\n", nombreOccurence);

	infoPssmCourante = (int**)malloc(sizeof(int*)*nombreOccurence);
	for (i=0; i<nombreOccurence; i++)
	{
		infoPssmCourante[i] = (int*)malloc(sizeof(int)*2);
	}

	i = 0;
	ps = pk->firstSequence;
	while (ps != NULL)
	{
		po = ps->firstOccurence;
		while (po != NULL)
		{
			infoPssmCourante[i][0] = ps->numSequence;
			infoPssmCourante[i][1] = po->position;
			nextPo = po->nextOccurence;
			po = nextPo;
			i++;
		}
		nextPs = ps->nextSequence;
		ps = nextPs;
	}

	for (i=0; i<nombreOccurence; i++)
	{
		printf("%d %d\n",infoPssmCourante[i][0], infoPssmCourante[i][1]);
	}

	//Amélioration du motif

	ameliorerMotif(infoPssmCourante, pssm, p_score, tableauSequences, nombreSequences, nombreOccurence, k, l, p_k_merCandidat, motifDeFond); 

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

double **construirePSSM(k_mer *p_k_merCandidat, char **tableauSequences, int nombreSequences, int k){

	double **pssm = NULL;
	sequence *parcoureurSequence;
	occurence *parcoureurOccurence;
	int tailleMotif;

	int nombreDeMotif;

	int i, j, m, n;

	j = 0;

	nombreDeMotif = 0;

	pssm = (double**)malloc(sizeof(double*)*4);

	tailleMotif = (strlen(p_k_merCandidat->k_mer)+k);

	for (i=0; i<4; i++)
	{
		pssm[i]=(double*)malloc(sizeof(double)*tailleMotif);		
	}

	parcoureurSequence = p_k_merCandidat->firstSequence;

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
		//printf("k_mer : %s. Nombre de sequence : %d. Nombre de motifs : %d\n",p_k_merCandidat->k_mer, j, nombreDeMotif);

	for (i=0;i<4;i++)
	{
		for(j=0;j<(strlen(p_k_merCandidat->k_mer)+k);j++)
		{
			pssm[i][j] /= nombreDeMotif;
		}
	}

	return pssm;
}

double calculDuScore(k_mer *p_k_merCandidat, char **tableauSequences, int nombreSequences, int k, double **pssm, double *motifDeFond){

	double score;

	score = 0;

	score = calculScoreK_mer(p_k_merCandidat, tableauSequences, k, pssm, motifDeFond);
	//calculProbaMotifDeFond(p_k_merCandidat, tableauSequences, k, motifDeFond);

	return score;
}

double calculScoreK_mer(k_mer *p_k_merCandidat, char **tableauSequences, int k, double **pssm, double *motifDeFond){

	double probaPSSM;
	double probaMotifDeFond;
	double scoreK_mer;
	double scoreOccurence;
	int i;
	int longueurMotif;

	sequence *parcoureurSequence;
	occurence *parcoureurOccurence;

	scoreK_mer = 0;

	longueurMotif = (strlen(p_k_merCandidat->k_mer)+k);
	
	parcoureurSequence = p_k_merCandidat->firstSequence;

	while (parcoureurSequence != NULL)
	{
		parcoureurOccurence = parcoureurSequence->firstOccurence;
		while (parcoureurOccurence != NULL)
		{

			probaPSSM = 1;
			probaMotifDeFond = 1;
			for (i=0; i<longueurMotif; i++)
			{
				
				switch (tableauSequences[parcoureurSequence->numSequence][i+parcoureurOccurence->position])
				{
					case 'A' : probaPSSM *= pssm[0][i]; probaMotifDeFond *= motifDeFond[0]; break;
					case 'T' : probaPSSM *= pssm[1][i]; probaMotifDeFond *= motifDeFond[0]; break;
					case 'C' : probaPSSM *= pssm[2][i]; probaMotifDeFond *= motifDeFond[0]; break;
					case 'G' : probaPSSM *= pssm[3][i]; probaMotifDeFond *= motifDeFond[0]; break;
					default : printf("calculScoreK_mer : format de base incorrect.\n"); break;
				}

			}
			//printf("proba PSSM : %e\n", probaPSSM);
			//printf("proba motif : %e\n", probaMotifDeFond);
			scoreOccurence = probaPSSM/probaMotifDeFond;
			scoreK_mer += log(scoreOccurence);
			//printf("Quotient : %e\n\n", scoreK_mer);
			parcoureurOccurence = parcoureurOccurence->nextOccurence;
		}
		parcoureurSequence = parcoureurSequence->nextSequence;
	}
	
	return scoreK_mer;
}

void ameliorerMotif(int **infoPssmCourante, double **pssmCourante, double *p_scoreCourant, char **tableauSequences,int nombreSequences, int nombreOccurence, int k, int l, k_mer *p_k_merCandidat, double *motifDeFond){

	double critereDeConvergene;
	int randomOccurence;
	int randomPosition;
	double **pssmNouvelle;
	double scoreNouveau;
	int i,j, compteur;

	critereDeConvergene = 2;
	compteur = 0;


	while (critereDeConvergene > 1.5)
	{
		critereDeConvergene = 0;
		compteur++;

		randomOccurence = rand()%(nombreOccurence);
		printf("randomOccurence : %d\n", randomOccurence);
		do
		{
			randomPosition = rand()%(strlen(tableauSequences[infoPssmCourante[randomOccurence][0]])-l);
			printf("randomPosition : %d\n", randomPosition);
		}
		while (randomPosition == infoPssmCourante[randomOccurence][1]);

		pssmNouvelle = (double**)malloc(sizeof(double*)*4);
		for (i = 0; i < 4; i++)
		{
			pssmNouvelle[i]=(double*)malloc(sizeof(double)*l);
			for (j = 0; j < l; j++)
			{
				pssmNouvelle[i][j] = 0;
			}
		}

		for (i = 0; i < nombreOccurence; i++)
		{
			if (i != randomOccurence)
			{
				for (j = 0; j < l; j++)
				{
					switch (tableauSequences[infoPssmCourante[i][0]][infoPssmCourante[i][1]+j])
					{
						case 'A' : pssmNouvelle[0][j] += 1; break;
						case 'T' : pssmNouvelle[1][j] += 1; break;
						case 'C' : pssmNouvelle[2][j] += 1; break;
						case 'G' : pssmNouvelle[3][j] += 1; break;
						default : printf("PSSM : format de base incorect. %c\n", tableauSequences[infoPssmCourante[i][0]][infoPssmCourante[i][1]+j]); break;
					}
				}
			}
			else
			{
				for (j = 0; j < l; j++)
				{
					switch (tableauSequences[infoPssmCourante[i][0]][randomPosition+j])
					{
						case 'A' : pssmNouvelle[0][j] += 1; break;
						case 'T' : pssmNouvelle[1][j] += 1; break;
						case 'C' : pssmNouvelle[2][j] += 1; break;
						case 'G' : pssmNouvelle[3][j] += 1; break;
						default : printf("PSSM : format de base incorect. %c\n", tableauSequences[infoPssmCourante[i][0]][randomPosition+j]); break;
					}
				}
			}
		}

		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < l; j++)
			{
				pssmNouvelle[i][j] /= nombreOccurence;
				critereDeConvergene += ((pssmNouvelle[i][j]-pssmCourante[i][j])*(pssmNouvelle[i][j]-pssmCourante[i][j]));
			}
		}

		printf("Critere de convergence : %f\n", critereDeConvergene);

		scoreNouveau=calculDuScore(p_k_merCandidat, tableauSequences, nombreSequences, k, pssmNouvelle, motifDeFond);

		if (scoreNouveau <= *p_scoreCourant)
		{
			printf("Arret amelioration apres %d iteration car le score n'est pas ameliore.\n", compteur);
			return ;
		}

		*p_scoreCourant = scoreNouveau;

		copieProfondePSSM(&pssmCourante, pssmNouvelle, 4, l);

		infoPssmCourante[randomOccurence][1]=randomPosition;
	}

	printf("Arret amelioration apres %d iteration car le critere de convergence n'est pas depasse.\n", compteur);

}