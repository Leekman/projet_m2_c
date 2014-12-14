#include "../lib/motif.h"

 void recherche_motif (int *masque, int l, int k, double **pssm, int ***infoPssmCourante, char **tableauSequences, int nombreSequences, dictionnaire **p_p_dictionnaire, double *p_score, double *motifDeFond, char ***p_ensembleT, double ***p_motifConsensusPSSM, char **p_motifConsensus, int *nbSequenceDuMotifCandidat, int *p_scoreMasque, int nbErreurMax, int iteration){
 
	int i,j,m;
	int nbSequenceDuMotif;
	int longueurSequencesCourante;
	char *k_merCourant = NULL;
	//double quorum;
	double scoreCourant;
	double **pssmCourante = NULL;
	k_mer *p_k_merCandidat = NULL;
	int nombreOccurence;

k_mer *nextPk = NULL;	
	k_mer *pk = NULL;
	sequence *parcoureurSequence = NULL;
	sequence *ps = NULL;
	occurence *po = NULL;
	sequence *nextPs = NULL;
	occurence *nextPo = NULL;

	//quorum = 0;
	*p_score = 0;

	////////////////////////////////////////////////////
	/*CONSTRUIT LE DICTIONNAIRE AVEC LE MASQUE COURANT*/
	////////////////////////////////////////////////////
	for (i=0; i<nombreSequences; i++) //Sequence 
	{
		longueurSequencesCourante=(strlen(tableauSequences[i])-1);

		for (j=0; j<((longueurSequencesCourante-l)); j++) //Position dans la séquence
		{
			k_merCourant=(char*)malloc(sizeof(char)*((l-k)+1));
			for (m=0; m<(l-k); m++) //Position dans le k_mer
			{
				k_merCourant[m]=tableauSequences[i][j+masque[m]];								
			}
			k_merCourant[m]='\0';
			ajoutAuDictionnaire(p_p_dictionnaire, k_merCourant, i, j);
			free(k_merCourant);			
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	/*CALCUL DU QUORUM ET DU SCORE DES MOTIFS COMMUNS, SELECTION DU K MER AVEC LE SCORE LE PLUS HAUT*/
	//////////////////////////////////////////////////////////////////////////////////////////////////

	nbSequenceDuMotif=0;

	pk=(*p_p_dictionnaire)->firstK_mer;
	while (pk != NULL)
	{	
		nbSequenceDuMotif=0;				
		parcoureurSequence = pk->firstSequence;
		while (parcoureurSequence != NULL)
		{
			parcoureurSequence = parcoureurSequence->nextSequence;
			nbSequenceDuMotif++;
		}
		if ((double)nbSequenceDuMotif/(double)nombreSequences > 0.30) //Si le motif est présent dans plus de deux séquences ???????
		{
			pssmCourante=construirePSSM(pk, tableauSequences, nombreSequences, k);			
			scoreCourant=calculDuScore(pk, tableauSequences, k, pssmCourante, motifDeFond);
			if (scoreCourant > *p_score)
			{
				//quorum = (double)nbSequenceDuMotif/(double)nombreSequences;
				*p_score = scoreCourant;
				copieProfondePSSM(&pssm, pssmCourante, 4, l);
				*nbSequenceDuMotifCandidat = nbSequenceDuMotif;
				p_k_merCandidat = pk;
			}
			liberationMemoirePSSM(pssmCourante);		
		}
		pk = pk->nextK_mer;
	}

	if (p_k_merCandidat == NULL)
	{
		//printf("Aucun motif commun trouvé avec le masque : %d\n", iteration);
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
			free(pk->k_mer);
			free(pk);
			pk = nextPk;
		}
		free(*p_p_dictionnaire);
		(*p_p_dictionnaire)=NULL;
		//liberationDictionnaire(*p_p_dictionnaire);
		return ;
	}

	///////////////////////////////////////////////
	/*AFFICHAGE DES INFORMATIONS DU KMER CANDIDAT*/
	///////////////////////////////////////////////

	//printf("\nInformations sur le kmer Candidat : %s, quorum : %f, score : %e\n", p_k_merCandidat->k_mer, quorum, *p_score);

	/////////////////////////////////////////////////////////
	/*RECUPERATION DES INFORMATIONS DU KMER DANS UN TABLEAU*/
	/////////////////////////////////////////////////////////

	nombreOccurence = recupNombreOccurence(p_k_merCandidat);

	infoPssmCourante[0] = allocIntDeuxDim(infoPssmCourante[0], nombreOccurence, 2);

	i = 0;

	pk = p_k_merCandidat;
	ps = pk->firstSequence;
	while (ps != NULL)
	{
		po = ps->firstOccurence;
		while (po != NULL)
		{
			infoPssmCourante[0][i][0] = ps->numSequence;
			infoPssmCourante[0][i][1] = po->position;
			nextPo = po->nextOccurence;
			po = nextPo;
			i++;
		}
		nextPs = ps->nextSequence;
		ps = nextPs;
	}

	/////////////////////////
	/*AMELIORATION DU MOTIF*/
	/////////////////////////

	for (i = 0; i < 20; i++)
	{
		ameliorerMotif(infoPssmCourante, pssm, p_score, tableauSequences, nombreSequences, nombreOccurence, k, l, p_k_merCandidat, motifDeFond); 
	}

	/////////////////////
	/*AFFINAGE DU MOTIF*/
	/////////////////////

	affinerMotif(p_ensembleT, infoPssmCourante, tableauSequences, nombreOccurence, l, pssm, motifDeFond, nbSequenceDuMotifCandidat);
	creerMotifConsensus(p_motifConsensusPSSM, p_motifConsensus, *p_ensembleT, l, nbSequenceDuMotifCandidat);

	///////////////////
	/*SCORE DU MASQUE*/
	///////////////////

	*p_scoreMasque = calculScoreMasque(nbErreurMax, l, *p_motifConsensus, *p_ensembleT, nbSequenceDuMotifCandidat);
	
	//////////////////////////////
	/*LIBERATION DU DICTIONNAIRE*/
	//////////////////////////////

	//liberationDictionnaire(*p_p_dictionnaire);

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
			free(pk->k_mer);
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
	tailleMotif = (strlen(p_k_merCandidat->k_mer)+k);

	pssm = allocDoubleDeuxDim(pssm, 4, tailleMotif);

	

	parcoureurSequence = p_k_merCandidat->firstSequence;

	while (parcoureurSequence != NULL)
		{
			parcoureurOccurence = parcoureurSequence->firstOccurence;
			m=0;
			while (parcoureurOccurence != NULL)
			{
				for (n = 0; n < tailleMotif; n++)
				{
					switch (tableauSequences[parcoureurSequence->numSequence][(parcoureurOccurence->position)+n])
					{
						case 'A' : pssm[0][n] += 1; break;
						case 'T' : pssm[1][n] += 1; break;
						case 'C' : pssm[2][n] += 1; break;
						case 'G' : pssm[3][n] += 1; break;
						default : printf("PSSM : Format de base incorrect. %c\n", tableauSequences[parcoureurSequence->numSequence][(parcoureurOccurence->position)+n]); break;
					}					
				}			
				parcoureurOccurence = parcoureurOccurence->nextOccurence;
				m++;
			}
			nombreDeMotif += m;
			parcoureurSequence = parcoureurSequence->nextSequence;
			j++;
		}

	for (i = 0; i < 4; i++)
	{
		for(j = 0; j < (strlen(p_k_merCandidat->k_mer)+k); j++)
		{
			pssm[i][j] /= nombreDeMotif;
		}
	}
	return pssm;
}

double calculDuScore(k_mer *p_k_merCandidat, char **tableauSequences, int k, double **pssm, double *motifDeFond){

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
					case 'T' : probaPSSM *= pssm[1][i]; probaMotifDeFond *= motifDeFond[1]; break;
					case 'C' : probaPSSM *= pssm[2][i]; probaMotifDeFond *= motifDeFond[2]; break;
					case 'G' : probaPSSM *= pssm[3][i]; probaMotifDeFond *= motifDeFond[3]; break;
					default : printf("calculDuScore : format de base incorrect.\n"); break;
				}

			}
			scoreOccurence = probaPSSM/probaMotifDeFond;
			scoreK_mer += log(scoreOccurence);
			parcoureurOccurence = parcoureurOccurence->nextOccurence;
		}
		parcoureurSequence = parcoureurSequence->nextSequence;
	}
	return scoreK_mer;
}


void ameliorerMotif(int ***infoPssmCourante, double **pssmCourante, double *p_scoreCourant, char **tableauSequences,int nombreSequences, int nombreOccurence, int k, int l, k_mer *p_k_merCandidat, double *motifDeFond){

	double critereDeConvergence;
	int randomOccurence;
	int randomPosition;
	double **pssmNouvelle;
	double scoreNouveau;
	int i,j, compteur;

	critereDeConvergence = 2;
	compteur = 0;


	while (critereDeConvergence > 0.005)
	{
		critereDeConvergence = 0;
		compteur++;

		randomOccurence = rand()%(nombreOccurence);
		do
		{
			randomPosition = rand()%(strlen(tableauSequences[infoPssmCourante[0][randomOccurence][0]])-l);
		}
		while (randomPosition == infoPssmCourante[0][randomOccurence][1]);

		pssmNouvelle = allocDoubleDeuxDim (pssmNouvelle, 4, l);

		for (i = 0; i < nombreOccurence; i++)
		{
			if (i != randomOccurence)
			{
				for (j = 0; j < l; j++)
				{
					switch (tableauSequences[infoPssmCourante[0][i][0]][infoPssmCourante[0][i][1]+j])
					{
						case 'A' : pssmNouvelle[0][j] += 1; break;
						case 'T' : pssmNouvelle[1][j] += 1; break;
						case 'C' : pssmNouvelle[2][j] += 1; break;
						case 'G' : pssmNouvelle[3][j] += 1; break;
						default : printf("PSSM : format de base incorect. %c\n", tableauSequences[infoPssmCourante[0][i][0]][infoPssmCourante[0][i][1]+j]); break;
					}
				}
			}
			else
			{
				for (j = 0; j < l; j++)
				{
					switch (tableauSequences[infoPssmCourante[0][i][0]][randomPosition+j])
					{
						case 'A' : pssmNouvelle[0][j] += 1; break;
						case 'T' : pssmNouvelle[1][j] += 1; break;
						case 'C' : pssmNouvelle[2][j] += 1; break;
						case 'G' : pssmNouvelle[3][j] += 1; break;
						default : printf("PSSM : format de base incorect. %c\n", tableauSequences[infoPssmCourante[0][i][0]][randomPosition+j]); break;
					}
				}
			}
		}

		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < l; j++)
			{
				pssmNouvelle[i][j] /= nombreOccurence;
				critereDeConvergence += ((pssmNouvelle[i][j]-pssmCourante[i][j])*(pssmNouvelle[i][j]-pssmCourante[i][j]));
			}
		}

		scoreNouveau=calculDuScore(p_k_merCandidat, tableauSequences, k, pssmNouvelle, motifDeFond);

		if (scoreNouveau <= *p_scoreCourant)
		{
			liberationMemoirePSSM(pssmNouvelle);
			return ;
		}

		*p_scoreCourant = scoreNouveau;

		copieProfondePSSM(&pssmCourante, pssmNouvelle, 4, l);

		infoPssmCourante[0][randomOccurence][1]=randomPosition;

		printf("Le motif à été amélioré.\n");
	}
}

void affinerMotif(char ***p_ensembleT, int ***infoPssmCourante, char **tableauSequences, int nombreOccurence, int l, double **pssm, double *motifDeFond, int *nbSequenceDuMotifCandidat){

	int i, j, m;
	double probaPSSM, probaMotifDeFond, scoreMotifCourant, scoreMax;
	int meilleurMotif;

	p_ensembleT[0] = (char**)malloc(sizeof(char*)*(*nbSequenceDuMotifCandidat+1));
	for (i = 0; i < *nbSequenceDuMotifCandidat; i++)
	{
		p_ensembleT[0][i] = (char*)malloc(sizeof(char)*(l+1));
	}

	for (i = 0; i < nombreOccurence; i++)
	{
		scoreMax = 0;

		///////////////////////////////////////////////
		/*IDENTIFICATION DU MOTIF MAXIMISANT LE SCORE*/
		///////////////////////////////////////////////

		for (j = 0; j < strlen(tableauSequences[infoPssmCourante[0][i][0]])-l; j++)
		{
			probaPSSM = 1;
			probaMotifDeFond = 1;
			for (m = 0; m < l; m++)
			{
				switch (tableauSequences[infoPssmCourante[0][i][0]][j+m])
				{
					case 'A' : probaPSSM *= pssm[0][m]; probaMotifDeFond *= motifDeFond[0]; break;
					case 'T' : probaPSSM *= pssm[1][m]; probaMotifDeFond *= motifDeFond[1]; break;
					case 'C' : probaPSSM *= pssm[2][m]; probaMotifDeFond *= motifDeFond[2]; break;
					case 'G' : probaPSSM *= pssm[3][m]; probaMotifDeFond *= motifDeFond[3]; break;
					default : printf("calculDuScore : format de base incorrect.\n"); break;
				}

			}
			scoreMotifCourant = probaPSSM/probaMotifDeFond;				

			if (scoreMotifCourant > scoreMax)
			{
				scoreMax = scoreMotifCourant;
				meilleurMotif = j;
				infoPssmCourante[0][i][1] = meilleurMotif;
			}
		}

		///////////////////////////////
		/*REMPLISSAGE DE L'ENSEMBLE T*/
		///////////////////////////////

		for (j = 0; j < l; j++)
		{
			p_ensembleT[0][i][j]=tableauSequences[infoPssmCourante[0][i][0]][meilleurMotif+j];	
		}
		p_ensembleT[0][i][j] = '\0';
	}
	liberationMemoirePSSM(pssm);
}

void creerMotifConsensus(double ***p_motifConsensusPSSM, char **p_motifConsensus, char **ensembleT, int l, int *nbSequenceDuMotifCandidat){

	int i, j, max;

	p_motifConsensus[0] = (char*)malloc(sizeof(char)*(l+1));

	p_motifConsensusPSSM[0] = allocDoubleDeuxDim (p_motifConsensusPSSM[0], 4, l); 

	for (i = 0; i < *nbSequenceDuMotifCandidat; i++)
	{
		for (j = 0; j < l; j++)
		{
			switch (ensembleT[i][j])
			{
				case 'A' : p_motifConsensusPSSM[0][0][j] +=1; break;
				case 'T' : p_motifConsensusPSSM[0][1][j] +=1; break;
				case 'C' : p_motifConsensusPSSM[0][2][j] +=1; break;
				case 'G' : p_motifConsensusPSSM[0][3][j] +=1; break;
			}
		}		
	}
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < l; j++)
		{
			p_motifConsensusPSSM[0][i][j] /= *nbSequenceDuMotifCandidat;
		}
	}

	for (i = 0; i < l; i++)
	{
		max = 0;
		for (j = 0; j < 4; j++)
		{
			if (p_motifConsensusPSSM[0][j][i] > p_motifConsensusPSSM[0][max][i])
			{
				max = j;
			}
		}
		switch (max)
		{
			case 0 : p_motifConsensus[0][i] = 'A'; break;
			case 1 : p_motifConsensus[0][i] = 'T'; break;
			case 2 : p_motifConsensus[0][i] = 'C'; break;
			case 3 : p_motifConsensus[0][i] = 'G'; break;
		}
	}
	p_motifConsensus[0][i] = '\0';
}


int calculScoreMasque(int nbErreurMax, int l, char *motifConsensus, char **ensembleT, int *nbSequenceDuMotifCandidat){

	int i, j;
	int nbDifferences;
	int score;

	score = 0;

	for (i = 0; i < *nbSequenceDuMotifCandidat; i++)
	{
		nbDifferences = 0;

		for (j = 0; j < l; j++)
		{
			if (ensembleT[i][j] != motifConsensus[j])
			{
				nbDifferences += 1;
			}
		}
		if (nbDifferences > nbErreurMax)
		{
			score += 1;
		}
	}

	return score;
}