#include "../lib/includes.h"
#include "../lib/param.h"


/******************/
/* Affiche l'aide */
/******************/
void notice(){
	printf("Utilisation : creation_sequence -d entier -n entier -t entier -m motif\n\
Description :\n\
Exemple : creation_sequence -d 2 -n 50 -t 100 -m TATATGCA\n\
OPTIONS :\n\
	-d --nbErreur	entier : nombre de substitutions maximum autorisees\n\
	-n --nbSeq	entier : nombre de sequences a creer\n\
	-t --tailleSeq	entier : taille des sequences a creer\n\
	-m --motif	chaine de caractere : motif a inserer dans la sequence\n\
	\n\
	-h --help	affiche cette aide\n\
	\n");
}



/*********************/
/* liste des options */
/*********************/

struct option long_options[] = {
		{"nbErreur", required_argument, NULL, 'd'},
		{"nbSeq", required_argument, NULL, 'n'},
		{"tailleSeq", required_argument, NULL, 't'},
		{"motif", required_argument, NULL, 'm'},
		{"help", no_argument, 0,  'h' }
};

///////////////////////////////////////////////////
/*Fonction permettant de récupérer les paramètres*/
///////////////////////////////////////////////////

//idée d'amélioration faire une fonction a la place de la procédure qui renverrai le motif. Permet de faire l'allocation dynamique dans la fonction et d'éviter les fuites de mémoires


char* getParam (int* nbErreur, int* nbSeq, int* tailleSeq, int argc, char *argv[]){
	char *motif = NULL;
	int opt = 0;
	int long_index = 0; /* index des options */
	if (argc<2)
	{
		fprintf(stderr, "ERREUR\n"); 
		notice(); 
		exit(1);
	}
	regex_t preg;
	regcomp(&preg,"[^ATCG]+",REG_EXTENDED|REG_NOSUB);
	while ((opt = getopt_long(argc, argv, "d:n:t:m:h", long_options, &long_index )) != -1) {
		switch (opt) {
			case 'd' : *nbErreur = atoi(optarg);
				break;
			case 'n' : *nbSeq = atoi(optarg);
				break;
			case 't' : *tailleSeq = atoi(optarg);
				break;
			case 'm' : 	motif=(char*)malloc(((strlen(optarg))+1)*sizeof(char));
						strcpy(motif, optarg);
				break;
			case 'h' : notice(); exit(0);
				break;
			default  : fprintf(stderr, "Mauvais argument\n"); notice(); exit(1);
				break;
		}
	}
	if (!regexec(&preg,motif,0,0,0))
	{
		printf("[ERREUR] Le motif ne peux contenir que les caracteres A,T,C ou G\n");
		notice();
		exit(1);
	}
	regfree(&preg);
	return motif;
}

