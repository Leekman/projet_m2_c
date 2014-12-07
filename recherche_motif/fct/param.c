#include "../lib/param.h"


/******************/
/* Affiche l'aide */
/******************/
void notice(){
	printf("Utilisation : ./bin/recherche_motif -l entier -k entier -c chaine de caractere\n\
Description :\n\
Exemple : ./bin/recherche_motif -l 4 -k 2 -c input/sequences.fasta\n\
OPTIONS :\n\
	-l --longueurMotif	entier : Longueur du motif commun a identifier\n\
	-k --nbFenetre	entier : nombre de fenetres des masques utilises\n\
	-c --chemin	chaine de caractere : chemin vers le répertoire contenant le fichier fasta a analyser\n\
	\n\
	-h --help	affiche cette aide\n\
	\n");
}


/*********************/
/* liste des options */
/*********************/

struct option long_options[] = {
		{"longueurMotif", required_argument, NULL, 'l'},
		{"nbFenetre", required_argument, NULL, 'k'},
		{"chemin", required_argument, NULL, 'c'},
		{"help", no_argument, 0,  'h' }
};


///////////////////////////////////////////////////
/*Fonction permettant de récupérer les paramètres*/
///////////////////////////////////////////////////

void getParam (char **chemin, int* longueurMotif, int* nbFenetre, int argc, char *argv[]){
	FILE *verifFichier = NULL;
	int opt = 0;
	int long_index = 0; /* index des options */
	regex_t preg;
	regcomp(&preg,"\\w*\\S*\\.fa(sta)?",REG_EXTENDED|REG_NOSUB);
	if (argc<2)
	{
		fprintf(stderr, "ERREUR\n"); 
		notice(); 
		exit(1);
	}
	while ((opt = getopt_long(argc, argv, "l:k:c:h", long_options, &long_index )) != -1) {
		switch (opt) {
			case 'l' : *longueurMotif = atoi(optarg);
				break;
			case 'k' : *nbFenetre = atoi(optarg);
				break;
			case 'c' : 	*chemin=(char*)malloc(((strlen(optarg))+1)*sizeof(char));
						strcpy(*chemin, optarg);
				break;
			case 'h' : notice(); exit(0);
				break;
			default  : fprintf(stderr, "Mauvais argument\n"); notice(); exit(1);
				break;
		}
	}

	////////////////////////////////////
	/*VERIFICATION BON TYPE DE FICHIER*/
	////////////////////////////////////
	if (regexec(&preg,*chemin,0,0,0))
	{
		printf("[ERREUR] Le fichier doit être de type fasta (.fa ou .fasta)\n");
		notice();
		exit(1);
	}
	///////////////////////////////
	/*VERIFICATION FICHIER EXISTE*/
	///////////////////////////////

	verifFichier = fopen (*chemin, "r+");
	if (verifFichier == NULL)
	{
		printf("[ERREUR] Le chemin specifie est incorrect\n");
		notice();
		fclose(verifFichier);
		exit(1);
	}

	if (*longueurMotif == 0 || *nbFenetre == 0 || *nbFenetre > *longueurMotif)
	{
		printf("[ERREUR] La longueur du motif ne peut etre inferieure au nombre de fenetres.\n");
		notice();
		exit(1);
	}
	regfree(&preg);
	fclose(verifFichier);
}