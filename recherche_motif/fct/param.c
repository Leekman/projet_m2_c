#include "../lib/param.h"


//////////////////
/*AFFICHE L'AIDE*/
//////////////////
void notice(){
	printf("Utilisation : ./bin/recherche_motif -d entier -l entier -k entier -c chaine de caractere -o chaine de caractere -i entier\n\
Exemple : ./bin/recherche_motif -d 2 -l 4 -k 2 -c input/sequences.fasta -o output/resultats_sequences.txt -i 3\n\
OPTIONS :\n\
	-d --nbErreurMax	entier : Nombre maximum de substitutions\n\
	-l --longueurMotif	entier : Longueur du motif commun a identifier\n\
	-k --nbFenetre	entier : nombre de fenetres des masques utilises\n\
	-c --cheminEntree	chaine de caractere : chemin vers le fichier fasta a analyser\n\
	-o --output	chaine de caractere : chemin vers le fichier de résultat à écrire\n\
	-i --iterations entier : nombre de masques à créer\n\
	-h --help	affiche cette aide\n\
	\n");
}


/////////////////////
/*LISTE DES OPTIONS*/
/////////////////////

struct option long_options[] = {
		{"nbErreurMax", required_argument, NULL, 'd'},
		{"longueurMotif", required_argument, NULL, 'l'},
		{"nbFenetre", required_argument, NULL, 'k'},
		{"cheminEntree", required_argument, NULL, 'c'},
		{"output", required_argument, NULL, 'o'},
		{"iterations", required_argument, NULL, 'i'},
		{"help", no_argument, 0,  'h' }
};


///////////////////////////////////////////////////
/*FONCTION PERMETTANT DE RECUPERER LES PARAMETRES*/
///////////////////////////////////////////////////

void getParam (char **cheminEntree, char **output, int *nbIterations, int *nbErreurMax, int* longueurMotif, int* nbFenetre, int argc, char *argv[]){
	FILE *verifFichierEntree = NULL;
	FILE *verifFichierSortie = NULL; 
	int opt = 0;
	int long_index = 0; /* index des options */
	regex_t preg;
	regcomp(&preg,"\\w*\\S*\\.fa(sta)?",REG_EXTENDED|REG_NOSUB);

	//////////////////////////////////////////////////////////////////////////////////////////
	/*VERIFICATION DU BON NOMBRE D'ARGUMENTS /!\ LE NOM DU PROGRAMME COMPTE POUR UN ARGUMENT*/
	//////////////////////////////////////////////////////////////////////////////////////////
	
	
	if (argc < 2 || argc > 13)
	{
		fprintf(stderr, "\n\n[ERREUR] Nombre d'arguments incorrect\n\n"); 
		notice(); 
		regfree(&preg);
		exit(1);
	}

	//////////////////////////////////////////
	/*LECTURE ET RECUPERATION DES PARAMETRES*/
	//////////////////////////////////////////
	
	while ((opt = getopt_long(argc, argv, "d:l:k:c:o:i:h", long_options, &long_index )) != -1) {
		switch (opt) {
			case 'd' : *nbErreurMax = atoi(optarg);
						break;
			
			case 'l' : *longueurMotif = atoi(optarg);
						break;
			
			case 'k' : *nbFenetre = atoi(optarg);
						break;
			
			case 'c' : 	*cheminEntree=(char*)malloc(((strlen(optarg))+1)*sizeof(char));
						strcpy(*cheminEntree, optarg);
						break;
			
			case 'o' : *output=(char*)malloc(((strlen(optarg))+1)*sizeof(char));
						strcpy(*output, optarg);
						break;

			case 'i' : *nbIterations = atoi(optarg);
						break;
			
			case 'h' : 	notice(); exit(0);
						break;

			default  : 	printf("Mauvais argument\n"); notice(); exit(1);
						break;
		}
	}
	////////////////////////////////////
	/*VERIFICATION BON TYPE DE FICHIER*/
	////////////////////////////////////
	
	if (regexec(&preg,*cheminEntree,0,0,0))
	{
		printf("\n\n[ERREUR] Le fichier doit être de type fasta (.fa ou .fasta)\n\n");
		regfree(&preg);
		notice();
		exit(1);
	}
	//////////////////////////////////
	/*VERIFICATIONS FICHIERS EXISTES*/
	//////////////////////////////////
	
	verifFichierEntree = fopen (*cheminEntree, "r");
	if (verifFichierEntree == NULL)
	{
		printf("\n\n[ERREUR] Le chemin d'entrée spécifié est incorrect\n\n");
		notice();
		regfree(&preg);
		exit(1);
	}
	fclose(verifFichierEntree);
	verifFichierSortie = fopen (*output, "a");
	if (verifFichierSortie == NULL)
	{
		printf("\n\n[ERREUR] Le chemin de sortie spécifié est incorrect\n\n");
		notice();
		regfree(&preg);
		exit(1);
	}
	fclose(verifFichierSortie);
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*VERIFICATION DE LA RECUPERATION DES PARAMETRES AINSI QUE DE LA CONDITION NOMBRE FENETRE > LONGUEUR DU MOTIF*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	if (*longueurMotif == 0 || *nbFenetre == 0 || *nbIterations == 0|| *nbFenetre > *longueurMotif)
	{
		printf("\n\n[ERREUR] Erreur lors de la saisie du nombre de fenetre, de la longueur du motif ou du nombre d'iterations.\n\n");
		notice();
		regfree(&preg);
		exit(1);
	}

	if (*nbFenetre > *longueurMotif)
	{
		printf("\n\n[ERREUR] Le nombre de fenetres ne peut être supérieur a la longueur du motif ! \n");
		notice();
		regfree(&preg);
		exit(1);
	}

	regfree(&preg);
}