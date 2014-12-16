#include "../lib/includes.h"
#include "../lib/param.h"

//////////////////
/*AFFICHE LA NOTICE UTILISATEUR*/
//////////////////
void notice(){
	printf("Utilisation : creation_sequence -d entier -n entier -t entier -m chaine de caracteres -v entier -o chaine de caracteres\n\
Description :\n\
Exemple : ./bin/generation_sequence -d 2 -n 50 -t 100 -m TATATGCA -v 50 -o ../recherche_motif/input/sequences.fasta -i output/info.txt \n\
OPTIONS :\n\
	-d --nbErreur	entier : nombre de substitutions maximum autorisees\n\
	-n --nbSeq	entier : nombre de sequences a creer\n\
	-t --tailleSeq	entier : taille des sequences a creer\n\
	-m --motif	chaine de caracteres : motif a inserer dans la sequence\n\
	-v --variable entier : la taille des séquences est variable de x %%\n\
	-o --output chaine de caracteres : chemin de sortie du fichier fasta\n\
	-i --info chaine de caracteres : chemin de sortie du fichier d'informations\n\
	-h --help	affiche cette aide\n\
	\n");
}


/////////////////////
/*LISTE DES OPTIONS*/
/////////////////////

struct option long_options[] = {
		{"nbErreur", required_argument, NULL, 'd'},
		{"nbSeq", required_argument, NULL, 'n'},
		{"tailleSeq", required_argument, NULL, 't'},
		{"motif", required_argument, NULL, 'm'},
		{"variable", required_argument, NULL, 'v'},
		{"output", required_argument, NULL, 'o'},
		{"info", required_argument, NULL, 'i'},
		{"help", no_argument, 0,  'h' }
};

///////////////////////////////////////////////////
/*Fonction permettant de récupérer les paramètres*/
///////////////////////////////////////////////////

void getParam (int *nbErreur, int *nbSeq, int *tailleSeq, char **motif, char **cheminFasta, char** cheminInfo, int *variation, int argc, char *argv[]){
	
	FILE *verifSortieFasta = NULL;
	FILE *verifSortieInfo = NULL;
	int opt = 0;
	int long_index = 0; /* index des options */
	////////////////////////////////////////////////
	/*VERIFICATION QU'ON A AU MOINS DEUX ARGUMENTS*/
	////////////////////////////////////////////////
	
	if (argc<2)
	{
		fprintf(stderr, "ERREUR\n"); 
		notice(); 
		exit(1);
	}
	regex_t preg;
	regcomp(&preg,"[^ATCG]+",REG_EXTENDED|REG_NOSUB);//expression régulière pour ATCG uniquement pour le motif
	                                       
	while ((opt = getopt_long(argc, argv, "d:n:t:m:v:o:i:h", long_options, &long_index )) != -1) {
		switch (opt) {
			
			case 'd' : 	*nbErreur = atoi(optarg);
						break;
			
			case 'n' : 	*nbSeq = atoi(optarg);
						break;
			
			case 't' : 	*tailleSeq = atoi(optarg);
						break;
			
			case 'm' : 	*motif=(char*)malloc(((strlen(optarg))+1)*sizeof(char));
						strcpy(*motif, optarg);
						break;

			case 'v' : 	*variation = atoi(optarg);
						break;

			case 'o' : *cheminFasta = (char*)malloc(((strlen(optarg))+1)*sizeof(char));
						strcpy(*cheminFasta, optarg);
						break;

			case 'i' : *cheminInfo = (char*)malloc(((strlen(optarg))+1)*sizeof(char));
						strcpy(*cheminInfo, optarg);
						break;

			case 'h' : 	notice(); exit(0);
						break;
			
			default  : 	fprintf(stderr, "Mauvais argument\n"); notice(); exit(1);
						break;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	/*VERIFICATION QUE LES CHEMINS DE SORTIES SONT CORRECT ET QUE LES FICHIERS PEUVENT S'OUVRIR*/
	/////////////////////////////////////////////////////////////////////////////////////////////

	verifSortieFasta = fopen(*cheminFasta, "w");
	if (verifSortieFasta == NULL)
	{
		printf("\n\n[ERREUR] Le chemin pour le fichier fasta est incorrect\n\n");
		notice();
		regfree(&preg);
		exit(1);
	}
	fclose(verifSortieFasta);


	verifSortieInfo = fopen(*cheminInfo, "w");
	if (verifSortieInfo == NULL)
	{
		printf("\n\n[ERREUR] Le chemin pour le fichier d'informations est incorrect\n\n");
		notice();
		regfree(&preg);
		exit(1);
	}
	fclose(verifSortieInfo);

	///////////////////////////////////////////////////////////////////////////////////////
	/*VERIFICATION QUE LA TAILLE DE LA SEQUENCE AINSI QUE LE NOMBRE ONT BIEN ETE SPECIFIE*/
	///////////////////////////////////////////////////////////////////////////////////////

	if (*tailleSeq == 0 || *nbSeq < 2)   
    {
        printf("\nERREUR Saisissez une taille ou un nombre de sequences correct\n");
        notice();
        regfree(&preg);
        exit(1);
    }

    ////////////////////////////////////////////////////////////
    /*TAILLE DU MOTIF NE PEUT EXCEDER LA TAILLE DE LA SEQUENCE*/
    ////////////////////////////////////////////////////////////

    if (strlen(*motif) > *tailleSeq)
    {
    	printf("ERREUR La taille du motif ne peut excéder la taille de la séquence !");
    	notice();
    	regfree(&preg);
    	exit(1);
    }
    ///////////////////////////////////////
    /*VERIFICATION QU'ON EST BIEN EN ATCG*/
    ///////////////////////////////////////

	if (!regexec(&preg,*motif,0,0,0))
	{
		printf("[ERREUR] Le motif ne peux contenir que les caracteres A,T,C ou G\n");
		notice();
		regfree(&preg);
		exit(1);
	}
	regfree(&preg);
}

