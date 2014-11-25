#include "../lib/masque.h"


int compare (void const *a, void const *b){

   int const *pa = a;
   int const *pb = b;

   return *pa - *pb;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool estDoublon(int n, int *doublon, int count){

	int i;

	for (i=0; i<count; i++) 
	{
		if (n==doublon[i])
		{
			return true;
		}
	}
	return false;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int *generateurMasque(int l, int k){
	
	int *masque;
	int i;

	masque = (int*) malloc((l-k)*sizeof(int));

	for (i=0; i<(l-k); i++)
	{
		do
		{
		masque[i]=rand()%l;	
		} while (estDoublon(masque[i], masque, i));		
	}

	qsort(masque, (l-k), sizeof (int), compare);

	return masque;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Test du masque

/*int main(){
	int *masque;
	int l,k;
	int i;
	l=10;
	k=4;

	masque=generateurMasque(l,k);


	for (i=0;i<(l-k);i++)
		printf("masque[%d] = %d\n",i,masque[i]);


	free (masque);
	return 0;
}*/