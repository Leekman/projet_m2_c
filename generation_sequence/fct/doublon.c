#include "../lib/doublon.h"
#include "../lib/includes.h"

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