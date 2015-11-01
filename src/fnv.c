#include "../header/fnv.h"
#include <stdio.h>
#include <string.h>




uint64 fnv64Bit( unsigned char pBuffer[], int start, int end)
 {
   uint64 nHashVal    = 0xcbf29ce484222325ULL,
          nMagicPrime = 0x00000100000001b3ULL;

   int i = start;
   while( i <= end ) {
	   nHashVal ^= pBuffer[i++];
	   nHashVal *= nMagicPrime;
   }
   return nHashVal;
 }


/*
void fnv64Bit(char hashstring[], uint64 *hashv, int start, int end)
{
    //unsigned char *s = (unsigned char *)hashstring;
    while (start < end) {
        *hashv ^= (unsigned char)hashstring[start++];
        *(hashv) += (*(hashv) << 1) + (*(hashv) << 4) + (*(hashv) << 5) +
                      (*(hashv) << 7) + (*(hashv) << 8) + (*(hashv) << 40);
    }
}*/







