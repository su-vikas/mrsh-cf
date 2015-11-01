#ifndef timing_h
#define timing_h

#include <time.h>
//MACRO FOR PRINTING TIME TO EXECUTE PART OF CODE
clock_t startm, stopm;
#define START if ( (startm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define STOP if ( (stopm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define PRINTTIME printf( "[*] time:%6.3f \n", ((double)stopm-startm)/CLOCKS_PER_SEC);
#define TIME ((double)stopm-startm)/CLOCKS_PER_SEC


#endif


