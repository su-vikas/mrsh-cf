/* 
 * File:   config.h
 * Author: Frank Breitinger
 *
 * Created on 1. Mai 2012, 12:15
 */

#ifndef CONFIG_H
#define	CONFIG_H

#define ROLLING_WINDOW          7
#define BLOCK_SIZE              320// 256
#define SHIFTOPS                11
#define MASK                    0x7FF
#define FILTERSIZE              256
#define SUBHASHES               5
#define BLOOMFILTERBITSIZE      (FILTERSIZE * 8)
#define MAXBLOCKS               192 //128
#define MINBLOCKS				10 //if a Bloom filter has less than MINBLOCKS it is skipped
#define SKIPPED_BYTES           BLOCK_SIZE/4
#define PROBABILITY             0.99951172 //Attention: 1 - ( 1/BLOOMFILTERBITSIZE )

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)

typedef unsigned long long  uint64; 
typedef unsigned char       uchar;
typedef unsigned int        uint32;

/*typedef short bool;*/
//#define true 1
//#define false 0

typedef struct{
    bool compare;
    bool file_comparison;
    bool print;
    bool gen_compare;
    bool compareLists;
    bool helpmessage;
    bool recursive;
    bool path_list_compare;
    short threshold;
} MODES;


extern MODES *mode; //= {.compare = false}


#endif	/* CONFIG_H */

