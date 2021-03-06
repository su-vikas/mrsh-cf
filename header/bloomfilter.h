/* 
 * File:   bloom.h
 * Author: mustafakarabat
 *
 * Created on 17. April 2012, 23:34
 */

#ifndef BLOOMFILTER_H
#define	BLOOMFILTER_H

/* 
 * We define a struct BLOOM, with all the properties a BLoom-Filter needs.
 */
typedef struct {
    // For a filter_size of 256 Bytes.
    unsigned char array[FILTERSIZE];
    
    // We store the number of blocks we add to each filter in count_added_blocks
    short int amount_of_blocks;
    
    // Pointer to next Bloomfilter
    struct BLOOMFILTER *next;
    
}BLOOMFILTER;



BLOOMFILTER     *init_empty_BF();
void            destroy_bf(BLOOMFILTER *bf);

void            bloom_set_bit(unsigned char *bit_array, unsigned short value);
unsigned short  count_bits_set_to_one_of_BF(uint64 *filter);
unsigned short  bloom_common_bits(uint64 *bit_array_one, uint64 *bit_array_two);

void            add_hash_to_bloomfilter(BLOOMFILTER *bf, uint64 hash_value);
void            convert_hex_binary(const unsigned char *hex_string, BLOOMFILTER *bf);


#endif	/* BLOOM_H */




