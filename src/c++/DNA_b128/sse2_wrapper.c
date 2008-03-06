//--------------------------------------------------
//                                        
// File: sse2_wrapper.c
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: sse2_wrapper.c,v 1.2 2006/11/27 14:20:18 isaac Exp $                                 
//
//--------------------------------------------------


#include <stdio.h>
#include <stdlib.h>
#include "sse2_wrapper.h"



//---------------------------
// MEMORY ALLOC/FREE
// Allocates num_b128*16 bytes memory that is 16 bytes aligned.
//
// Returns a pointer to memory that is aligned according to align_size.
// The memory looks as follows:
//
//  X                       Y
//   ---------------------------------------------
//  |            |address X | start aligned memory
//   ---------------------------------------------
//  
// The return pointer is Y which is an address on the form 0x?????????0.
// (i.e. the last hex digit is zero)
//
// To be able to free the memory the pointer to X is saved in the int just before
// the aligned memory.
b128 *
alloc_b128(int num_b128) {

  char *ptr,*ptr2,*aligned_ptr;

  ptr=(char *)malloc(num_b128*sizeof(b128) + 16 + sizeof(int));
  if( ptr == NULL ) return NULL;
  //printf("ptr %x\n",ptr);
  
  //Compute the aligned pointer and make sure that there is
  //enough room for the pointer to the beginning of the allocated block.
  ptr2 = ptr + sizeof(int);
  aligned_ptr = ptr2 + (16 - ((size_t)ptr2 & 0xF));
  //printf("align_ptr %x\n",aligned_ptr);
 
  //Save the pointer to the beginning of the allocated block.
  ptr2 = aligned_ptr - sizeof(int);
  *((int *)ptr2) = (int)(aligned_ptr - ptr);
  //printf("ptr2pos %x\n",(aligned_ptr - ptr));
  
  return (b128 *) (aligned_ptr);
}
//--------
// CALLOC (exactly the same as alloc_b128 except for the calloc part.
b128 *
calloc_b128(int num_b128) {

  char *ptr,*ptr2,*aligned_ptr;

  ptr=(char *)calloc(num_b128*sizeof(b128) + 16 + sizeof(int),sizeof(char));
  if( ptr == NULL ) return NULL;
  //printf("ptr %x\n",ptr);
  
  //Compute the aligned pointer and make sure that there is
  //enough room for the pointer to the beginning of the allocated block.
  ptr2 = ptr + sizeof(int);
  aligned_ptr = ptr2 + (16 - ((size_t)ptr2 & 0xF));
  //printf("align_ptr %x\n",aligned_ptr);
 
  //Save the pointer to the beginning of the allocated block.
  ptr2 = aligned_ptr - sizeof(int);
  *((int *)ptr2) = (int)(aligned_ptr - ptr);
  //printf("ptr2pos %x\n",(aligned_ptr - ptr));
  
  return (b128 *) (aligned_ptr);
}


void 
free_b128(b128 *_ptr) {
  char *ptr = (char *)_ptr; 
  //printf("free ptr %x\n",ptr);
  //Collect the pointer to the beginning of the allocated block.
  int *ptr2 = ((int *)ptr) - 1;
  //printf("ptr2pos %x\n",*ptr2);
  ptr -= *ptr2;
  //printf("free %x\n",ptr);
  free(ptr);
}





//---------------------------
// PRINTING

void
print_bits_32(int a){
  int i;
  for (  i = 31 ; i > -1 ; i-- ){
    printf("%d",((a>>i)&0x1));
    if ( i%8 == 0 )
      printf(" ");
  }
}

void
print_bits_b128(b128 a){

  print_bits_32(get_int_3_b128(a));
  printf(",");
  print_bits_32(get_int_2_b128(a));
  printf(",");
  print_bits_32(get_int_1_b128(a));
  printf(",");
  print_bits_32(get_int_0_b128(a));
  printf(",");
  

  printf("\n");
}


void
print_ints_b128(b128 a){

  printf("%d %d %d %d\n",get_int_3_b128(a),get_int_2_b128(a),get_int_1_b128(a),get_int_0_b128(a));
}



//Divides a into blocks of size block_size and prints the integer in
//each such block. block_size has to be a power of 2.
void
print_blocks_b128(b128 a, int block_size){
  int working_int;
  int shift, b_int;
  //have to special handle 32 since we can not shift an int 32 bits
  int block_mask = (block_size == 32 ? -1 : (1<<block_size)-1);
  int blocks_per_int = 32/block_size;
  
  for ( working_int = 3 ; working_int > -1 ; working_int--){
    b_int = get_int_b128(a,working_int);
    
    for ( shift = block_size*(blocks_per_int-1) ; shift > -1 ; shift-=block_size ){
      printf("%d ", ((b_int>>shift) & block_mask));
    }
    printf(",");
  }

  printf("\n");

}
