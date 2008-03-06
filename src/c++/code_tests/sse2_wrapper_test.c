
#include "sse2_wrapper.h"
#include <stdio.h>
#include <assert.h>


int
main(int argc,
  char **argv){

  b128 a,b;
  b128 *ptr;
  int i;
  
  printf("Creating zero\n");
  print_ints_b128(set_zero_b128());

  printf("Setting first to 2\n");
  print_ints_b128(set_first_int_b128(2));

  printf("Setting all 3 2 1 0\n");
  print_ints_b128(set_ints_b128(3,2,1,0));


  printf("Adding {3 2 1 0}+{3 3 3 3}\n");
  print_ints_b128(add_b128(set_ints_b128(3,2,1,0),set_ints_b128(3,3,3,3)));
  assert(equal_b128(add_b128(set_ints_b128(3,2,1,0),set_ints_b128(3,3,3,3)),set_ints_b128(6,5,4,3)));
  
  printf("Subtracting {3 2 1 0}-{3 3 3 3}\n");
  print_ints_b128(sub_b128(set_ints_b128(3,2,1,0),set_ints_b128(3,3,3,3)));
  assert(equal_b128(sub_b128(set_ints_b128(3,2,1,0),set_ints_b128(3,3,3,3)),set_ints_b128(0,-1,-2,-3)));

  
  printf("Immediate access\n");
  a = set_ints_b128(3,2,1,0);
  printf("%d  =   %d\n",get_int_0_b128(a),get_int_b128(a,0));
  printf("%d  =   %d\n",get_int_1_b128(a),get_int_b128(a,1));
  printf("%d  =   %d\n",get_int_2_b128(a),get_int_b128(a,2));
  printf("%d  =   %d\n",get_int_3_b128(a),get_int_b128(a,3));
  i = get_immediate_int16_b128(a,6);
  printf("%d  =   %d\n",i,get_int_b128(a,3));

  a = set_int_0_b128(a,3);
  a = set_int_1_b128(a,2);
  a = set_int_2_b128(a,1);
  a = set_int_3_b128(a,0);
  assert(equal_b128(a,set_ints_b128(0,1,2,3)));
  

  
  printf("Setting bits 1 50 100 127\n");
  a = set_zero_b128();
  a = set_bit_b128(a,1,1);
  a = set_bit_b128(a,50,1);
  a = set_bit_b128(a,100,1);
  a = set_bit_b128(a,127,1);
  print_bits_b128(a);
  assert(get_bit_b128(a,1));
  assert(get_bit_b128(a,50));
  assert(get_bit_b128(a,100));
  assert(get_bit_b128(a,127));


  printf(" xor \n");
  a = set_zero_b128();
  a = set_bit_b128(a,1,1);
  a = set_bit_b128(a,50,1);
  a = set_bit_b128(a,100,1);
  a = set_bit_b128(a,127,1);
  b = set_zero_b128();
  b = set_bit_b128(b,1,1);
  b = set_bit_b128(b,49,1);
  b = set_bit_b128(b,102,1);
  b = set_bit_b128(b,127,1);
  print_bits_b128(a);
  print_bits_b128(b);
  print_bits_b128(xor_b128(a,b));
  assert(get_bit_b128(xor_b128(a,b),127)==0);
  assert(get_bit_b128(xor_b128(a,b),102)==1);

  printf(" and \n");
  a = set_zero_b128();
  a = set_bit_b128(a,1,1);
  a = set_bit_b128(a,32,1);
  a = set_bit_b128(a,100,1);
  a = set_bit_b128(a,127,1);
  b = set_zero_b128();
  b = set_bit_b128(b,1,1);
  b = set_bit_b128(b,49,1);
  b = set_bit_b128(b,102,1);
  b = set_bit_b128(b,127,1);
  print_bits_b128(a);
  print_bits_b128(b);
  print_bits_b128(and_b128(a,b));
  assert(get_bit_b128(and_b128(a,b),127)==1);
  assert(get_bit_b128(and_b128(a,b),102)==0);

  printf(" or \n");
  a = set_zero_b128();
  a = set_bit_b128(a,1,1);
  a = set_bit_b128(a,32,1);
  a = set_bit_b128(a,100,1);
  a = set_bit_b128(a,127,1);
  b = set_zero_b128();
  b = set_bit_b128(b,1,1);
  b = set_bit_b128(b,31,1);
  b = set_bit_b128(b,102,1);
  b = set_bit_b128(b,127,1);
  print_bits_b128(a);
  print_bits_b128(b);
  print_bits_b128(or_b128(a,b));
  assert(get_bit_b128(or_b128(a,b),127)==1);
  assert(get_bit_b128(or_b128(a,b),102)==1);

  printf(" negate \n");
  a = negate_b128(set_ints_b128(0x12345678,0x12345678,0x12345678,0x12345678));
  b = set_ints_b128(~0x12345678,~0x12345678,~0x12345678,~0x12345678);
  print_bits_b128(a);
  print_bits_b128(b);
  assert(equal_b128(a,b));


  printf(" SHIFT BYTES LEFT\n");
  a = set_zero_b128();
  a = set_bit_b128(a,1,1);
  a = set_bit_b128(a,32,1);
  a = set_bit_b128(a,100,1);
  a = set_bit_b128(a,127,1);
  b = shift_bytes_left_b128(a,1);
  print_bits_b128(a);
  print_bits_b128(b);
  assert(get_bit_b128(b,9) == 1);
  assert(get_bit_b128(b,40) == 1);
  assert(get_bit_b128(b,127) == 0);
  printf(" SHIFT BYTES RIGHT\n");
  b = shift_bytes_right_b128(a,2);
  print_bits_b128(a);
  print_bits_b128(b);
  assert(get_bit_b128(b,16) == 1);
  
  printf(" SHIFT EACH BITS RIGHT\n");
  a = set_zero_b128();
  a = set_bit_b128(a,3,1);
  a = set_bit_b128(a,32,1);
  a = set_bit_b128(a,63,1);
  a = set_bit_b128(a,100,1);
  a = set_bit_b128(a,127,1);
  b = shift_each32_bits_right_b128(a,set_first_int_b128(3));
  print_bits_b128(a);
  print_bits_b128(b);
  assert(get_bit_b128(b,0) == 1);
  assert(get_bit_b128(b,63) == 0);


  printf("Shift 16 bits\n");
  a = set_ints_b128(0x00020003,0x00010002,0 ,0x00010001 );
  print_blocks_b128(a,16);
  b = shift_each32_bits_right_b128(a,set_first_int_b128(16));
  print_blocks_b128(b,16);
  print_blocks_b128(b,32);
  
  printf("MEM CHECKS\n");
  ptr = alloc_b128(10);
  set_bit_in_mem_b128(ptr,129,1);
  assert(get_bit_in_mem_b128(ptr,129) == 1);
  print_bits_b128(get_from_mem_b128(ptr,1));
  print_bits_b128(ptr[1]);
  
  free_b128(ptr);

  printf("Calloc CHECKS\n");
  ptr = calloc_b128(100000);
  print_bits_b128(ptr[0]);
  for (i = 0 ; i < 100000 ; i++ ){
    assert( equal_b128(ptr[i], set_zero_b128()) );
    if ( !equal_b128(ptr[i], set_zero_b128() )){
      printf("Error: not equal calloc\n");
    }
  }
  free_b128(ptr);

  
  printf("Block printing\n");
  a = set_ints_b128(0x1234,0x4321,0xf0f,0xf0f0f0f0);
  print_blocks_b128(a,4);


  return 1;
}
