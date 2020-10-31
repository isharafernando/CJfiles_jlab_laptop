#include <stdio.h>
#include <stdint.h>
#include "caen775dec.h"

void 
caen775DataDecode(uint32_t data)
{
  caen775_generic_word_t gword;

  gword.raw = data;
  
  switch(gword.bf.word_type)
    {
    case CAEN775_DATATYPE_HEADER:
      {
	caen775_header_t d; d.raw = data;

	printf("%8X - HEADER - slot = %2d  crate = %3d  data_count = %d\n",
	       d.raw,
	       d.bf.geo,
	       d.bf.crate,
	       d.bf.data_count);
	break;
      }

    case CAEN775_DATATYPE_DATA:
      {
	caen775_data_t d; d.raw = data;

	printf("%8X - DATA - slot = %2d  channel = %2d V = %d UN = %d OV = %d  ADC = %d\n",
	       d.raw,
	       d.bf.geo,
	       d.bf.channel_number,
	       d.bf.valid,
	       d.bf.underflow,
	       d.bf.overflow,
	       d.bf.adc);
	break;
      }
      
    case CAEN775_DATATYPE_END_OF_BLOCK:
      {
	caen775_end_of_block_t d; d.raw = data;

	printf("%8X - EOB - slot = %2d  event counter = %d\n",
	       d.raw,
	       d.bf.geo,
	       d.bf.event_counter);
	break;
      }
      
    default:
      {
	printf("%8X - UNDEFINED TYPE = %d\n",
	       gword.raw,
	       gword.bf.word_type);
      }


    }
}
