#include <stdio.h>
#include <stdint.h>
#include "f1v1dec.h"

void
f1v1DataDecode(uint32_t data)
{
  f1_data_t f1word;

  f1word.raw = data;

  printf("%8X - Slot %2u Err %d%d%d",
	 f1word.raw,
	 f1word.bf.slot_id,
	 f1word.bf.res_locked,
	 f1word.bf.output_fifo_overflow,
	 f1word.bf.hit_fifo_overflow);

  switch (f1word.bf.chip_word_type)
    {
    case 0:			/* CHIP HEADER / TRAILER */
      {
	f1_chip_header_t d;
	d.raw = data;
	printf(" - H/T Chip=%u Chan=%u Err=%d XOR=%d  Trignum=%2d  Time=%3u\n",
	       d.bf.channel_address,	// chip
	       d.bf.chip_address,	// chan
	       d.bf.trigger_fifo_overflow,	// ERR
	       d.bf.Xor_setup_register, // XOR
	       d.bf.trigger_number,     // TRIGNUM
	       d.bf.trigger_time);      // TIME

	break;
      }

    case 1:			/* CHIP DATA */
      {
    	f1_chip_data_t d;
    	d.raw = data;
    	printf(" - DATA   Chip=%u   Chan=%u   Time=%4u\n",
    	       d.bf.chip_address,	// chip
    	       d.bf.channel_address,	// chan
    	       d.bf.time_measurement);	// t

    	break;
      }

    default:
      {
	printf(" - UNDEFINED TYPE = %d\n",
	       f1word.bf.chip_word_type);
	break;
      }
    }
}
