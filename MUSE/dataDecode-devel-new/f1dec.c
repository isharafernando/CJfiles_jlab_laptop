#include <stdio.h>
#include <stdint.h>
#include "jlabdec.h"
#include "f1dec.h"

void
f1DataDecode(uint32_t data)
{
  static uint32_t type_last = 15;	/* initialize to type FILLER WORD */
  static uint32_t time_last = 0;
  static int new_type = 0;
  int type_current = 0;
  generic_data_word_t gword;

  gword.raw = data;

  if (gword.bf.data_type_defining)	/* data type defining word */
    {
      new_type = 1;
      type_current = gword.bf.data_type_tag;
    }
  else
    {
      new_type = 0;
      type_current = type_last;
    }

  switch (type_current)
    {
    case 0:			/* BLOCK HEADER */
      {
	block_header_t d;
	d.raw = data;

	printf
	  ("%8X - BLOCK HEADER - slot = %d  modID = %d   n_evts = %d   n_blk = %d\n",
	   d.raw, d.bf.slot_number, d.bf.module_ID,
	   d.bf.number_of_events_in_block, d.bf.event_block_number);

	break;
      }

    case 1:			/* BLOCK TRAILER */
      {
	f1_block_trailer_t d;
	d.raw = data;

	printf("%8X - BLOCK TRAILER - slot = %d   n_words = %d\n",
	       d.raw, d.bf.slot_number, d.bf.words_in_block);

	break;
      }

    case 2:			/* EVENT HEADER */
      {
	event_header_t d;
	d.raw = data;

	printf("%8X - EVENT HEADER - trig num = %d\n",
	       d.raw, d.bf.event_number);
	break;


	break;
      }

    case 3:			/* TRIGGER TIME */
      {
	if (new_type)
	  {
	    f1_trigger_time_1_t d;
	    d.raw = data;

	    printf("%8X - TRIGGER TIME 1 - time = %06x\n",
		   d.raw, (d.bf.T_C << 16) | (d.bf.T_D << 8) | (d.bf.T_E));

	    time_last = 1;
	  }
	else
	  {
	    if (time_last == 1)
	      {
		f1_trigger_time_2_t d;
		d.raw = data;
		printf("%8X - TRIGGER TIME 2 - time = %04x\n",
		       d.raw, (d.bf.T_A << 8) | (d.bf.T_B));

	      }
	    else
	      printf("%8X - TRIGGER TIME - (ERROR)\n", data);
	  }
	break;
      }

    case 7:			/* EVENT DATA */
      {
	f1_time_measurement_t d;
	d.raw = data;
	printf("%8X - TIME MEASUREMENT   ERR=%X  fake = %u  chip=%u  chan=%u  t = %u\n",
	       d.raw,
	       (d.bf.resolution_locked << 2) |
	       (d.bf.output_fifo_overflow << 1) |
	       d.bf.hit_fifo_overflow,	// ERR
	       d.bf.fake,	// Fake data flag
	       d.bf.chip_number,	// chip
	       d.bf.chip_channel_number,	// chan
	       d.bf.time_measurement);	// t

	break;
      }

    case 8:			/* CHIP HEADER */
      {
	f1_chip_header_t d;
	d.raw = data;
	printf("%8X - CHIP HEADER (%X,%u)  chip=%u  chan=%u  trig = %u  t = %3u\n",
	       d.raw,
	       (d.bf.resolution_locked << 2) |
	       (d.bf.output_fifo_overflow << 1) |
	       d.bf.hit_fifo_overflow,	// ERR
	       d.bf.setup_register_tag, d.bf.chip_number,	// chip
	       d.bf.chip_channel_number,	// chan
	       d.bf.trigger_number, d.bf.trigger_time);

	break;
      }

    case 14:			/* DATA NOT VALID (no data available) */
      {
	data_not_valid_t d;
	d.raw = data;

	printf("%8X - DATA NOT VALID = %d\n", d.raw, d.bf.data_type_tag);

	break;
      }

    case 15:			/* FILLER WORD */
      {
	filler_word_t d;
	d.raw = data;

	printf("%8X - FILLER WORD = %d\n", d.raw, d.bf.data_type_tag);

	break;
      }

    case 4:			/* UNDEFINED TYPE */
    case 5:			/* UNDEFINED TYPE */
    case 6:			/* UNDEFINED TYPE */
    case 9:			/* UNDEFINED TYPE */
    case 10:			/* UNDEFINED TYPE */
    case 11:			/* UNDEFINED TYPE */
    case 12:			/* UNDEFINED TYPE */
    case 13:			/* UNDEFINED TYPE */
    default:
      {
	printf("%8X - UNDEFINED TYPE = %d\n",
	       gword.raw, gword.bf.data_type_tag);
	break;
      }
    }

  type_last = type_current;	/* save type of current data word */
}
