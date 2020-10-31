#include <stdio.h>
#include <stdint.h>
#include "fa125dec.h"
#include "jlabdec.h"


void
fa125DataDecode(uint32_t data)
{
  static uint32_t type_last = 15;	/* initialize to type FILLER WORD */
  static uint32_t time_last = 0;
  static int new_type = 0;
  int type_current = 0;
  static int isample = 0;
  static int ipk = 0;
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

  switch( type_current )
    {
    case 0:			/* BLOCK HEADER */
      {
	block_header_t d; d.raw = data;
	
	printf("%8X - BLOCK HEADER - slot = %d  modID = %d   n_evts = %d   n_blk = %d\n",
	       d.raw,
	       d.bf.slot_number,
	       d.bf.module_ID,
	       d.bf.number_of_events_in_block,
	       d.bf.event_block_number);
	
	break;
      }

    case 1:			/* BLOCK TRAILER */
      {
	block_trailer_t d; d.raw = data;
	
	printf("%8X - BLOCK TRAILER - slot = %d   n_words = %d\n",
		     d.raw,
		     d.bf.slot_number,
		     d.bf.words_in_block);
	break;
      }

    case 2:			/* EVENT HEADER */
      {
	event_header_t d; d.raw = data;

	printf("%8X - EVENT HEADER - slot = %d  evt_num = %d\n",
	       d.raw, 
	       d.bf.slot_number,
	       d.bf.event_number);
	break;
      }

    case 3:			/* TRIGGER TIME */
      {
	if (new_type)
	  {
	    fa125_trigger_time_1_t d; d.raw = data;
	    printf("%8X - TRIGGER TIME 1 - time = 0x%06x\n",
		   d.raw,
		   (d.bf.T_A << 16) | (d.bf.T_B << 8) | d.bf.T_C);

	    time_last = 1;
	  }
	else
	  {
	    if (time_last == 1)
	      {
		fa125_trigger_time_2_t d; d.raw = data;
		printf("%8X - TRIGGER TIME 2 - time = 0x%06x\n",
		       d.raw,
		       (d.bf.T_D << 16) | (d.bf.T_E << 8) | d.bf.T_F);
	      }
	    else
	      printf("%8X - TRIGGER TIME - (ERROR)\n",
		     data);

	    time_last = 0;
	  }
	break;
      }

    case 4:			/* WINDOW RAW DATA */
      {
	if (new_type)
	  {
	    fa125_window_raw_data_1_t d; d.raw = data;

	    printf("%8X - WINDOW RAW DATA - chan = %2d   nsamples = %d\n", 
		   d.raw,
		   d.bf.channel_number,
		   d.bf.PTW);

	    isample = 0;
	  }
	else
	  {
	    fa125_window_raw_data_n_t d; d.raw = data;

	    printf("%8X - RAW SAMPLES (%3d) - valid = %d  adc = %4d   valid = %d  adc = %4d\n", 
		   d.raw,
		   isample,
		   !d.bf.invalid_1,
		   d.bf.adc_sample_1,
		   !d.bf.invalid_2,
		   d.bf.adc_sample_2);		   


	    isample += 2;
	  }
	break;
      }

    case 5:			/* PULSE DATA, CDC */
      {
	if (new_type)
	  {
	    fa125_pulse_cdc_1_t d; d.raw = data;
	    printf("%8X - PULSE DATA (CDC IT) - chan = %2d  LE time = %d Q = %d OVF = %d\n",
		   d.raw,
		   d.bf.channel_number,
		   d.bf.leading_edge_time,
		   d.bf.time_quality,
		   d.bf.overflow_count);

	  }
	else
	  {
	    fa125_pulse_cdc_2_t d; d.raw = data;
	    printf("%8X - PULSE DATA (CDC IT) - ped = %d  integral = %d  firstmax ampl = %d\n",
		   d.raw,
		   d.bf.pedestal,
		   d.bf.integral,
		   d.bf.first_max_amplitude);

	  }
	break;
      }

    case 6:			/* PULSE DATA, FDC - Integral and Time */
      {
	if (new_type)
	  {
	    ipk = 0;
	    fa125_pulse_fdc_integral_1_t d; d.raw = data;
	    printf("%8X - PULSE DATA (FDC IT) - chan = %2d  NPK = %d  LE time = %d Q = %d OVF = %d\n",
		   d.raw,
		   d.bf.channel_number,
		   d.bf.NPK,
		   d.bf.leading_edge_time,
		   d.bf.time_quality,
		   d.bf.overflow_count);
	  }
	else
	  {
	    ipk++;
	    fa125_pulse_fdc_integral_n_t d; d.raw = data;
	    printf("%8X - PULSE DATA (FDC IT) %d - ped = %d  integral = %d  peak time = %d\n",
		   d.raw,
		   ipk,
		   d.bf.pedestal,
		   d.bf.integral,
		   d.bf.peak_time);
	  }
	break;
      }

    case 9:			/* PULSE DATA, FDC - Peak Ampl and Time */
      {
	if (new_type)
	  {
	    ipk = 0;
	    fa125_pulse_fdc_amplitude_1_t d; d.raw = data;
	    printf("%8X - PULSE DATA (FDC AT) - chan = %2d  NPK = %d  LE time = %d Q = %d OVF = %d\n",
		   d.raw,
		   d.bf.channel_number,
		   d.bf.NPK,
		   d.bf.leading_edge_time,
		   d.bf.time_quality,
		   d.bf.overflow_count);
	  }
	else
	  {
	    ipk++;
	    fa125_pulse_fdc_amplitude_n_t d; d.raw = data;
	    printf("%8X - PULSE DATA (FDC AT) %d - Ampl = %d  Time = %d  Pedestal = %d\n",
		   d.raw,
		   ipk,
		   d.bf.amplitude,
		   d.bf.peak_time,
		   d.bf.pedestal);
	  }
	break;
      }

    case 13:			/* END OF EVENT */
      {
	event_trailer_t d;
	d.raw = data;

	printf("%8X - END OF EVENT = %d\n", d.raw, d.bf.data_type_tag);
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

    case 7:
    case 8:
    case 10:
    case 11:
    case 12:			/* UNDEFINED TYPE */
    default:
      {
	printf("%8X - UNDEFINED TYPE = %d\n",
	       gword.raw, gword.bf.data_type_tag);
	break;
      }

    }

  type_last = type_current;	/* save type of current data word */


}
