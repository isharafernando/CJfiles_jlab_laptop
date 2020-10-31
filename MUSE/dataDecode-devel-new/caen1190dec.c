#include <stdio.h>
#include <stdint.h>
#include "caen1190dec.h"

void
caen1190DataDecode(uint32_t data)
{
  caen1190_generic_word_t gword;

  gword.raw = data;

  switch(gword.bf.word_type)
    {
    case CAEN1190_DATATYPE_GLOBAL_HEADER:
      {
	caen1190_global_header_t d; d.raw = data;

	printf("%8X - GLOBAL HEADER - slot = %2d  event count = %d\n",
	       d.raw,
	       d.bf.geo,
	       d.bf.event_count);
	break;
      }

    case CAEN1190_DATATYPE_TDC_HEADER:
      {
	caen1190_tdc_header_t d; d.raw = data;

	printf("%8X - TDC HEADER - tdc = %d  event ID = %4d  bunch ID = %4d\n",
	       d.raw,
	       d.bf.tdc,
	       d.bf.event_id,
	       d.bf.bunch_id);
	break;
      }

    case CAEN1190_DATATYPE_MEASUREMENT:
      {
	caen1190_tdc_measurement_t d; d.raw = data;

	printf("%8X - MEASUREMENT - edge = %d  channel = %3d  time = %6d\n",
	       d.raw,
	       d.bf.edge,
	       d.bf.channel,
	       d.bf.measurement);
	break;
      }

    case CAEN1190_DATATYPE_TDC_TRAILER:
      {
	caen1190_tdc_trailer_t d; d.raw= data;

	printf("%8X - TDC TRAILER - tdc = %d  event ID = %4d  word count = %4d\n",
	       d.raw,
	       d.bf.tdc,
	       d.bf.event_id,
	       d.bf.word_count);
	break;
      }

    case CAEN1190_DATATYPE_TDC_ERROR:
      {
	caen1190_tdc_error_t d; d.raw = data;

	printf("%8X - TDC ERROR - tdc = %d  error flags = 0X%04X\n",
	       d.raw,
	       d.bf.tdc,
	       d.bf.error_flags);
	break;
      }

    case CAEN1190_DATATYPE_EXTENDED_TRIGGER_TIME:
      {
	caen1190_extended_trigger_time_t d; d.raw = data;

	printf("%8X - EXTENDED TRIGGER TIME - time = %d\n",
	       d.raw,
	       d.bf.trigger_time);
	break;
      }

    case CAEN1190_DATATYPE_GLOBAL_TRAILER:
      {
	caen1190_global_trailer_t d; d.raw = data;

	printf("%8X - GLOBAL TRAILER - slot = %2d  word count = %5d  Error Flags = %d%d%d\n",
	       d.raw,
	       d.bf.geo,
	       d.bf.word_count,
	       d.bf.tdc_error,
	       d.bf.output_buffer_overflow,
	       d.bf.trigger_lost);
	break;
      }

    case CAEN1190_DATATYPE_FILLER:
      {
	caen1190_filler_t d; d.raw = data;
	printf("%8X - FILLER WORD = %d\n",
	       d.raw,
	       d.bf.word_type);
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
