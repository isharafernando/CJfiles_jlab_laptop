#include <stdio.h>
#include <stdint.h>
#include "jlabdec.h"
#include "tidec.h"

void
tiDataDecode(uint32_t data)
{
  static uint8_t blocklevel = 0;
  static uint8_t ievent = 0;
  static uint16_t event_length = 0;
  static uint32_t ievent_word = 0;
  generic_data_word_t gword;
  ti_event_header_t ehead;
  trigger_bank_t tbank;
  
  gword.raw = data;
  /* Check if the word is a block header */
  if((gword.bf.data_type_defining && (gword.bf.data_type_tag == 0))
     && !((ievent + 1) < blocklevel))
    {

      block_header_t d; d.raw = data;
      
      printf("%8X - BLOCK HEADER - slot = %d  modID = %d   n_evts = %d   n_blk = %d\n",
	     d.raw,
	     d.bf.slot_number,
	     d.bf.module_ID,
	     d.bf.number_of_events_in_block,
	     d.bf.event_block_number);
      
      blocklevel = d.bf.number_of_events_in_block;
      ievent_word = 0;
      event_length = 0;
      return;
    }
  
  tbank.raw = data;
  /* Check if the word is a trigger bank header */
  if(((tbank.bf.tag & 0xff10) == 0xff10)
     && (tbank.bf.type == 0x20)
     && !((ievent + 1) < blocklevel))
    {
      printf("%8X - TRIGGER BANK HEADER - timestamp = %d  blocklevel = %d\n",
	     tbank.raw,
	     tbank.bf.tag & 0xf,
	     tbank.bf.number_of_events);
      
      blocklevel = tbank.bf.number_of_events;
      ievent_word = 0;
      event_length = 0;
      ievent = 0;
      return;
    }
  
  ehead.raw = data;
  /* Check if the word is an event header */
  if(ehead.bf.type == 0x01)
    {
      printf("%8X - EVENT HEADER (%3d / %3d) - trigtype = %d  len = %d\n",
	     ehead.raw,
	     ievent + 1,
	     blocklevel,
	     ehead.bf.trigger_type,
	     ehead.bf.length);
      
      ievent_word = 0;
      ievent++;
      event_length = ehead.bf.length;
      return;
    }
  
  /* Assume anything else is trigger/event data */
  if((blocklevel > 0) && (event_length > 0))
    {
      if(ievent_word < event_length) /* Not sure if this needs to be <= */
	{
	  printf("%8X - EVENT DATA (%3d / %3d)\n",
		 data,
		 ievent_word+1, event_length);
	  ievent_word++;
	  
	  if(ievent_word == event_length)
	    {
	      ievent_word = 0;
	      event_length = 0;
	    }
	}
    }
  else
    {
      printf("%8X - EVENT DATA (\?\?)\n",
	     data);
    }
  
}
