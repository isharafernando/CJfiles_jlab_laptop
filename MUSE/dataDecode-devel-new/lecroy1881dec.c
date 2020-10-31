#include <stdio.h>
#include <stdint.h>
#include "lecroy1881dec.h"

void
lecroy1881DataDecode(uint32_t data)
{
  static uint32_t nwords = 0, iword = 0;

  if((nwords == 0)
     || ((iword + 1) == nwords))
    {
      /* Assume we're at the beginning of module output.  data should
	 be the header */
      lecroy1881_header_t d; d.raw = data;

      printf("%8X - HEADER - G = %2d  P = %d  Buff = %d  Word Count = %d\n",
	     d.raw,
	     d.bf.geo_address,
	     d.bf.word_parity,
	     d.bf.buffer_number,
	     d.bf.word_count);

      iword = 0;
      nwords = d.bf.word_count;
    }
  else
    {
      lecroy1881_data_t d; d.raw = data;

      printf("%8X - DATA (%4d/%4d) - G = %2d  P = %d  Buff = %d  Chan %3d  data = %5d\n",
	     d.raw,
	     iword + 2,
	     nwords,
	     d.bf.geo_address,
	     d.bf.word_parity,
	     d.bf.buffer_number,
	     d.bf.channel_number,
	     d.bf.data);
      iword++;
    }
  
}
