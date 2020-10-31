#ifndef __CAEN775DEC__
#define __CAEN775DEC__
#include <stdint.h>

#define CAEN775_DATATYPE_HEADER        2
#define CAEN775_DATATYPE_DATA          0
#define CAEN775_DATATYPE_END_OF_BLOCK  4
#define CAEN775_DATATYPE_INVALID       6

typedef struct
{
  uint32_t undef:24;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen775_generic_word;

typedef union
{
  uint32_t raw;
  caen775_generic_word bf;
} caen775_generic_word_t;

typedef struct
{
  uint32_t undef0:8;
  uint32_t data_count:6;
  uint32_t undef1:2;
  uint32_t crate:8;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen775_header;

typedef union
{
  uint32_t raw;
  caen775_header bf;
} caen775_header_t;

typedef struct
{
  uint32_t adc:12;
  uint32_t overflow:1;
  uint32_t underflow:1;
  uint32_t valid:1;
  uint32_t undef:1;
  uint32_t channel_number:5;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen775_data;

typedef union
{
  uint32_t raw;
  caen775_data bf;
} caen775_data_t;

typedef struct
{
  uint32_t event_counter:24;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen775_end_of_block;

typedef union
{
  uint32_t raw;
  caen775_end_of_block bf;
} caen775_end_of_block_t;

void caen775DataDecode(uint32_t data);

#endif /* __CAEN775DEC__ */
