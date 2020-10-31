#ifndef __CAEN792DEC__
#define __CAEN792DEC__
#include <stdint.h>

#define CAEN792_DATATYPE_HEADER        2
#define CAEN792_DATATYPE_DATA          0
#define CAEN792_DATATYPE_END_OF_BLOCK  4
#define CAEN792_DATATYPE_INVALID       6

typedef struct
{
  uint32_t undef:24;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen792_generic_word;

typedef union
{
  uint32_t raw;
  caen792_generic_word bf;
} caen792_generic_word_t;

typedef struct
{
  uint32_t undef0:8;
  uint32_t data_count:6;
  uint32_t undef1:2;
  uint32_t crate:8;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen792_header;

typedef union
{
  uint32_t raw;
  caen792_header bf;
} caen792_header_t;

typedef struct
{
  uint32_t adc:12;
  uint32_t overflow:1;
  uint32_t underflow:1;
  uint32_t undef:2;
  uint32_t channel_number:5;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen792_data;

typedef union
{
  uint32_t raw;
  caen792_data bf;
} caen792_data_t;

typedef struct
{
  uint32_t event_counter:24;
  uint32_t word_type:3;
  uint32_t geo:5;
} caen792_end_of_block;

typedef union
{
  uint32_t raw;
  caen792_end_of_block bf;
} caen792_end_of_block_t;

void caen792DataDecode(uint32_t data);
#endif /* __CAEN792DEC__ */
