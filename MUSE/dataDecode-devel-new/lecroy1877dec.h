#ifndef __LECROY1877DEC__
#define __LECROY1877DEC__
#include <stdint.h>

typedef struct
{
  uint32_t word_count:11;
  uint32_t buffer_number:3;
  uint32_t undef0:1;
  uint32_t word_parity:1;
  uint32_t undef1:11;
  uint32_t geo_address:5;
} lecroy1877_header;

typedef union
{
  uint32_t raw;
  lecroy1877_header bf;
} lecroy1877_header_t;

typedef struct
{
  uint32_t data:16;
  uint32_t leading_edge:1;
  uint32_t channel_number:7;
  uint32_t buffer_number:2; /* Hit Count (%4) in 1877s */
  uint32_t word_parity:1;
  uint32_t geo_address:5;
} lecroy1877_data;

typedef union
{
  uint32_t raw;
  lecroy1877_data bf;
} lecroy1877_data_t;

void lecroy1877DataDecode(uint32_t data);

#endif /* __LECROY1877DEC__ */
