#ifndef __LECROY1881DEC__
#define __LECROY1881DEC__
#include <stdint.h>

typedef struct
{
  uint32_t word_count:7;
  uint32_t buffer_number:6;
  uint32_t undef0:2;
  uint32_t word_parity:1;
  uint32_t undef1:11;
  uint32_t geo_address:5;
} lecroy1881_header;

typedef union
{
  uint32_t raw;
  lecroy1881_header bf;
} lecroy1881_header_t;

typedef struct
{
  uint32_t data:14;
  uint32_t undef0:3;
  uint32_t channel_number:7;
  uint32_t buffer_number:2; /* %4 */
  uint32_t word_parity:1;
  uint32_t geo_address:5;
} lecroy1881_data;

typedef union
{
  uint32_t raw;
  lecroy1881_data bf;
} lecroy1881_data_t;

void lecroy1881DataDecode(uint32_t data);

#endif /* __LECROY1881DEC__ */
