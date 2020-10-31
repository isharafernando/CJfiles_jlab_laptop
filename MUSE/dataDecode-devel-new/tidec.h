#ifndef __TIDEC__
#define __TIDEC__
#include <stdint.h>

typedef struct
{
  uint32_t number_of_events:8;
  uint32_t type:8;
  uint32_t tag:16;
} trigger_bank;

typedef union
{
  uint32_t raw;
  trigger_bank bf;
} trigger_bank_t;

typedef struct
{
  uint32_t length:16;
  uint32_t type:8;
  uint32_t trigger_type:8;
} ti_event_header;

typedef union
{
  uint32_t raw;
  ti_event_header bf;
} ti_event_header_t;

void tiDataDecode(uint32_t data);

#endif /* __TIDEC__ */
