#ifndef __F1V1DEC__
#define __F1V1DEC__
#include <stdint.h>

/* Module data word */
typedef struct
{
  uint32_t chip_word:23;
  uint32_t chip_word_type:1;
  uint32_t hit_fifo_overflow:1;
  uint32_t output_fifo_overflow:1;
  uint32_t res_locked:1;
  uint32_t slot_id:5;
} f1_data;

typedef union
{
  uint32_t raw;
  f1_data bf;
} f1_data_t;

/* F1 Chip Header/Trailer */
typedef struct
{
  uint32_t channel_address:3;
  uint32_t chip_address:3;
  uint32_t Xor_setup_register:1;
  uint32_t trigger_time:9;
  uint32_t trigger_number:6;
  uint32_t trigger_fifo_overflow:1;
} f1_chip_header;

typedef union
{
  uint32_t raw;
  f1_chip_header bf;
} f1_chip_header_t;

/* F1 Chip Data */
typedef struct
{
  uint32_t time_measurement:16;
  uint32_t channel_address:3;
  uint32_t chip_address:3;
  uint32_t zero:1;
} f1_chip_data;

typedef union
{
  uint32_t raw;
  f1_chip_data bf;
} f1_chip_data_t;


void f1v1DataDecode(uint32_t data);

#endif /* __F1DEC__ */
