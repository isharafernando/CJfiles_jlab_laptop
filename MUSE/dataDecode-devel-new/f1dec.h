#ifndef __F1DEC__
#define __F1DEC__
#include <stdint.h>

/* 1: BLOCK TRAILER */
typedef struct
{
  uint32_t words_in_block:22;
  uint32_t slot_number:5;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} f1_block_trailer;

typedef union
{
  uint32_t raw;
  f1_block_trailer bf;
} f1_block_trailer_t;

/* 3: TRIGGER TIME */
typedef struct
{
  uint32_t T_E:8;
  uint32_t T_D:8;
  uint32_t T_C:8;
  uint32_t undef:3;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} f1_trigger_time_1;

typedef union
{
  uint32_t raw;
  f1_trigger_time_1 bf;
} f1_trigger_time_1_t;

typedef struct
{
  uint32_t T_B:8;
  uint32_t T_A:8;
  uint32_t undef:15;
  uint32_t data_type_defining:1;
} f1_trigger_time_2;

typedef union
{
  uint32_t raw;
  f1_trigger_time_2 bf;
} f1_trigger_time_2_t;

/* 7: TIME MEASUREMENT */
typedef struct
{
  uint32_t time_measurement:16;
  uint32_t chip_channel_number:3;
  uint32_t chip_number:3;
  uint32_t fake:1;
  uint32_t undef:1;
  uint32_t hit_fifo_overflow:1;
  uint32_t output_fifo_overflow:1;
  uint32_t resolution_locked:1;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} f1_time_measurement;

typedef union
{
  uint32_t raw;
  f1_time_measurement bf;
} f1_time_measurement_t;


/* 8: CHIP HEADER */
typedef struct
{
  uint32_t chip_channel_number:3;
  uint32_t chip_number:3;
  uint32_t setup_register_tag:1;
  uint32_t trigger_time:9;
  uint32_t trigger_number:6;
  uint32_t trigger_fifo_overflow:1;
  uint32_t undef:1;
  uint32_t hit_fifo_overflow:1;
  uint32_t output_fifo_overflow:1;
  uint32_t resolution_locked:1;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} f1_chip_header;

typedef union
{
  uint32_t raw;
  f1_chip_header bf;
} f1_chip_header_t;



void f1DataDecode(uint32_t data);

#endif /* __F1DEC__ */
