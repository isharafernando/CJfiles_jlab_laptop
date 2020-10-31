#ifndef __CAEN1190DEC__
#define __CAEN1190DEC__
#include <stdint.h>

#define CAEN1190_DATATYPE_GLOBAL_HEADER          8
#define CAEN1190_DATATYPE_TDC_HEADER             1
#define CAEN1190_DATATYPE_MEASUREMENT            0
#define CAEN1190_DATATYPE_TDC_TRAILER            3
#define CAEN1190_DATATYPE_TDC_ERROR              4
#define CAEN1190_DATATYPE_EXTENDED_TRIGGER_TIME 17
#define CAEN1190_DATATYPE_GLOBAL_TRAILER        16
#define CAEN1190_DATATYPE_FILLER                24

typedef struct
{
  uint32_t undef:27;
  uint32_t word_type:5;
} caen1190_generic_word;

typedef union
{
  uint32_t raw;
  caen1190_generic_word bf;
} caen1190_generic_word_t;

typedef struct
{
  uint32_t geo:5;
  uint32_t event_count:22;
  uint32_t word_type:5;
} caen1190_global_header;

typedef union
{
  uint32_t raw;
  caen1190_global_header bf;
} caen1190_global_header_t;

typedef struct
{
  uint32_t bunch_id:12;
  uint32_t event_id:12;
  uint32_t tdc:2;
  uint32_t undef:1;
  uint32_t word_type:5;
} caen1190_tdc_header;

typedef union
{
  uint32_t raw;
  caen1190_tdc_header bf;
} caen1190_tdc_header_t;

typedef struct
{
  uint32_t measurement:19;
  uint32_t channel:7;
  uint32_t edge:1;
  uint32_t word_type:5;
} caen1190_tdc_measurement;

typedef union
{
  uint32_t raw;
  caen1190_tdc_measurement bf;
} caen1190_tdc_measurement_t;

typedef struct
{
  uint32_t word_count:12;
  uint32_t event_id:12;
  uint32_t tdc:2;
  uint32_t undef:1;
  uint32_t word_type:5;
} caen1190_tdc_trailer;

typedef union
{
  uint32_t raw;
  caen1190_tdc_trailer bf;
} caen1190_tdc_trailer_t;

typedef struct
{
  uint32_t error_flags:15;
  uint32_t undef0:9;
  uint32_t tdc:2;
  uint32_t undef1:1;
  uint32_t word_type:5;
} caen1190_tdc_error;

typedef union
{
  uint32_t raw;
  caen1190_tdc_error bf;
} caen1190_tdc_error_t;

typedef struct
{
  uint32_t trigger_time:27;
  uint32_t word_type:5;
} caen1190_extended_trigger_time;

typedef union
{
  uint32_t raw;
  caen1190_extended_trigger_time bf;
} caen1190_extended_trigger_time_t;

typedef struct
{
  uint32_t geo:5;
  uint32_t word_count:16;
  uint32_t undef:3;
  uint32_t tdc_error:1;
  uint32_t output_buffer_overflow:1;
  uint32_t trigger_lost:1;
  uint32_t word_type:5;
} caen1190_global_trailer;

typedef union
{
  uint32_t raw;
  caen1190_global_trailer bf;
} caen1190_global_trailer_t;

typedef struct
{
  uint32_t undef:27;
  uint32_t word_type:5;
} caen1190_filler;

typedef union
{
  uint32_t raw;
  caen1190_filler bf;
} caen1190_filler_t;

void caen1190DataDecode(uint32_t data);
#endif /* __CAEN1190DEC__ */
