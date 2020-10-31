#ifndef __FA125DEC__
#define __FA125DEC__
#include <stdint.h>

/* 2: EVENT HEADER */
typedef struct
{
  uint32_t trigger_number:22;
  uint32_t slot_number:5;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;  
} fa125_event_header;

typedef union
{
  uint32_t raw;
  fa125_event_header bf;
} fa125_event_header_t;

/* 3: TRIGGER TIME */
typedef struct
{
  uint32_t T_C:8;
  uint32_t T_B:8;
  uint32_t T_A:8;
  uint32_t undef:3;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} fa125_trigger_time_1;

typedef union
{
  uint32_t raw;
  fa125_trigger_time_1 bf;
} fa125_trigger_time_1_t;

typedef struct
{
  uint32_t T_F:8;
  uint32_t T_E:8;
  uint32_t T_D:8;
  uint32_t undef:7;
  uint32_t data_type_defining:1;
} fa125_trigger_time_2;

typedef union
{
  uint32_t raw;
  fa125_trigger_time_2 bf;
} fa125_trigger_time_2_t;


/* 4: WINDOW RAW DATA */
typedef struct
{
  uint32_t PTW:12;
  uint32_t undef:3;
  uint32_t slot_number:5;
  uint32_t channel_number:7;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} fa125_window_raw_data_1;

typedef union
{
  uint32_t raw;
  fa125_window_raw_data_1 bf;
} fa125_window_raw_data_1_t;

typedef struct
{
  uint32_t adc_sample_2:13;
  uint32_t invalid_2:1;
  uint32_t undef2:2;
  uint32_t adc_sample_1:13;
  uint32_t invalid_1:1;
  uint32_t undef1:1;
  uint32_t data_type_defining:1;
} fa125_window_raw_data_n;

typedef union
{
  uint32_t raw;
  fa125_window_raw_data_n bf;
} fa125_window_raw_data_n_t;

/* 5: PULSE DATA, CDC FORMAT */
typedef struct
{
  uint32_t overflow_count:3;
  uint32_t time_quality:1;
  uint32_t leading_edge_time:11;
  uint32_t NPK:5;
  uint32_t channel_number:7;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} fa125_pulse_cdc_1;

typedef union
{
  uint32_t raw;
  fa125_pulse_cdc_1 bf;
} fa125_pulse_cdc_1_t;

typedef struct
{
  uint32_t first_max_amplitude:9;
  uint32_t integral:14;
  uint32_t pedestal:8;
  uint32_t data_type_defining:1;
} fa125_pulse_cdc_2;

typedef union
{
  uint32_t raw;
  fa125_pulse_cdc_2 bf;
} fa125_pulse_cdc_2_t;

/* 6: PULSE DATA, FDC FORMAT - Integral + Time */
typedef struct
{
  uint32_t overflow_count:3;
  uint32_t time_quality:1;
  uint32_t leading_edge_time:11;
  uint32_t NPK:5;
  uint32_t channel_number:7;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} fa125_pulse_fdc_integral_1;

typedef union
{
  uint32_t raw;
  fa125_pulse_fdc_integral_1 bf;
} fa125_pulse_fdc_integral_1_t;

typedef struct
{
  uint32_t pedestal:11;
  uint32_t peak_time:8;
  uint32_t integral:12;
  uint32_t data_type_defining:1;
} fa125_pulse_fdc_integral_n;

typedef union
{
  uint32_t raw;
  fa125_pulse_fdc_integral_n bf;
} fa125_pulse_fdc_integral_n_t;

/* 9: PULSE DATA, FDC FORMAT - Amplitude + Time */
typedef struct
{
  uint32_t overflow_count:3;
  uint32_t time_quality:1;
  uint32_t leading_edge_time:11;
  uint32_t NPK:5;
  uint32_t channel_number:7;
  uint32_t data_type_tag:4;
  uint32_t data_type_defining:1;
} fa125_pulse_fdc_amplitude_1;

typedef union
{
  uint32_t raw;
  fa125_pulse_fdc_amplitude_1 bf;
} fa125_pulse_fdc_amplitude_1_t;

typedef struct
{
  uint32_t pedestal:11;
  uint32_t peak_time:8;
  uint32_t amplitude:12;
  uint32_t data_type_defining:1;
} fa125_pulse_fdc_amplitude_n;

typedef union
{
  uint32_t raw;
  fa125_pulse_fdc_amplitude_n bf;
} fa125_pulse_fdc_amplitude_n_t;

void  fa125DataDecode(uint32_t data);


#endif /* __FA125DEC__ */
