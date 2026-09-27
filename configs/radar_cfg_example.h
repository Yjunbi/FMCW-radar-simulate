/* Entirely synthetic teaching parameters. No physical device is represented.
 * Frequency: kHz for start frequency, MHz for bandwidth; times: microseconds.
 * This minimal initializer is input data for the Python examples only.
 */
#define SYNTHETIC_TEACHING_WAVE { \
    .chirp_startfreq_khz = 60000000U, \
    .chirp_bandwidth = 300.0f, \
    .chirp_rise = 40.0f, \
    .chirp_fall = 4U, \
    .chirp_idle = 6.0f, \
    .tx_phase_start = {0x10U, 0x50U, 0xA0U, 0xE0U}, \
    .tx_phase_step = {0x00U, 0x30U, 0x70U, 0xD0U}, \
    .adc_sample_rate = RADAR_ADC_RATE_10M, \
    .adc_chirp_size = 256U, \
    .adc_sample_delay = 3U, \
    .range_fft_size = 256U, \
    .range_win_en = 1U, \
    .range_win_type = WIN_BLACKMAN, \
    .num_rangebin = 128U, \
    .doppler_fft_size = 256U, \
    .doppler_win_en = 1U, \
    .doppler_win_type = WIN_BLACKMAN, \
}
