#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

#include <libusb.h>
#include "rtl-sdr.h"

#include "bmp.c"


//---------------------------------------------------------------
// benoit defs

#define uchar unsigned char
//---------------------------------------------------------------

typedef struct {
    float   real;
    float   imag;
} complex;

//---------------------------------------------------------------

void set_complex_int(complex *c, int r, int i)
{
    c->real = r / 32768.0;
    c->imag = i / 32768.0;
}

//---------------------------------------------------------------

void set_complex_float(complex *c, float r, float i)
{
    c->real = r;
    c->imag = i;
}

//---------------------------------------------------------------

void complex_from_angle(complex *c, float angle)
{
    c->real = cos(angle);
    c->imag = sin(angle);
}

//---------------------------------------------------------------

void complex_mult(complex *a, complex *b, complex *r)
{
    r->real = a->real * b->real - a->imag * b->imag;
    r->imag = a->imag * b->real + a->real * b->imag;
}

//---------------------------------------------------------------

float modulus(complex *v)
{
    return sqrt(v->real * v->real + v->imag * v->imag);
}

//---------------------------------------------------------------

float modulus_square(complex *v)
{
    return v->real * v->real + v->imag * v->imag;
}

//---------------------------------------------------------------

void normalize(complex *v)
{
    float norm  = 1.95 - ((v->real * v->real) + (v->imag * v->imag));
    
    v->real *= norm;
    v->imag *= norm;
}

//---------------------------------------------------------------

void complex_copy(complex *v1, complex *v2)
{
    v2->real = v1->real;
    v2->imag = v1->imag;
}

//---------------------------------------------------------------

void demod(short *input, float *output, int cnt, float frequency, float sample_rate)
{
    complex val;
    complex carrier;
    complex tmp;
    complex rotation;
    
    float   angle_per_sample;
    
    angle_per_sample = 2.0 * 3.1415926535 * -frequency / sample_rate;
    complex_from_angle(&rotation, angle_per_sample);
    
    set_complex_float(&carrier, 1.0, 0.0);
    
    for (int i = 0; i < cnt; i++) {
        short   i, q;

        i = *input++;
        q = *input++;
        
        set_complex_int(&val, i, q);
        
        complex_mult(&val, & carrier, &tmp);
        
        *output++ = modulus(&tmp);
        
        complex_mult(&carrier, &rotation, &tmp);
        normalize(&tmp);
        complex_copy(&tmp, &carrier);
        
    }
}

//---------------------------------------------------------------

void demodr(short *input, short *output, int cnt, float frequency, float sample_rate)
{
    complex val;
    complex carrier;
    complex tmp;
    complex rotation;
    
    float   angle_per_sample;
    
    angle_per_sample = 2.0 * 3.1415926535 * -frequency / sample_rate;
    complex_from_angle(&rotation, angle_per_sample);
    
    set_complex_float(&carrier, 1.0, 0.0);
    
    for (int i = 0; i < cnt; i++) {
        short   i, q;
        
        //printf("%f %f\n", carrier.real, carrier.imag);
        i = *input++;
        q = *input++;
        
        set_complex_int(&val, i, q);
        
        complex_mult(&val, & carrier, &tmp);
        
        *output++ = tmp.real*32768;
        *output++ = tmp.imag*32768;
        
        
        complex_mult(&carrier, &rotation, &tmp);
        normalize(&tmp);
        complex_copy(&tmp, &carrier);
    }
}

//---------------------------------------------------------------

int16_t* Sinewave;
int N_WAVE, LOG2_N_WAVE;
int16_t *fft_buf;

//---------------------------------------------------------------

void sine_table(int size) {
	int i;
	double d;
	LOG2_N_WAVE = size;
	N_WAVE = 1 << LOG2_N_WAVE;
	Sinewave = malloc(sizeof(int16_t) * N_WAVE * 3 / 4);
	for (i = 0; i < N_WAVE * 3 / 4; i++) {
		d = (double) i * 2.0 * M_PI / N_WAVE;
		Sinewave[i] = (int) round(32767 * sin(d));
	}
}

//---------------------------------------------------------------

int16_t FIX_MPY(int16_t a, int16_t b) {
	int c = ((int) a * (int) b) >> 14;
	b = c & 0x01;
	return (c >> 1) + b;
}

//---------------------------------------------------------------

int fix_fft(int16_t iq[], int m)
/* interleaved iq[], 0 <= n < 2**m, changes in place */
{
	int mr, nn, i, j, l, k, istep, n, shift;
	int16_t qr, qi, tr, ti, wr, wi;
	n = 1 << m;
	if (n > N_WAVE) {
		return -1;
	}
	mr = 0;
	nn = n - 1;
	/* decimation in time - re-order data */
	for (m = 1; m <= nn; ++m) {
		l = n;
		do {
			l >>= 1;
		} while (mr + l > nn);
		mr = (mr & (l - 1)) + l;
		if (mr <= m) {
			continue;
		}
		// real = 2*m, imag = 2*m+1
		tr = iq[2 * m];
		iq[2 * m] = iq[2 * mr];
		iq[2 * mr] = tr;
		ti = iq[2 * m + 1];
		iq[2 * m + 1] = iq[2 * mr + 1];
		iq[2 * mr + 1] = ti;
	}
	l = 1;
	k = LOG2_N_WAVE - 1;
	while (l < n) {
		shift = 1;
		istep = l << 1;
		for (m = 0; m < l; ++m) {
			j = m << k;
			wr = Sinewave[j + N_WAVE / 4];
			wi = -Sinewave[j];
			if (shift) {
				wr >>= 1;
				wi >>= 1;
			}
			for (i = m; i < n; i += istep) {
				j = i + l;
				tr = FIX_MPY(wr, iq[2 * j]) - FIX_MPY(wi, iq[2 * j + 1]);
				ti = FIX_MPY(wr, iq[2 * j + 1]) + FIX_MPY(wi, iq[2 * j]);
				qr = iq[2 * i];
				qi = iq[2 * i + 1];
				if (shift) {
					qr >>= 1;
					qi >>= 1;
				}
				iq[2 * j] = qr - tr;
				iq[2 * j + 1] = qi - ti;
				iq[2 * i] = qr + tr;
				iq[2 * i + 1] = qi + ti;
			}
		}
		--k;
		l = istep;
	}
	return 0;
}

//---------------------------------------------------------------

#define LSIZE	(8)
#define	SIZE	(1<<LSIZE)

//---------------------------------------------------------------


char i_to_s(int x) {
	static char intensity[] = "..:*@@@@@@@@@$$$$$$$$$$$$$$$$$$$$$$$";

	x = (x + 10) / 200;

	if (x > 15)
		x = 15;

	return intensity[x];
}

//---------------------------------------------------------------
// util junk
//---------------------------------------------------------------

int find_max(short *array) {
	float	max = 0;
	int	max_i = 0;

	for (int i = 0; i < SIZE / 2; i++) {
		int x = abs(array[i * 2]);
		if (x > max) {
			max = x;
			max_i = i;
		}
	}

	return max_i;
}

void display(short *array) {
	for (int i = 0; i < SIZE / 2; i++) {
		int x = abs(array[i * 2]);
		printf("%c", i_to_s(x));
	}
	printf("\n");
}

//---------------------------------------------------------------

void print_time() {
	struct timeval t;
	struct timezone z;
    
	gettimeofday(&t, &z);

    printf("%ld ", t.tv_sec);
}

//---------------------------------------------------------------

double nanotime() {
	struct timeval t;
	struct timezone z;

	gettimeofday(&t, &z);

	return t.tv_usec;
}

//---------------------------------------------------------------
// globals
//---------------------------------------------------------------

short data[SIZE * 2];
uchar input[SIZE * 2];
FILE *fp_in;
rtlsdr_dev_t *dev = NULL;


#define DEFAULT_SAMPLE_RATE             2048000

void init_sdr(int dev_index)
{
    char vendor[256], product[256], serial[256];
    

    int device_count = rtlsdr_get_device_count();
	if (!device_count) {
		printf("No supported devices found.\n");
		exit(1);
	}
    
    printf("Found %d device(s):\n", device_count);
    for (int i = 0; i < device_count; i++) {
        rtlsdr_get_device_usb_strings(i, vendor, product, serial);
        printf("  %d:  %s, %s, SN: %s\n", i, vendor, product, serial);
    }
    printf("\n");
    printf("Using device %d: %s\n", dev_index, rtlsdr_get_device_name(dev_index));
    
    int r = rtlsdr_open(&dev, dev_index);
    if (r < 0) {
		printf("devices index not found.\n");
		exit(1);
    }
    r = rtlsdr_set_sample_rate(dev, DEFAULT_SAMPLE_RATE);

    r = rtlsdr_set_center_freq(dev, 42500000);
    r = rtlsdr_set_tuner_gain_mode(dev, 1);
    r = rtlsdr_set_tuner_gain(dev, 456);
    r = rtlsdr_reset_buffer(dev);
 }

void init_iq_source(char *filename)
{
	fp_in = fopen(filename, "rb");
}

//---------------------------------------------------------------

int get_sdr_data(uchar *buffer, int icnt) {
    int n_read;
    
    int r = rtlsdr_read_sync(dev, buffer, 2*icnt, &n_read);
    return n_read;
}

//---------------------------------------------------------------

int get_iq_data_file(uchar *buffer, int icnt) {
	int cnt = fread(buffer, 1, icnt * 2, fp_in);
    fseek(fp_in,-(icnt * 2) + 12*128 , SEEK_CUR);
	return cnt;
}

//---------------------------------------------------------------

void iq_to_unsigned(uchar *in_buffer, short *cvt)
{
    for (int i = 0; i < SIZE; i++) {
        float v;
        
        v = in_buffer[i * 2];
        v = v - 127;
        v *= 200;
        cvt[i * 2] = v;
        v = in_buffer[1 + i * 2];
        v = v - 127;
        v *= 200;
        cvt[1 + i * 2] = v;
    }
}

//---------------------------------------------------------------

void clean_fft(short *data)
{
	data[0] = 0;		//clear dc
}


//---------------------------------------------------------------
float   freq = 0;
uchar    rgb[(SIZE*SIZE)/2 * 3];

//---------------------------------------------------------------

void set_line(int line, float *fft_val)
{
    uchar *tmp;
    
    tmp = &rgb[line * (SIZE/2) * 3];
    
    for (int i = 0; i < SIZE/2; i++) {
        int v = sqrt(abs(fft_val[i * 2] / (4.0)));
        if (v > 255) v = 255;
        *tmp++ = v;
        *tmp++ = v;
        *tmp++ = v;
    }

}

//---------------------------------------------------------------

int gmain(int argc, char **argv) {
    
    int     line = 0;
    float   fft_sum[SIZE];
    int     buffer_count = 0;

	printf("start\n");
	if (argc < 1) {
		printf("usage is fft <filename>\n");
		return -1;
	}
	sine_table(LSIZE);
    
	init_iq_source(argv[1]);
    
	int cnt = 0;
    buffer_count = 0;
    for (int i = 0; i < SIZE; i++) {
        fft_sum[i] = 0;
    }
    
	do {
        cnt = get_iq_data_file(input, SIZE);
        iq_to_unsigned(input, data);
		if (cnt == SIZE * 2) {
			fix_fft(data, LSIZE);
			clean_fft(data);
            for (int i = 0; i < SIZE/2; i++) {
                    fft_sum[i*2] += abs(data[i*2]);
            }
            buffer_count++;
            
            if (buffer_count == (512)) {
                set_line(line, fft_sum);
                line++;
                buffer_count = 0;
                for (int i = 0; i < SIZE; i++) {
                    fft_sum[i] = 0;
                }
            }
		}
	} while (line < SIZE/2 && cnt == SIZE * 2);
    
	fclose(fp_in);
    write_bmp("test.bmp", SIZE/2, SIZE/2, (char*)&rgb);

    return 0;
}


//---------------------------------------------------------------


//---------------------------------------------------------------

int main(int argc, char **argv) {
    
    int     line = 0;
    float   fft_sum[SIZE];
    int     buffer_count = 0;
    
    print_time();
    
	sine_table(LSIZE);
    init_sdr(0);
    
	int cnt = 0;
    buffer_count = 0;
    for (int i = 0; i < SIZE; i++) {
        fft_sum[i] = 0;
    }
    
	do {
        cnt = get_sdr_data(input, SIZE);
        iq_to_unsigned(input, data);
		if (cnt == SIZE * 2) {
			fix_fft(data, LSIZE);
			clean_fft(data);
            for (int i = 0; i < SIZE/2; i++) {
                fft_sum[i*2] += abs(data[i*2]);
            }
            buffer_count++;
            
            if (buffer_count == (8192)) {
                buffer_count = 0;
                print_time();
                for (int i = 0; i < SIZE/2; i++) {
                    printf("%f\t", sqrt(fft_sum[i*2]));
                }
                printf("\n");
                
                for (int i = 0; i < SIZE; i++) {
                    fft_sum[i] = 0;
                }
            }
		}
	} while (cnt == SIZE * 2);
    
	fclose(fp_in);
    
    return 0;
}


//---------------------------------------------------------------