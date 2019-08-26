#ifndef __PPM_H__
#define __PPM_H__

#ifdef __cplusplus
extern "C" {
#endif

unsigned char *load_ppm(char *filename, int *width, int *height);
unsigned char *load_pgm(char *filename, int *width, int *height);

void store_ppm(char *filename, int width, int height, unsigned char *image);
void store_pgm(char *filename, int width, int height, unsigned char *image);

float ***load_ppm_normalized_rgb(char *filename, int *width, int *height);
float **load_pgm_normalized_grayscale(char *filename, int *width, int *height);
unsigned char **load_pgm_uchar2d(char *filename, int *width, int *height);
unsigned char ***load_ppm_uchar3d(char *filename, int *width, int *height);

void store_ppm_normalized_rgb(char *filename, int width, int height, float ***image);
void store_pgm_normalized_grayscale(char *filename, int width, int height, float **image);
void store_pgm_uchar2d(char *filename, int width, int height, unsigned char **image);
void store_ppm_uchar3d(char *filename, int width, int height, unsigned char ***image);

#ifdef __cplusplus
}
#endif

#endif
