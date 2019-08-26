#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "util.h"
#include "ppm.h"

unsigned char *load_ppm(char *filename, int *width, int *height)
{
    FILE *fp;
    char buf[1024];
    unsigned char *img;
    int i = 0, j, w, h;
    
    if((fp = fopen(filename, "rb")) == NULL)
    {
        fprintf(stderr, "failed to open file: %s\n", filename); return NULL;
    }

    while(fgets(buf, 1024, fp) != NULL) /* this is not a good parser */
    {
        if(buf[0] == '#' || buf[0] == '\n') ;
        else if(i == 0)
        {
            if(buf[0] == 'P' && buf[1] == '6') i++;
            else { fprintf(stderr, "file must be P6 ppm\n"); return NULL; }
        }
        else if(i == 1)
        {
            if(sscanf(buf, "%d %d", &w, &h) == 2)
            {
                i++;
                if((img = (unsigned char *)malloc(w * h * 3 * sizeof(unsigned char))) == NULL)
                {
                    fprintf(stderr, "failed to allocate memory in load_ppm\n"); return NULL;
                }
            }
            else { fprintf(stderr, "there must be width & height\n"); return NULL; }
        }
        else if(i == 2)
        {
            if(sscanf(buf, "%d", &j) == 1 && j == 255) { i++; break; }
            else { fprintf(stderr, "max level must be 255\n"); return NULL; }
        }
    }
    
    fread(img, sizeof(unsigned char), w * h * 3, fp);
    fclose(fp);

    *width = w;
    *height = h;
    return img;
}

unsigned char *load_pgm(char *filename, int *width, int *height)
{
    FILE *fp;
    char buf[1024];
    unsigned char *img;
    int i = 0, j, w, h;
    
    if((fp = fopen(filename, "rb")) == NULL)
    {
        fprintf(stderr, "failed to open file: %s\n", filename); return NULL;
    }

    while(fgets(buf, 1024, fp) != NULL) /* this is not a good parser */
    {
        if(buf[0] == '#' || buf[0] == '\n') ;
        else if(i == 0)
        {
            if(buf[0] == 'P' && buf[1] == '5') i++;
            else { fprintf(stderr, "file must be P5 pgm\n"); return NULL; }
        }
        else if(i == 1)
        {
            if(sscanf(buf, "%d %d", &w, &h) == 2)
            {
                i++;
                if((img = (unsigned char *)malloc(w * h * sizeof(unsigned char))) == NULL)
                {
                    fprintf(stderr, "failed to allocate memory in load_pgm\n"); return NULL;
                }
            }
            else { fprintf(stderr, "there must be width & height\n"); return NULL; }
        }
        else if(i == 2)
        {
            if(sscanf(buf, "%d", &j) == 1 && j == 255) { i++; break; }
            else { fprintf(stderr, "max level must be 255\n"); return NULL; }
        }
    }
    
    fread(img, sizeof(unsigned char), w * h, fp);
    fclose(fp);

    *width = w;
    *height = h;
    return img;
}


float ***load_ppm_normalized_rgb(char *filename, int *width, int *height)
{
    int i, j, k;
    unsigned char *tmp = load_ppm(filename, width, height);
    if(tmp == NULL) return NULL;
    
    float ***img = malloc_float3d(3, *width, *height);

    for(k = 0; k < 3; k++)
        for(i = 0; i < *width; i++)
            for(j = 0; j < *height; j++) img[k][i][j] = (float)tmp[(j * *width + i) * 3 + k] / 255.0;
    free(tmp);
    return img;
}

float **load_pgm_normalized_grayscale(char *filename, int *width, int *height)
{
    int i, j;
    
    unsigned char *tmp = load_pgm(filename, width, height);
    if(tmp == NULL) return NULL;
    
    float **img = malloc_float2d(*width, *height);

    for(i = 0; i < *width; i++)
        for(j = 0; j < *height; j++) img[i][j] = (float)tmp[j * *width + i] / 255.0;
    free(tmp);
    return img;
}

unsigned char **load_pgm_uchar2d(char *filename, int *width, int *height)
{
    int i, j;
    
    unsigned char *tmp = load_pgm(filename, width, height);
    if(tmp == NULL) return NULL;
    
    unsigned char **img = malloc_uchar2d(*width, *height);
    
    for(i = 0; i < *width; i++)
        for(j = 0; j < *height; j++) img[i][j] = tmp[j * *width + i];
    free(tmp);
    return img;
}

unsigned char ***load_ppm_uchar3d(char *filename, int *width, int *height)
{
    int i, j, k;
    
    unsigned char *tmp = load_ppm(filename, width, height);
    if(tmp == NULL) return NULL;
    
    unsigned char ***img = malloc_uchar3d(3, *width, *height);

    for(k = 0; k < 3; k++)
        for(i = 0; i < *width; i++)
            for(j = 0; j < *height; j++) img[k][i][j] = tmp[(j * *width + i) * 3 + k];
    free(tmp);
    return img;
}


void store_ppm(char *filename, int width, int height, unsigned char *image)
{
    FILE *fp;

    if((fp = fopen(filename, "wb")) == NULL)
    {
        fprintf(stderr, "failed to open file: %s\n", filename); exit(1);
    }

    fprintf(fp, "P6\n%d %d\n255\n", width, height);
    fwrite(image, sizeof(unsigned char), width * height * 3, fp);
    fclose(fp);
}

void store_pgm(char *filename, int width, int height, unsigned char *image)
{
    FILE *fp;

    if((fp = fopen(filename, "wb")) == NULL)
    {
        fprintf(stderr, "failed to open file: %s\n", filename); exit(1);
    }

    fprintf(fp, "P5\n%d %d\n255\n", width, height);
    fwrite(image, sizeof(unsigned char), width * height, fp);
    fclose(fp);
}


void store_ppm_normalized_rgb(char *filename, int width, int height, float ***image)
{
    int i, j, k, n;
    unsigned char *tmp;

    if((tmp = (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char))) == NULL)
    {
        fprintf(stderr, "failed to allocate memory\n"); exit(1);
    }
    for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
    for(k = 0; k < 3; k++)
    {
        n = (int)(image[k][i][j] * 255.0);
        if(n < 0) n = 0;
        else if(n > 255) n = 255;
        tmp[(j * width + i) * 3 + k] = (unsigned char)n;
    }
    store_ppm(filename, width, height, tmp);
    free(tmp);
}

void store_pgm_normalized_grayscale(char *filename, int width, int height, float **image)
{
    int i, j, n;
    unsigned char uc, *tmp;

    if((tmp = (unsigned char *)malloc(width * height * sizeof(unsigned char))) == NULL)
    {
        fprintf(stderr, "failed to allocate memory\n"); exit(1);
    }
    for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
    {
        n = (int)(image[i][j] * 255.0);
        if(n < 0) n = 0;
        else if(n > 255) n = 255;
        uc = (unsigned char)n;
        tmp[j * width + i] = uc;
    }
    store_pgm(filename, width, height, tmp);
    free(tmp);
}

void store_pgm_uchar2d(char *filename, int width, int height, unsigned char **image)
{
    int i, j;
    unsigned char *tmp;

    if((tmp = (unsigned char *)malloc(width * height * sizeof(unsigned char))) == NULL)
    {
        fprintf(stderr, "failed to allocate memory\n"); exit(1);
    }
    for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
    {
        tmp[j * width + i] = image[i][j];
    }
    store_pgm(filename, width, height, tmp);
    free(tmp);
}

void store_ppm_uchar3d(char *filename, int width, int height, unsigned char ***image)
{
    int i, j, k;
    unsigned char *tmp;

    if((tmp = (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char))) == NULL)
    {
        fprintf(stderr, "failed to allocate memory\n"); exit(1);
    }
    for(k = 0; k < 3; k++)
    for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
    {
        tmp[(j * width + i) * 3 + k] = image[k][i][j];
    }
    store_ppm(filename, width, height, tmp);
    free(tmp);
}
