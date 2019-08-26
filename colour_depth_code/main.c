#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "util.h"
#include "ppm.h"
#include "ii.h"
#include "mrf.h"
#include "GCoptimization.h"

#define PRECISION (8388608.0)

static void make_offset(int w, int h, int channel, int offset, int in_x, int in_y, int *out_x, int *out_y);
static void graph_cut(int w, int h, int num_label, float lambda, float ***cost, float ***img, int **result);
static float ***downsample_from_ii(int w, int h, ullong ***ii_rgb, int size);

static void print_help()
{
    fprintf(stderr, "main in.ppm size disp_min disp_max\n");
    fprintf(stderr, "-d int  : downsample size\n");
    fprintf(stderr, "-s float: smoothness weight\n");
    exit(0);
}

int	main(int argc, char *argv[])
{
    int i, j, k, w, h;
    float ***img;
    
	if(argc < 5) print_help();
    int n = 5;
    int downsample_size = 1;
    float smooth_weight = 1.0;
    while(n < argc)
    {
        if(strcmp(argv[n], "-d") == 0 && n+1 < argc) { downsample_size = atoi(argv[n+1]); n+=2; }
        else if(strcmp(argv[n], "-s") == 0 && n+1 < argc) { smooth_weight = (float)atof(argv[n+1]); n+=2; }
        else print_help();
    }

    const int size = atoi(argv[2]);
    const int disp_min = atoi(argv[3]);
    const int disp_max = atoi(argv[4]);
    const int num_label = disp_max - disp_min + 1;
    const int max_abs_disp = disp_max > -disp_min ? disp_max : -disp_min;

    if((img = load_ppm_normalized_rgb(argv[1], &w, &h)) == NULL) exit(1);
    const int ds_w = w / downsample_size;
    const int ds_h = h / downsample_size;

    int ***rgb   = malloc_int3d(3, w, h);
    int ***rgb2  = malloc_int3d(3, w, h);
    int ***cross = malloc_int3d(3, w, h);
    ullong ***ii_rgb   = malloc_ullong3d(3, w, h);
    ullong ***ii_rgb2  = malloc_ullong3d(3, w, h);
    ullong ***ii_cross = malloc_ullong3d(3, w, h);

    for(k = 0; k < 3; k++)
    {
        for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
        {
            rgb[k][i][j]   = (int)(PRECISION * img[k][i][j]);
            rgb2[k][i][j]  = (int)(PRECISION * img[k][i][j] * img[k][i][j]);
        }

        generate_ii(w, h, rgb[k],  ii_rgb[k]);
        generate_ii(w, h, rgb2[k], ii_rgb2[k]);
    }

    float ***err = malloc_float3d(ds_h, ds_w, num_label);
    memset(err[0][0], 0, ds_h * ds_w * num_label * sizeof(float));
    float **min_err = malloc_float2d(w, h);
    int **min_label = malloc_int2d(w, h);
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        min_err[i][j] = 1.0e+20;
        min_label[i][j] = INT_MAX;
    }

    const double area = (double)((size * 2 + 1) * (size * 2 + 1));
    const double denom = 1.0 / (PRECISION * area);

    int ofs;
    for(ofs = disp_min; ofs <= disp_max; ofs++)
    {
        for(k = 0; k < 3; k++)
        {
            for(i = 0; i < w; i++)
            for(j = 0; j < h; j++)
            {
                int x, y;
                make_offset(w, h, k, ofs, i, j, &x, &y);

                int m = (k + 1) % 3;
                int u, v;
                make_offset(w, h, m, ofs, i, j, &u, &v);
                
                cross[k][i][j] = (int)(PRECISION * img[k][x][y] * img[m][u][v]);
            }
            
            generate_ii(w, h, cross[k], ii_cross[k]);
        }

        for(i = size + max_abs_disp; i < w - size - max_abs_disp; i++)
        for(j = size + max_abs_disp; j < h - size - max_abs_disp; j++)
        {
            int x1 = i - size;
            int y1 = j - size;
            int x2 = i + size;
            int y2 = j + size;
            int u1[3], v1[3], u2[3], v2[3];
            for(k = 0; k < 3; k++)
            {
                make_offset(w, h, k, ofs, x1, y1, &u1[k], &v1[k]);
                make_offset(w, h, k, ofs, x2, y2, &u2[k], &v2[k]);
            }

            double c[3], c2[3], cr[3];
            for(k = 0; k < 3; k++)
            {
                c[k]  = sample_ii(w, h, ii_rgb[k], u1[k], v1[k], u2[k], v2[k]) * denom;
                c2[k] = sample_ii(w, h, ii_rgb2[k], u1[k], v1[k], u2[k], v2[k]) * denom;
                cr[k] = sample_ii(w, h, ii_cross[k], x1, y1, x2, y2) * denom;
            }

            double cov[3][3];
            cov[0][0] = c2[0] - c[0] * c[0];
            cov[1][1] = c2[1] - c[1] * c[1];
            cov[2][2] = c2[2] - c[2] * c[2];
            cov[0][1] = cov[1][0] = cr[0] - c[0] * c[1];
            cov[0][2] = cov[2][0] = cr[2] - c[0] * c[2];
            cov[1][2] = cov[2][1] = cr[1] - c[1] * c[2];
            for(k = 0; k < 3; k++)
                if(cov[k][k] < 0) cov[k][k] = 0;

            double det = cov[0][0] * cov[1][1] * cov[2][2]
                       + cov[0][1] * cov[1][2] * cov[2][0]
                       + cov[0][2] * cov[1][0] * cov[2][1]
                       - cov[0][2] * cov[1][1] * cov[2][0]
                       - cov[0][1] * cov[1][0] * cov[2][2]
                       - cov[0][0] * cov[1][2] * cov[2][1];
            float e = (float)(det / (cov[0][0] * cov[1][1] * cov[2][2] + 1.0e-10));

            if(min_err[i][j] > e)
            {
                min_err[i][j] = e;
                min_label[i][j] = ofs;
            }
            err[j/downsample_size][i/downsample_size][ofs - disp_min] += e;
        }
    }

    float ***tmp = malloc_float3d(3, w, h);
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        if(min_label[i][j] == INT_MAX || min_label[i][j] == disp_min || min_label[i][j] == disp_max)
        {
            for(k = 0; k < 3; k++) tmp[k][i][j] = (k == 2 ? 0.5 : 0);
            for(k = 0; k < num_label; k++) err[j/downsample_size][i/downsample_size][k] = 0;
        }
        else
        {
            for(k = 0; k < 3; k++) tmp[k][i][j] = (min_label[i][j] - disp_min) / (float)(num_label - 1);
        }
    }
    store_ppm_normalized_rgb("depth_local.ppm", w, h, tmp);

    float ***ds_img = downsample_from_ii(w, h, ii_rgb, downsample_size);
    int **result = malloc_int2d(ds_w, ds_h);
    graph_cut(ds_w, ds_h, num_label, smooth_weight, err, ds_img, result);

    unsigned char **cut = malloc_uchar2d(ds_w, ds_h);
    for(i = 0; i < ds_w; i++)
    for(j = 0; j < ds_h; j++)
    {
        cut[i][j] = (unsigned char)(result[i][j] * 255 / (num_label - 1));
    }
    store_pgm_uchar2d("depth_global.pgm", ds_w, ds_h, cut);

    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        int disparity = (int)(result[i/downsample_size][j/downsample_size] + 0.5) + disp_min;
        int x, y;
        for(k = 0; k < 3; k++)
        {
            make_offset(w, h, k, disparity, i, j, &x, &y);
            tmp[k][i][j] = img[k][x][y];
        }
    }
    store_ppm_normalized_rgb("restored.ppm", w, h, tmp);

    free_float3d(img);
    free_int3d(rgb);
    free_int3d(rgb2);
    free_int3d(cross);
    free_ullong3d(ii_rgb);
    free_ullong3d(ii_rgb2);
    free_ullong3d(ii_cross);    
    return 0;
}
// R, G, and B planes are shifted to rightward, upward, and leftward
static void make_offset(int w, int h, int channel, int offset, int in_x, int in_y, int *out_x, int *out_y)
{
    int x, y;
    switch(channel)
    {
    case 0: // R
        x = in_x + offset;
        y = in_y;
        break;
        
    case 1: // G
        x = in_x;
        y = in_y - offset;
        break;
        
    default: // B
        x = in_x - offset;
        y = in_y;
        break;
    }
    
    if(x < 0) x = 0; else if(x >= w) x = w - 1;
    if(y < 0) y = 0; else if(y >= h) y = h - 1;

    *out_x = x;
    *out_y = y;
}
// Sum of the square disparity of 3 channels
static float dist2(float *rgb1, float *rgb2)
{
    int i;
    float d = 0;
    for(i = 0; i < 3; i++) d += (rgb1[i] - rgb2[i]) * (rgb1[i] - rgb2[i]);
    return d;
}
static float set_weight(float *rgb1, float *rgb2, float scale)
{
    const float e = 0.1;
    return (e + expf(-scale * dist2(rgb1, rgb2))) / (1.0 + e);
}
// the standard energy minimization framework
static void graph_cut(int w, int h, int num_label, float lambda, float ***cost, float ***img, int **result)
{
    int i, j;
    float e, prev_e;
    float t, tot_t;

    int n = 0;
    float sigma2 = 0;
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        if(i < w-1) { sigma2 += dist2(img[i][j], img[i+1][j]); n++; }
        if(j < h-1) { sigma2 += dist2(img[i][j], img[i][j+1]); n++; }
    }
    const float scale = (float)n / (2.0 * sigma2);

    float **h_weights = malloc_float2d(h, w);
    float **v_weights = malloc_float2d(h, w);
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        if(i < w-1) h_weights[j][i] = set_weight(img[i][j], img[i+1][j], scale);
        if(j < h-1) v_weights[j][i] = set_weight(img[i][j], img[i][j+1], scale);
    }
    
    DataCost *data         = new DataCost(cost[0][0]);
    SmoothnessCost *smooth = new SmoothnessCost(1, 1.0, lambda, h_weights[0], v_weights[0]);
    EnergyFunction *eng    = new EnergyFunction(data, smooth);
    MRF* mrf               = new Expansion(w, h, num_label, eng);
    mrf->initialize();
    mrf->clearAnswer();

    tot_t = 0;
    e = 1.0e+20;
    do
    {
        mrf->optimize(1, t);
        tot_t = tot_t + t;
        prev_e = e;
        e = mrf->totalEnergy();
        printf("energy = %e (%f secs)\n", e, tot_t);
    }
    while(e < prev_e);

    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++) result[i][j] = mrf->getLabel(j * w + i);

    delete mrf;
}
//down sample
static float ***downsample_from_ii(int w, int h, ullong ***ii_rgb, int size)
{
    int i, j, k;
    int ds_w = w / size;
    int ds_h = h / size;
    float ***img = malloc_float3d(ds_w, ds_h, 3);
    const float denom = 1.0 / (PRECISION * (float)(size * size));
    for(i = 0; i < ds_w; i++)
    for(j = 0; j < ds_h; j++)
    for(k = 0; k < 3; k++)
    {
        img[i][j][k] = sample_ii(w, h, ii_rgb[k], i * size, j * size, (i + 1) * size - 1, (j + 1) * size - 1) * denom;
    }
    return img;
}
