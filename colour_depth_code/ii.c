#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "ii.h"

void generate_ii(int w, int h, int **in, ullong **out)
{
    int i, j;

    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        out[i][j] = (j == 0 ? 0 : out[i][j-1]) + in[i][j];
    }
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        out[i][j] += (i == 0 ? 0 : out[i-1][j]);
    }
}

int sample_ii(int w, int h, ullong **ii, int x1, int y1, int x2, int y2)
{
    if(x1 < 0 || y1 < 0 || x2 >= w || y2 >= h || x1 > x2 || y1 > y2)
    {
        fprintf(stderr, "invalid area (%d, %d) - (%d, %d)\n", x1, y1, x2, y2); exit(1);
    }
    
    x1--;
    y1--;

    ullong ret = ii[x2][y2];
    if(x1 >= 0)
    {
        if(y1 >= 0) ret += ii[x1][y1];
        ret -= ii[x1][y2];
    }
    if(y1 >= 0) ret -= ii[x2][y1];

    return (int)ret;
}
