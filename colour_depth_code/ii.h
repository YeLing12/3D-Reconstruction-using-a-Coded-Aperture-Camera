#ifndef __II_H__
#define __II_H__

#ifdef __cplusplus
extern "C" {
#endif

void generate_ii(int w, int h, int **in, ullong **out);
int sample_ii(int w, int h, ullong **ii, int x1, int y1, int x2, int y2);

#ifdef __cplusplus
}
#endif

#endif
