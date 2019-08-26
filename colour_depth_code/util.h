#ifndef __UTIL_H__
#define __UTIL_H__

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


typedef unsigned char uchar;
typedef unsigned long long ullong;


#ifdef __cplusplus
extern "C" {
#endif

#define MyMalloc(type)   type   *malloc_ ## type ## 1d(int n);                \
                         type  **malloc_ ## type ## 2d(int w, int h);         \
                         type ***malloc_ ## type ## 3d(int w, int h, int d);  \
                         void free_ ## type ## 1d(type   *buf);               \
                         void free_ ## type ## 2d(type  **buf);               \
                         void free_ ## type ## 3d(type ***buf);

MyMalloc(uchar)
MyMalloc(int)
MyMalloc(float)
MyMalloc(double)
MyMalloc(ullong)
#undef MyMalloc


#ifdef __cplusplus
}
#endif

#endif
