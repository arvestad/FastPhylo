
#ifndef STD_C_UTILS_H
#define STD_C_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#define MAX(A,B) ((A)<(B) ? (B) : (A))
#define MIN(A,B) ((A)<(B) ? (A) : (B))


int randomInt(int low, int high);
float randomFloat(float low, float high);


#ifdef __cplusplus
}
#endif

#endif
