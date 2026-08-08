
#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#define MAX(A, B) ((A) < (B) ? (B) : (A))
#define MIN(A, B) ((A) < (B) ? (A) : (B))

    int randomInt(int low, int high);
    float randomFloat(float low, float high);

#ifdef __cplusplus
}
#endif
