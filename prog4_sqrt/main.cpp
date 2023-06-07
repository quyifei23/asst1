#include <stdio.h>
#include <algorithm>
#include <pthread.h>
#include <math.h>
#include <immintrin.h>
#include<avx2intrin.h>
#include "CycleTimer.h"
#include "sqrt_ispc.h"

using namespace ispc;

extern void sqrtSerial(int N, float startGuess, float* values, float* output);

static void verifyResult(int N, float* result, float* gold) {
    for (int i=0; i<N; i++) {
        if (fabs(result[i] - gold[i]) > 1e-4) {
            printf("Error: [%d] Got %f expected %f\n", i, result[i], gold[i]);
        }
    }
}

void sqrtAVX2(int N,float startGuess,float* values,float* output)
{
    //__m256 one=_mm256_set1_ps(1.f);
    __m256 Ofive=_mm256_set1_ps(0.5f);
    __m256 three=_mm256_set1_ps(3.f);
    __m256 minierror=_mm256_set1_ps(0.9999f);
    __m256 maxerror=_mm256_set1_ps(1.0001f);
    __m256 guess,tmp,mask1,mask2,mask,value;
    float judge[8];
    for(int i=0;i<N;i+=8)
    {
        guess=_mm256_set1_ps(startGuess);
        value=_mm256_loadu_ps(values+i);
        tmp=_mm256_mul_ps(guess,value);
        tmp=_mm256_mul_ps(guess,tmp);
        mask1=_mm256_cmp_ps(tmp,maxerror,_CMP_GT_OS);
        mask2=_mm256_cmp_ps(tmp,minierror,_CMP_LT_OS);
        mask=_mm256_or_ps(mask1,mask2);
        _mm256_storeu_ps(judge,mask);
        while(judge[0]||judge[1]||judge[2]||judge[3]||judge[4]||judge[5]||judge[6]||judge[7])
        {
            tmp=_mm256_sub_ps(three,tmp);
            guess=_mm256_mul_ps(tmp,guess);
            guess=_mm256_mul_ps(Ofive,guess);

            tmp=_mm256_mul_ps(guess,value);
            tmp=_mm256_mul_ps(guess,tmp);
            mask1=_mm256_cmp_ps(tmp,maxerror,_CMP_GT_OS);
            mask2=_mm256_cmp_ps(tmp,minierror,_CMP_LT_OS);
            mask=_mm256_or_ps(mask1,mask2);
            _mm256_storeu_ps(judge,mask);
        }
        value=_mm256_mul_ps(guess,value);
        _mm256_storeu_ps(output+i,value);
    }
}
int main() {

    const unsigned int N = 20 * 1000 * 1000;
    const float initialGuess = 1.0f;

    float* values = new float[N];
    float* output = new float[N];
    float* gold = new float[N];

    for (unsigned int i=0; i<N; i++)
    {
        // TODO: CS149 students.  Attempt to change the values in the
        // array here to meet the instructions in the handout: we want
        // to you generate best and worse-case speedups
        
        // starter code populates array with random input values
        values[i] = .001f + 2.998f * static_cast<float>(rand()) / RAND_MAX;
        
        //maxmize the speedup
        //values[i]=2.99f;

        //minimize the speedup
        //if(i%8)
        //values[i]=1.f;
        //else
        //values[i]=2.99f;
    }

    // generate a gold version to check results
    for (unsigned int i=0; i<N; i++)
        gold[i] = sqrt(values[i]);

    //
    // And run the serial implementation 3 times, again reporting the
    // minimum time.
    //
    double minSerial = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrtSerial(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minSerial = std::min(minSerial, endTime - startTime);
    }

    printf("[sqrt serial]:\t\t[%.3f] ms\n", minSerial * 1000);

    verifyResult(N, output, gold);


    double minAVX2 = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrtAVX2(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minAVX2 = std::min(minAVX2, endTime - startTime);
    }

    printf("[sqrt AVX2]:\t\t[%.3f] ms\n", minAVX2 * 1000);

    verifyResult(N, output, gold);
    //
    // Compute the image using the ispc implementation; report the minimum
    // time of three runs.
    //
    double minISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrt_ispc(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minISPC = std::min(minISPC, endTime - startTime);
    }

    printf("[sqrt ispc]:\t\t[%.3f] ms\n", minISPC * 1000);

    verifyResult(N, output, gold);

    // Clear out the buffer
    for (unsigned int i = 0; i < N; ++i)
        output[i] = 0;

    //
    // Tasking version of the ISPC code
    //
    double minTaskISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrt_ispc_withtasks(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minTaskISPC = std::min(minTaskISPC, endTime - startTime);
    }

    printf("[sqrt task ispc]:\t[%.3f] ms\n", minTaskISPC * 1000);

    verifyResult(N, output, gold);
    
    
    printf("\t\t\t\t(%.2fx speedup from AVX2)\n", minSerial/minAVX2);
    printf("\t\t\t\t(%.2fx speedup from ISPC)\n", minSerial/minISPC);
    printf("\t\t\t\t(%.2fx speedup from task ISPC)\n", minSerial/minTaskISPC);

    delete [] values;
    delete [] output;
    delete [] gold;

    return 0;
}
