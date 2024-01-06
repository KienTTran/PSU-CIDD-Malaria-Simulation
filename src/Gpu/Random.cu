//
// Created by kient on 6/17/2023.
//

#include <thrust/execution_policy.h>
#include <thrust/count.h>
#include "Random.cuh"
#include "Utils.cuh"
#include "Model.h"
#include "Core/Config/Config.h"

#if __CUDA_ARCH__ > 200
#define MAXXGRID 2147483647
#else
#define MAXXGRID 65535
#endif

#define SYNC_STREAM cudaStreamDefault

GPU::Random::Random() {
    d_states = nullptr;
    n_threads = 1024;
}

GPU::Random::~Random() {
    cudaFree(d_states);
}

__global__ void setup(int num,curandState *state, long seed)
{
    auto id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id < num) curand_init(seed, id, 0, &state[id]);
}

struct setupCurandStates
{
    long seed_;
    setupCurandStates(long seed) : seed_(seed) {}
    __device__
    curandState operator()(int id){
        curandState s;
        curand_init(seed_, id, 0, &s);
        return s;
    }
};

void GPU::Random::init_curand_states(ThrustTVectorDevice<curandState> &d_curand_states, int size, long seed){
    thrust::transform(thrust::counting_iterator<int>(0),
                      thrust::counting_iterator<int>(size),
                      d_curand_states.begin(),
                      setupCurandStates(seed));
}

void GPU::Random::init(int n, long seed) {
    cudaMalloc((void **) &d_states, sizeof(curandState) * n);
    n_blocks = (n + n_threads + 1) / n_threads;
    setup<<<n_blocks,n_threads>>>(n,d_states, seed);
    check_cuda_error(cudaDeviceSynchronize());
    check_cuda_error(cudaPeekAtLastError());

    d_location_curand_states = ThrustTVectorDevice<curandState>(Model::CONFIG->number_of_locations());
    init_curand_states(d_location_curand_states, Model::CONFIG->number_of_locations(),seed);
}

void GPU::Random::free() {
    cudaFree(d_states);
}

/*
 * Binomial are from
 * https://stackoverflow.com/questions/23561551/a-efficient-binomial-random-number-generator-code-in-java/23574723#23574723
 * to get O(Np) when p is small instead of O(N) in curand_binomial_naive
 * Not need to use Rejection because N is not very large, as mentioned in
 * https://peterchng.com/blog/2020/10/23/building-binomial-and-multinomial-samplers-in-java/
 */
// Function to generate random numbers from a binomial distribution
__device__ unsigned int curand_binomial(curandState *state, double p, unsigned int N)
{
    double log_q = log(1.0 - p);
    int count = 0;
    double sum = 0;
    while(true) {
        sum += log(curand_uniform_double(state)) / (N - count);
        if(sum < log_q) {
            return count;
        }
        count++;
    }
}

/*
 * Binomial naive from
 * https://peterchng.com/blog/2020/10/23/building-binomial-and-multinomial-samplers-in-java/
 */
__device__ unsigned int curand_binomial_naive(curandState *state, double p, unsigned int N)
{
    unsigned int count = 0;
    if(N <= 0) return 0;
    for (unsigned int i = 0; i < N; ++i)
    {
        double rand_num = curand_uniform_double(state);
        if (rand_num < p)
        {
            count++;
        }
    }
    return count;
}

/*
 * Multinomial are from gsl_multinomial
 * This is parallel multinomial on n_locations, with K draws at each location.
 * The output must have size n_locations*K
 * N is the number of trials, p is the probability of each outcome, n is the output
 */
__global__ void curand_multinomial(curandState *state, int n_locations, int K, int N[], double p[], unsigned int n[],
                                   double norm[], double sum_p[], unsigned int sum_n[])
{
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int index = thread_index; index < n_locations; index += stride){
        if(N[index] <= 0) {
            for (size_t i = 0; i < K; i++){
                n[index * K + i] = 0;
            }
            return;
        }
        for (size_t i = 0; i < K; i++)
        {
            norm[index] += p[index*K+i];
        }

        __syncthreads();

        // Calculate multinomial distribution
        for (size_t i = 0; i < K; i++)
        {
            if (p[index*K+i] > 0.0)
            {
                n[index*K+i] = curand_binomial_naive(state + index, p[index*K+i] / (norm[index] - sum_p[index]), N[index] - sum_n[index]);
            }
            else
            {
                n[index*K+i] = 0;
            }

            sum_p[index] += p[index*K+i];
            sum_n[index] += n[index*K+i];
        }
        __syncthreads();
    }
}

/*
 * Multinomial are from gsl_multinomial
 * This is parallel multinomial on n_locations, with K draws at each location.
 * The output must have size n_locations*K
 * N is the number of trials, p is the probability of each outcome, n is the output
 */
void GPU::Random::random_multinomial(int n_locations, int K, ThrustTVectorDevice<int> d_N, ThrustTVectorDevice<double> &d_p, ThrustTVectorDevice<unsigned int> &d_n){
    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int n_blocks = (n_locations + n_threads + 1) / n_threads;
    ThrustTVectorDevice<double> d_norm(n_locations,0.0);
    ThrustTVectorDevice<double> d_sum_p(n_locations,0.0);
    ThrustTVectorDevice<unsigned int> d_sum_n(n_locations,0);
    curand_multinomial<<<n_blocks, n_threads>>>(d_states,
                                                n_locations,
                                                K,
                                                thrust::raw_pointer_cast(d_N.data()),
                                                thrust::raw_pointer_cast(d_p.data()),
                                                thrust::raw_pointer_cast(d_n.data()),
                                                thrust::raw_pointer_cast(d_norm.data()),
                                                thrust::raw_pointer_cast(d_sum_p.data()),
                                                thrust::raw_pointer_cast(d_sum_n.data()));
    check_cuda_error(cudaDeviceSynchronize());
    check_cuda_error(cudaPeekAtLastError());
    d_norm.clear();
    d_sum_p.clear();
    d_sum_n.clear();
    ThrustTVectorDevice<double>().swap(d_norm);
    ThrustTVectorDevice<double>().swap(d_sum_p);
    ThrustTVectorDevice<unsigned int>().swap(d_sum_n);
}

void GPU::Random::sampling_multinomial(int n_locations, int K, ThrustTVectorDevice<int> N, ThrustTVectorDevice<double> &d_p, ThrustTVectorHost<unsigned int> &h_samples){
    ThrustTVectorDevice<unsigned int> d_n(n_locations*K,0);
    random_multinomial(n_locations, K, N, d_p, d_n);
}

/*
 * Below codes are from https://github.com/BIDData/BIDMat
 * https://github.com/BIDData/BIDMat/wiki/Random-number-generators
 * https://github.com/BIDData/BIDMat/blob/master/jni/src/Random.cu
 *
 */

__global__ void random_poisson(int n, float *A, int *B, curandState *rstates) {
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int nthreads = blockDim.x * gridDim.x;
    curandState rstate = rstates[id];
    for (int i = id; i < n; i += nthreads) {
        int cr = curand_poisson(&rstate, A[i]);
        B[i] = cr;
    }
}


__global__ void init_random(unsigned long long seed, unsigned long long offset, curandState *rstates) {
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    curand_init(seed, id, offset, &rstates[id]);
}

__forceinline__ __device__ int wait_time(curandState *prstate, float p, int n) {
    float q = - log(1-p);
    float X = 0;
    float sum = 0;
    int i = 0;
    while (i < 100 && sum <= q) {
        float E = - log(curand_uniform(prstate));  // safe since curand_uniform wont return 0
        sum += E / (n - X);
        X += 1;
        i += 1;
    }
    return X - 1;
}

int GPU::Random::randomPoisson(int n, float *A, int *B, int nthreads, unsigned long long seed, unsigned long long offset) {
    int nblocks = min(1024, max(1,nthreads/1024));
    int nth = min(n, 1024);
    curandState *rstates;
    int err;
    err = cudaMalloc(( void **)& rstates , nblocks * nth * sizeof(curandState));
    if (err > 0) {
        fprintf(stderr, "Error in cudaMalloc %d", err);
        return err;
    }
    cudaStreamSynchronize(SYNC_STREAM);
    init_random<<<nblocks,nth>>>(seed, offset, rstates);
    cudaStreamSynchronize(SYNC_STREAM);
    random_poisson<<<nblocks,nth>>>(n, A, B, rstates);
    cudaStreamSynchronize(SYNC_STREAM);
    cudaFree(rstates);
    err = cudaGetLastError();
    return err;
}


// Implements Devroye's rejection method from http://luc.devroye.org/chapter_ten.pdf

__global__ void random_binomial(int nrows, int ncols, float *A, int atype, int *C, int ctype, int *Out, curandState *rstates) {
    int nvals = nrows * ncols;
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    curandState rstate;
    float X, Y, V, p;
    int n;
    bool pflipped;
    if (id < nvals) {                            // initialize the RNGs
        rstate = rstates[id];
    }
    const float pi = 3.1415926f;
    for (int j = id; j < nvals; j += blockDim.x * gridDim.x) {
        int jcol = j / nrows;
        int jrow = j - jcol * nrows;
        switch (atype) {                           // atype and ctype determine whether these inputs are elements, rows, columns or matrices.
            case 0: p = A[0]; break;
            case 1: p = A[jrow]; break;
            case 2: p = A[jcol]; break;
            case 3: p = A[j];
        }
        switch (ctype) {
            case 0: n = C[0]; break;
            case 1: n = C[jrow]; break;
            case 2: n = C[jcol]; break;
            case 3: n = C[j];
        }
        if (p > 0.5f) {                            // flip p so that its less than 1/2.
            pflipped = true;
            p = 1.0f - p;
        } else {
            pflipped = false;
        }
        float np = n * p;
        if (np < 21) {
            X = wait_time(&rstate, p, n);           // Use a wait time method if small expected output
        } else {
            float oldp = p;
            p = floor(np) / n;                       // round np to an integral value for the rejection stage
            p = max(1e-7f, min(1 - 1e-7f, p));       // prevent divide-by-zeros
            np = n * p;
            float n1mp = n * (1-p);
            float pvar = np * (1-p);
            float delta1 = max(1.0f, floor(sqrt(pvar * log(128 * np / (81 * pi * (1-p))))));
            float delta2 = max(1.0f, floor(sqrt(pvar * log(128 * n1mp / (pi * p)))));
            float sigma1 = sqrt(pvar)*(1+delta1/(4*np));
            float sigma2 = sqrt(pvar)*(1+delta2/(4*n1mp));
            float sigma1sq = sigma1 * sigma1;
            float sigma2sq = sigma2 * sigma2;
            float c = 2 * delta1 / np;
            float a1 = 0.5f * exp(c) * sigma1 * sqrt(2*pi);
            float a2 = 0.5f * sigma2 * sqrt(2*pi);
            float a3 = exp(delta1/n1mp - delta1*delta1/(2*sigma1sq))*2*sigma1sq/delta1;
            float a4 = exp(-delta2*delta2/(2*sigma2sq))*2*sigma2sq/delta2;
            float s = a1 + a2 + a3 + a4;
            int i = 0;
            while (i < 100) {                            // Give up eventually
                i += 1;
                float U = s * curand_uniform(&rstate);
                float E1 = - log(curand_uniform(&rstate)); // safe since curand_uniform wont return 0
                if (U <= a1 + a2) {
                    float N = curand_normal(&rstate);
                    if (U <= a1) {
                        Y = sigma1 * abs(N);
                        if (Y >= delta1) continue;
                        X = floor(Y);
                        V = - E1 - N * N/2 + c;
                    } else {
                        Y = sigma2 * abs(N);
                        if (Y >= delta2) continue;
                        X = floor(-Y);
                        V = - E1 - N * N/2;
                    }
                } else {
                    float E2 = - log(curand_uniform(&rstate));
                    if (U <= a1 + a2 + a3) {
                        Y = delta1 + 2*sigma1sq*E1/delta1;
                        X = floor(Y);
                        V = - E2 - delta1*Y/(2*sigma1sq) + delta1/n1mp;
                    } else {
                        Y = delta2 + 2*sigma2sq*E1/delta2;
                        X = floor(-Y);
                        V = - E2 - delta2*Y/(2*sigma2sq);
                    }
                }
                if (X < - np || X > n1mp) continue;
                if (V > lgamma(np+1) + lgamma(n1mp+1) - lgamma(np+X+1) - lgamma(n1mp-X+1) + X*log(p/(1-p))) continue;
                break;
            }
            X += np;
            X += wait_time(&rstate, (oldp-p)/(1-p), n-X); // Now correct for rounding np to integer
        }
        if (pflipped) {                                  // correct for flipped p.
            X = n - X;
        }
        Out[j] = (int)X;
    }
}

int GPU::Random::randomBinomial(int nrows, int ncols, float *A, int atype, int *C, int ctype, int *Out, unsigned long long seed, unsigned long long offset) {
    int nvals = nrows * ncols;
    int nthreads = min(256, 32*(1+(nvals-1)/32));  // at least nvals, round up to a multiple of 32 l.e. 256
    int nblocks = min(128, 1 + (nvals-1)/nthreads);
    curandState *rstates;
    int err = cudaMalloc(( void **)& rstates , nthreads * nblocks * sizeof(curandState));
    if (err > 0) {
        fprintf(stderr, "Error in cudaMalloc %d", err);
        return err;
    }
    cudaStreamSynchronize(SYNC_STREAM);
    init_random<<<nblocks,nthreads>>>(seed, offset, rstates);
    cudaStreamSynchronize(SYNC_STREAM);
    random_binomial<<<nblocks,nthreads>>>(nrows, ncols, A, atype,  C, ctype, Out, rstates);
    cudaStreamSynchronize(SYNC_STREAM);
    cudaFree(rstates);
    err = cudaGetLastError();
    return err;
}

//
// Based on  Marsaglia and Tsang, <i>A Simple Method for Generating Gamma Variables.</i>
// ACM Transactions on Mathematical Software,  Volume 26 Issue 3, September, 2000
//
// Transformation for a < 1 is based on Stuart's theorem, section IV.6.4 of Devroye's book:
// http://luc.devroye.org/rnbookindex.html
//
__forceinline__ __device__ float random_gamma1(float a, curandState *prstate) {
    int small_int = 0;
    if (a < 1.0f) {
        small_int = 1;
        a += 1;
    }
    float x = 0;
    float d = a - 0.33333333f;
    float c = 0.33333333f / sqrt(d);
    while (1) {
        float z = curand_normal(prstate);
        float v0 = 1 + c * z;
        if (v0 <= 0) continue;
        float u = curand_uniform(prstate);
        float v = v0*v0*v0;
        x = d * v;
        float z2 = z * z;
        if (u < 1 - 0.0331f*z2*z2) break;
        if (log(u) < 0.5f*z2 + d - x + d*log(v)) break;
    }
    if (small_int) {
        a -= 1;
        float u = curand_uniform(prstate);
        x *= pow(u, 1/a);
    }
    return x;
}

__global__ void random_gamma(int nrows, int ncols, float *A, int atype, float *B, int btype, float *Out, curandState *rstates) {
    int nvals = nrows * ncols;
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    float a, b, x;
    curandState *rstate;
    if (id < nvals) {                            // initialize the RNGs
        rstate = &rstates[id];
    }
    x = 0;
    for (int j = id; j < nvals; j += blockDim.x * gridDim.x) {
        int jcol = j / nrows;
        int jrow = j - jcol * nrows;
        switch (atype) {                           // atype and ctype determine whether these inputs are elements, rows, columns or matrices.
            case 0: a = A[0]; break;
            case 1: a = A[jrow]; break;
            case 2: a = A[jcol]; break;
            case 3: a = A[j];
        }
        switch (btype) {
            case 0: b = B[0]; break;
            case 1: b = B[jrow]; break;
            case 2: b = B[jcol]; break;
            case 3: b = B[j];
        }
        x = random_gamma1(a, rstate);
        Out[j] = x * b;
    }
}

int GPU::Random::randomGamma(int nrows, int ncols, float *A, int atype, float *B, int btype, float *Out, unsigned long long seed, unsigned long long offset) {
    int nvals = nrows * ncols;
    int nthreads = min(256, 32*(1+(nvals-1)/32));  // at least nvals, round up to a multiple of 32 l.e. 256
    int nblocks = min(128, 1 + (nvals-1)/nthreads/256);
    curandState *rstates;
    int err = cudaMalloc(( void **)& rstates , nthreads * nblocks * sizeof(curandState));
    if (err > 0) {
        fprintf(stderr, "Error in cudaMalloc %d", err);
        return err;
    }
    cudaStreamSynchronize(SYNC_STREAM);
    init_random<<<nblocks,nthreads>>>(seed, offset, rstates);
    cudaStreamSynchronize(SYNC_STREAM);
    random_gamma<<<nblocks,nthreads>>>(nrows, ncols, A, atype,  B, btype, Out, rstates);
    cudaStreamSynchronize(SYNC_STREAM);
    cudaFree(rstates);
    err = cudaGetLastError();
    return err;
}

