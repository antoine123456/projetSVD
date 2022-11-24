/* Contains few benchmark utility functions to avoid too many repetitions */

#pragma once

#include <stdio.h>
#include <math.h>

// Benchmark return value structure
typedef struct {
    double mean;
    double std;
} benchmark_t;

#define NSAMPLES 10

// Computes b's mean and std
void Stats(benchmark_t *b, double times[NSAMPLES]);

// Benchmark display (terminal)
void PrintBenchmark(benchmark_t b, int n);