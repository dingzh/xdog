//
// by Jan Eric Kyprianidis <www.kyprianidis.com>
// Copyright (C) 2010-2012 Computer Graphics Systems Group at the
// Hasso-Plattner-Institut, Potsdam, Germany <www.hpi3d.de>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
#include "gpu_stgauss2.h"
#include "gpu_st.h"
#include "gpu_image.h"
// #include "gpu_sampler.h"


static texture<float, 2, cudaReadModeElementType> s_texSRC1;
static texture<float4, 2, cudaReadModeElementType> s_texSRC4;
<<<<<<< HEAD

inline __host__ __device__ texture<float,2>& texSRC1() { return s_texSRC1; }
inline __host__ __device__ texture<float4,2>& texSRC4() { return s_texSRC4; }

static texture<float4, 2, cudaReadModeElementType> s_texST;
inline __host__ __device__ texture<float4,2>& texST() { return s_texST; }

=======
static texture<float4, 2, cudaReadModeElementType> s_texST;

// inline __host__ __device__ texture<float , 2>& texSRC1() { return s_texSRC1; }
// inline __host__ __device__ texture<float4, 2>& texSRC4() { return s_texSRC4; }
// inline __host__ __device__ texture<float4, 2>& texST() { return s_texST; }


struct gpu_sampler_SRC1 {
    __host__ gpu_sampler_SRC1(const gpu_image<float>& img, cudaTextureFilterMode filter_mode=cudaFilterModePoint) {
        s_texSRC1.filterMode = filter_mode;
        GPU_SAFE_CALL(cudaBindTexture2D(0, s_texSRC1, img.ptr(), img.w(), img.h(), img.pitch()));
    }

    __host__ ~gpu_sampler_SRC1() {
        s_texSRC1.filterMode = cudaFilterModePoint;
        cudaUnbindTexture(s_texSRC1);
    }
    
    __device__ float operator()(float x, float y) const { 
        return tex2D(s_texSRC1, x, y); 
    } 
};


struct gpu_sampler_SRC4 {
    __host__ gpu_sampler_SRC4(const gpu_image<float4>& img, cudaTextureFilterMode filter_mode=cudaFilterModePoint) {
        s_texSRC4.filterMode = filter_mode;
        GPU_SAFE_CALL(cudaBindTexture2D(0, s_texSRC4, img.ptr(), img.w(), img.h(), img.pitch()));
    }

    __host__ ~gpu_sampler_SRC4() {
        s_texSRC4.filterMode = cudaFilterModePoint;
        cudaUnbindTexture(s_texSRC4);
    }
    
    __device__ float4 operator()(float x, float y) const { 
        return tex2D(s_texSRC4, x, y); 
    } 
};

struct gpu_sampler_ST {
    __host__ gpu_sampler_ST(const gpu_image<float4>& img, cudaTextureFilterMode filter_mode=cudaFilterModePoint) {
        s_texST.filterMode = filter_mode;
        GPU_SAFE_CALL(cudaBindTexture2D(0, s_texST, img.ptr(), img.w(), img.h(), img.pitch()));
    }

    __host__ ~gpu_sampler_ST() {
        s_texST.filterMode = cudaFilterModePoint;
        cudaUnbindTexture(s_texST);
    }
    
    __device__ float4 operator()(float x, float y) const { 
        return tex2D(s_texST, x, y); 
    } 
};

struct gpu_resampler {
    float2 s_;

    __host__ gpu_resampler(const gpu_image<float4>& img, float2 s, cudaTextureFilterMode filter_mode=cudaFilterModePoint) {
        s_ = s;
        s_texST.filterMode = filter_mode;
        GPU_SAFE_CALL(cudaBindTexture2D(0, s_texST, img.ptr(), img.w(), img.h(), img.pitch()));
    }

    __host__ ~gpu_resampler() {
        s_texST.filterMode = cudaFilterModePoint;
        cudaUnbindTexture(s_texST);
    }

    __device__ float4 operator()(float x, float y) const { 
        return tex2D(s_texST, s_.x * x, s_.x * y); 
    } 
};
>>>>>>>   fix compatibility to CUDA7.5

template <typename T, typename SRC>
struct stgauss2_filter {
     __device__ stgauss2_filter(const SRC& src, float sigma, float precision )
         : src_(src)
     { 
        radius_ = precision * sigma;
        twoSigma2_ = 2 * sigma * sigma;
        c_ = make_zero<T>();
        w_ = 0;
    }

    __device__ float radius() const {
        return radius_;
    }

    __device__ void operator()(float sign, float u, float2 p) {
        float k = __expf(-u * u / twoSigma2_);
        c_ += k * src_(p.x, p.y);
        w_ += k;
    }

    const SRC& src_;
    float radius_;
    float twoSigma2_;
    T c_;
    float w_;
};


template<int order, typename T, typename SRC, typename ST> 
__global__ void imp_stgauss2_filter( gpu_plm2<T> dst, SRC src, ST st, float sigma, float cos_max, 
                                     bool adaptive, float step_size, float precision ) 
{
    const int ix = __mul24(blockDim.x, blockIdx.x) + threadIdx.x;
    const int iy = __mul24(blockDim.y, blockIdx.y) + threadIdx.y;
    if(ix >= dst.w || iy >= dst.h) 
        return;

    float2 p0 = make_float2(ix + 0.5f, iy + 0.5f);
    if (adaptive) {
        float A = st2A(st(p0.x, p0.y));
        sigma *= 0.25f * (1 + A)*(1 + A);
    }
    stgauss2_filter<T,SRC> f(src, sigma, precision);
    if (order == 1) st_integrate_euler(p0, st, f, cos_max, dst.w, dst.h, step_size);
    if (order == 2) st_integrate_rk2(p0, st, f, cos_max, dst.w, dst.h, step_size);
    if (order == 4) st_integrate_rk4(p0, st, f, cos_max, dst.w, dst.h, step_size);
    dst(ix, iy) = f.c_ / f.w_;
}


gpu_image<float> gpu_stgauss2_filter( const gpu_image<float>& src, const gpu_image<float4>& st, 
                                      float sigma, float max_angle, bool adaptive,
                                      bool src_linear, bool st_linear, int order, float step_size,
                                      float precision )
{     
    if (sigma <= 0) return src;
    gpu_image<float> dst(src.size());

<<<<<<< HEAD
    gpu_sampler<float, texSRC1> src_sampler(src, src_linear? cudaFilterModeLinear : cudaFilterModePoint);
    float cos_max = cosf(radians(max_angle));

    if (src.size() == st.size()) {
        gpu_sampler<float4, texST> st_sampler(st, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
        if (order == 1) imp_stgauss2_filter<1,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 2) imp_stgauss2_filter<2,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 4) imp_stgauss2_filter<4,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
    } else {
        float2 s = make_float2((float)st.w() / src.w(), (float)st.h() / src.h());
        gpu_resampler<float4, texST> st_sampler(st, s, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
=======
    gpu_sampler_SRC1 src_sampler(src, src_linear? cudaFilterModeLinear : cudaFilterModePoint);
    float cos_max = cosf(radians(max_angle));

    if (src.size() == st.size()) {
        gpu_sampler_ST st_sampler(st, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
        if (order == 1) 
            imp_stgauss2_filter<1,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 2) 
            imp_stgauss2_filter<2,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 4) 
            imp_stgauss2_filter<4,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
    } else {
        float2 s = make_float2((float)st.w() / src.w(), (float)st.h() / src.h());
        gpu_resampler st_sampler(st, s, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
>>>>>>>   fix compatibility to CUDA7.5
        if (order == 1) imp_stgauss2_filter<1,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 2) imp_stgauss2_filter<2,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 4) imp_stgauss2_filter<4,float><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
    }
    GPU_CHECK_ERROR();
    return dst;
}


gpu_image<float4> gpu_stgauss2_filter( const gpu_image<float4>& src, const gpu_image<float4>& st, 
                                       float sigma, float max_angle, bool adaptive,
                                       bool src_linear, bool st_linear, int order, float step_size,
                                       float precision )
{     
    if (sigma <= 0) return src;
    gpu_image<float4> dst(src.size());

<<<<<<< HEAD
    gpu_sampler<float4, texSRC4> src_sampler(src, src_linear? cudaFilterModeLinear : cudaFilterModePoint);
    float cos_max = cosf(radians(max_angle));

    if (src.size() == st.size()) {
        gpu_sampler<float4, texST> st_sampler(st, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
=======
    gpu_sampler_SRC4 src_sampler(src, src_linear? cudaFilterModeLinear : cudaFilterModePoint);
    float cos_max = cosf(radians(max_angle));

    if (src.size() == st.size()) {
        gpu_sampler_ST st_sampler(st, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
>>>>>>>   fix compatibility to CUDA7.5
        if (order == 1) imp_stgauss2_filter<1,float4><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 2) imp_stgauss2_filter<2,float4><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 4) imp_stgauss2_filter<4,float4><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
    } else {
        float2 s = make_float2((float)st.w() / src.w(), (float)st.h() / src.h());
<<<<<<< HEAD
        gpu_resampler<float4, texST> st_sampler(st, s, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
=======
        gpu_resampler st_sampler(st, s, st_linear? cudaFilterModeLinear : cudaFilterModePoint);
>>>>>>>   fix compatibility to CUDA7.5
        if (order == 1) imp_stgauss2_filter<1,float4><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 2) imp_stgauss2_filter<2,float4><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
        else if (order == 4) imp_stgauss2_filter<4,float4><<<dst.blocks(), dst.threads()>>>(dst, src_sampler, st_sampler, sigma, cos_max, adaptive, step_size, precision);
    }
    GPU_CHECK_ERROR();
    return dst;
}
