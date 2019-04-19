/**
 * @file n_body.hpp
 * @author Salvatore Cardamone
 * @brief Basic N-Body implementation using CUDA.
 */
#ifndef N_BODY_CUDA_HPP__
#define N_BODY_CUDA_HPP__

// See https://stackoverflow.com/a/6978720
#ifdef __CUDACC__
#define CUDA_MEMBER __host__ __device__
#else
#define CUDA_MEMBER
#endif /* #ifdef __CUDACC__ */


class NBody {

public:

private:
  unsigned int nParticles_;
  std::vector<float4> particle_, acc_;

  CUDA_MEMBER void PairwiseInteraction( const unsigned int& i, const unsigned int& j );

};

#endif /* #ifndef N_BODY_CUDA_HPP__ */
