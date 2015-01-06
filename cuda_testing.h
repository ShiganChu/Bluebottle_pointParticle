/****h* Bluebottle/cuda_testing_kernel
 * NAME
 *  cuda_testing_kernel
 * FUNCTION
 *  Bluebottle testing CUDA kernel functions.
 ******
 */

#ifndef _CUDA_TESTING_H
#define _CUDA_TESTING_H

extern "C"
{
#include "bluebottle.h"
}

/****f* cuda_testing_kernel/PP_memcpy_p_test<<<>>>()
 * NAME
 *  PP_memcpy_p_test<<<>>>()
 * USAGE
 */
__global__ void PP_memcpy_p_test(real *dst, real *src, dom_struct *dom);
/*
 * FUNCTION
 *  Copy the computed _rhs_p to _p for testing output.
 * ARGUMENTS
 *  dst -- the destination _p
 *  src -- the source _rhs_p
 *  dom -- the dom_struct associated with this device
 ******
 */

#endif
