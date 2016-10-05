// cuda kernel to calculate the fresnel propagation kernel.
// Johannes Hagemann 2011 - 2016
#include "/z/johannes/SciPAL/include/base/CudaComplex.h"
#
__global__
void __generate_k_element(double2* out,
                          double2* in,
                          double F,
                          size_t const height,
                          size_t const width,
                          size_t const size)

{
    typedef double NumberType;
    typedef typename SciPAL::CudaComplex<NumberType> Complex;

    //Calculate the thread ID. The thread ID determines which pixel is calculated.
    size_t index = blockDim.x*blockIdx.x+threadIdx.x;

    //Prevents kernel to calculate something outside the image vector.
    if(index<size)
    {
    size_t x = index % width - (width/2);
    size_t y = index / width - (height/2);

    NumberType kx2 = x * x * 4*M_PI*M_PI / (width * width);
    NumberType ky2 = y * y * 4*M_PI*M_PI / (height * height);
    NumberType v = - 0.5 * (kx2+ky2) / (2 * M_PI * F);
    NumberType cosv;
    NumberType sinv;
    sincos(v, &sinv, &cosv);
    Complex e(cosv, sinv); //no scaling factor, matlab fft take care of that
    //we have our value, now find position in shifted array
    size_t new_x = (x + width)%width;
    size_t new_y = (y + height)%height;
    size_t new_index = new_y * width + new_x;

    out[new_index] = toNumberType2(e * Complex(in[new_index]));
    }
}
