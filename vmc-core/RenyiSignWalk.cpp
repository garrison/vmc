#include <cmath>

#include "RenyiSignWalk.hpp"

probability_t RenyiSignWalk::probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const
{
    return std::abs(phialpha1_ratio * phialpha2_ratio * phibeta1_ratio * phibeta2_ratio);
}
