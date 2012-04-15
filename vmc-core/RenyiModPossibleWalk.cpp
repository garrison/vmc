#include <cmath>

#include "RenyiModPossibleWalk.hpp"

probability_t RenyiModPossibleWalk::probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const
{
    (void) phibeta1_ratio;
    (void) phibeta2_ratio;
    return std::norm(phialpha1_ratio * phialpha2_ratio);
}
