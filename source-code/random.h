#include <random>
#include "SFMT.h"

using namespace std;

sfmt_t sfmt;

void SFMTinitialize() {
    sfmt_init_gen_rand(&sfmt, unsigned(time(0)));
}

unsigned int drand() {
    return sfmt_genrand_uint64(&sfmt);
}

double prand() {
    return sfmt_genrand_res53(&sfmt);
}