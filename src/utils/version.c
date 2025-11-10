/**
 * version.c - Version information
 */

#include "fdpricing.h"

uint32_t fdp_get_version(void)
{
    return FDP_VERSION;
}

int fdp_is_compatible(void)
{
    uint32_t version = fdp_get_version();
    uint32_t major = (version >> 16) & 0xFF;
    return major == FDP_VERSION_MAJOR;
}
