#include <sys/resource.h>
#include <stddef.h>

// Subroutine returns the peak resident set size in bytes
//    On linux systems we need to multiply by 1024 to convert from kB.
//    On Max systems, bytes are automatically returned.

void get_max_rss(long *max_rss_bytes) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        *max_rss_bytes = usage.ru_maxrss * 1024;
    } else {
        *max_rss_bytes = -1; // Error indicator
    }
}
