#include <stdio.h>
#include <stdlib.h>

void xermsg_(const char *librar, const char *subrou, const char *messg,
             const int *nerr, const int *level,
             int librar_len, int subrou_len, int messg_len)
{
    if (*nerr < -9999999 || *nerr > 99999999 || *nerr == 0 ||
        *level < -1 || *level > 2) {
	fprintf(stderr, "xermsg: invalid error number or level\n");
        abort();
    }
    const char *level_str = "fatal";
    if (*level <= 0) {
        level_str = "warning";
    } else if (*level == 1) {
        level_str = "error";
    }
    fprintf(stderr, "%.*s.%.*s: [%s] %.*s (%i)\n",
            librar_len, librar, subrou_len, subrou, level_str,
            messg_len, messg, nerr);
    fflush(stderr);

    /* TODO: allow user to trap non-fatal errors */
    if (*level > 1) {
        abort();
    }
}
