/* get GCC version info */
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

/* check for GCC 4.1.2 (40102) */
#if GCC_VERSION < 40404
#define FORT_MOD_NAME(NAME) __p3dfft__##NAME
#define FORTNAME(NAME) NAME##_

/* check for GCC 4.4.4 (40404) */
#elif GCC_VERSION >= 40404
#define FORT_MOD_NAME(NAME) __p3dfft_MOD_##NAME
#define FORTNAME(NAME) NAME##_
#endif
