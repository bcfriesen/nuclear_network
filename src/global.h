#ifndef GLOBAL
#define GLOBAL

#include <stddef.h>
#ifdef MAIN_FILE
// number of isotopes and also the number of ODEs to solve
const size_t nvar = 13;
#else
extern const size_t nvar;
#endif

#endif
