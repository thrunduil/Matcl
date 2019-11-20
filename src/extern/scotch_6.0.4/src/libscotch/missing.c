
#include <string.h>

int strcasecmp(const char *s1, const char *s2)
{
    return _stricmp(s1, s2);
};

int strncasecmp(const char *f1, const char *f2, size_t n )
{
    return _strnicmp(f1,f2,n);
};