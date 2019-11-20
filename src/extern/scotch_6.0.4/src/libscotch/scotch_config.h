#pragma once

#ifdef SCOTCH_EXPORTS
    #define SCOTCH_EXPORT _declspec(dllexport)
#else
    #define SCOTCH_EXPORT _declspec(dllimport)
#endif

#define restrict 
#define PETSC_USE_MICROSOFT_TIME 1