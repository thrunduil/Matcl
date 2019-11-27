#pragma once

#ifdef SUPERLU_EXPORTS
    #define SUPERLU_EXPORT __declspec(dllexport)
#else
    #define SUPERLU_EXPORT __declspec(dllimport)
#endif