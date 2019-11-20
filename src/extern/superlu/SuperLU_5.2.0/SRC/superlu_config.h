#pragma once

#ifdef SUPERLU_EXPORTS
    #define SUPERLU_EXPORT __declspec(dllexport)
#else
    #define SUPERLU_EXPORT __declspec(dllimport)
#endif

void matcl_report_error(const char* msg);
#define USER_ABORT(msg) matcl_report_error(msg);