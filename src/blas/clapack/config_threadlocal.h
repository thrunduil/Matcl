#pragma once

#ifdef _MSC_VER
#define CLAPACK_THREADLOCAL __declspec(thread)
#else
#define CLAPACK_THREADLOCAL __thread
#endif
