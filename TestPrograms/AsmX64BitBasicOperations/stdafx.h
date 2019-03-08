// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>

#include <intrin.h>

// TODO: reference additional headers your program requires here

#define	X64_API	__declspec(dllexport)

typedef signed __int64 sint64;
typedef unsigned __int64 uint64;
typedef unsigned char uint8;

const uint64 MaxValueUint64 = 0xFFFFFFFFFFFFFFFFULL;
const uint64 Bit63 = 0x8000000000000000ULL;

static uint64 zero[2] = { 0, 0 };
static int KaratsubaBreakPoint = 18;
static int FourierBreakPoint = 6 * 1024 + 512;
static int CarrylessKaratsubaBreakPoint = 14;
static int CarrylessFourierBreakPoint = 14 * 1024;
