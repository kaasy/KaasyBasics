#include "stdafx.h"
#include <math.h>

extern "C"
{
	X64_API uint8 __cdecl Add(uint64* input1, uint64* input2, uint64* result, uint8 carry, int n);
	X64_API uint8 __cdecl Subtract(uint64* input1, uint64* input2, uint64* result, uint8 borrow, int n);
	X64_API void __cdecl inverseFFTCore(uint64* number, int startIndex, int coreLength, uint64* wpower);
	X64_API void __cdecl butterflyCore(uint64* number, int startIndex, int coreLength);
	X64_API void __cdecl forwardFFTCore(uint64* number, int startIndex, int coreLength, uint64* wpower);
	X64_API void __cdecl submod(uint64* a, uint64* b);
	X64_API void __cdecl mulmod(uint64* a, uint64* b);
	X64_API void __cdecl PowMod(uint64 *a, uint64* power);
	X64_API void __cdecl Multiply(uint64* input1, uint64* input2, uint64* result, int n);
	X64_API void __cdecl Square(uint64* input, uint64* result, int n);
	X64_API void __cdecl Karatsuba(uint64* input1, uint64* input2, uint64* result, int n, uint64* temporaryBuffer = NULL);
	X64_API void __cdecl KaratsubaSquare(uint64* input, uint64* result, int n, uint64* temporaryBuffer = NULL);
	X64_API void __cdecl PropagateCarries(uint64* number, int bitsPerDigit, int numberLength);
	X64_API uint64 __cdecl MultiplyDigitAndAdd(uint64* input, uint64 digit, uint64* result, int n);
	X64_API uint64 __cdecl MultiplyDigitAndSubtract(uint64* input, uint64 digit, uint64* result, int n);
	X64_API uint64 __cdecl SubShl(uint64* inputOutput, uint64 *input2, uint8 shift, int n);
	X64_API	int GetKaratsubaCarrylessMultiplicationBufferSize(int n);
}

static bool FFTRootsInitialized = false;
static uint64 RootMultipliers[53 * 2];
static uint64 RootInverters[53 * 2];
const uint64 FFTPrimeModuloHigh = 0xFFFFFFFFFFFFFFFFULL;
const uint64 FFTPrimeModuloLow = 0xFFC0000000000001ULL;
const uint64 InverseLowP = 0x0040000000000001ULL;	// == FFTPrimeModuloLow ^ -1 mod 2^64
const uint64 NegativeInverseLowP = 0xFFBFFFFFFFFFFFFFULL;	// == -FFTPrimeModuloLow ^ -1 mod 2^64
static uint64 FFTPrimeModuloDigits[2] = { FFTPrimeModuloLow, FFTPrimeModuloHigh };


static void inline remodP(uint64 *value)
{
	if (value[1] == FFTPrimeModuloHigh && value[0] >= FFTPrimeModuloLow)
	{
		value[0] -= FFTPrimeModuloLow;
		value[1] = 0;
	}
}

static void InitializeFFTRoots()
{
	if (FFTRootsInitialized)
	{
		return;
	}

	uint64 roots[54 * 2];
	uint64 invertedRoots[54 * 2];
	uint64 pminus2[2], temporary[2];
	pminus2[1] = FFTPrimeModuloHigh;
	pminus2[0] = FFTPrimeModuloLow - 2ULL;
	invertedRoots[53 * 2 + 0] = roots[53 * 2 + 0] = 0x4087ea9b3d2c84a8ULL;
	invertedRoots[53 * 2 + 1] = roots[53 * 2 + 1] = 0xfaaf0d013359bbb3ULL;
	PowMod(&invertedRoots[53 * 2], pminus2);

	for (int i = 53 * 2; (i -= 2) >= 0; )
	{
		roots[i] = roots[i + 2];
		roots[i + 1] = roots[i + 3];
		mulmod(&roots[i], &roots[i]);

		invertedRoots[i] = invertedRoots[i + 2];
		invertedRoots[i + 1] = invertedRoots[i + 3];
		mulmod(&invertedRoots[i], &invertedRoots[i]);

		temporary[0] = roots[i];
		temporary[1] = roots[i + 1];
		mulmod(temporary, &roots[i + 2]);
		RootMultipliers[i] = FFTPrimeModuloLow;
		RootMultipliers[i + 1] = FFTPrimeModuloHigh;
		submod(&RootMultipliers[i], temporary);
		remodP(&RootMultipliers[i]);

		temporary[0] = invertedRoots[i];
		temporary[1] = invertedRoots[i + 1];
		mulmod(temporary, &invertedRoots[i + 2]);
		RootInverters[i] = FFTPrimeModuloLow;
		RootInverters[i + 1] = FFTPrimeModuloHigh;
		submod(&RootInverters[i], temporary);
		remodP(&RootInverters[i]);
	}

	FFTRootsInitialized = true;
}

static void fft(uint64 *list, int n)
{
	InitializeFFTRoots();

	uint64 wpower[2];
	for (int k = n; (k >>= 1) > 0;)
	{
		butterflyCore(list, 0, k);
		wpower[0] = 1ULL;
		wpower[1] = 0ULL;
		for (int j = k << 1, x = 0; j < n; x++, j += k << 1)
		{
			int pos = 0, bit = 1;
			while ((x & bit) != 0)
			{
				bit <<= 1;
				pos++;
			}
			mulmod(wpower, &RootMultipliers[pos << 1]);
			forwardFFTCore(list, j, k, wpower);
		}
	}
}

static void ifft(uint64* list, int n)
{
	InitializeFFTRoots();

	uint64 wpower[2];
	for (int k = 1; k < n; k <<= 1)
	{
		butterflyCore(list, 0, k);
		wpower[0] = 1ULL;
		wpower[1] = 0ULL;
		for (int j = k << 1, x = 0; j < n; x++, j += k << 1)
		{
			int pos = 0, bit = 1;
			while ((x & bit) != 0)
			{
				bit <<= 1;
				pos++;
			}
			mulmod(wpower, &RootInverters[pos << 1]);
			inverseFFTCore(list, j, k, wpower);
		}
	}
}

static int minimumInt(int a, int b)
{
	return a <= b ? a : b;
}

const double LOG2_E = 1.4426950408889634073599246810019;    // = log2(e) = 1 / ln(2)
const double LOG2_3 = 1.5849625007211561814537389439478;   // == log2(3) == ln(3) / ln(2)

static int CompareBaseBits(int bitsPerDigit, uint64 n)
{   //middle digit <= n * (base - 1) ^ 2 + [2^(128 - bitsPerDigit)] <= P - 1
	uint64 value[1];
	uint64 square[2];
	uint64 multiplier[2];
	uint64 result[4];
	uint64 adder[4];
	value[0] = (1ULL << bitsPerDigit) - 1;
	multiplier[0] = n; multiplier[1] = 0ULL;
	Square(value, square, 1);
	Multiply(square, multiplier, result, 2);
	adder[0] = adder[1] = adder[2] = adder[3] = 0ULL;
	if (128 - bitsPerDigit < 64)
	{
		adder[0] = 1ULL << (128 - bitsPerDigit);
	}
	else
	{
		adder[1] = 1ULL << (64 - bitsPerDigit);
	}
	Add(result, adder, result, 0, 4);

	if (result[2] > 0 || result[3] > 0 || (result[0] >= FFTPrimeModuloLow && result[1] == FFTPrimeModuloHigh))
	{
		return 1;
	}
	if (result[0] == FFTPrimeModuloLow - 1ULL && result[1] == FFTPrimeModuloHigh)
	{
		return 0;
	}
	return -1;
}


static void GetFFTParameters(int n, int &bitsPerDigit, int &logFFTSize)
{
	uint64 inputBits = n * 64ULL;
	bitsPerDigit = 63;
	uint64 fftSize = (inputBits + bitsPerDigit - 1) / bitsPerDigit;
	bitsPerDigit = minimumInt(63, (int)ceil(64 - log((double)fftSize) * (LOG2_E * 0.5)));  //use bitsPerDigit <= 63
	while (true)
	{
		fftSize = (inputBits + bitsPerDigit - 1) / bitsPerDigit;
		if (CompareBaseBits(bitsPerDigit, fftSize) <= 0)
		{
			logFFTSize = 1 + (int)ceil(log((double)fftSize) * LOG2_E);
			return;
		}
		bitsPerDigit--;
	}
}

static __int64 minimumInt64(__int64 a, __int64 b)
{
	return a <= b ? a : b;
}

static uint64 ReadBits(uint64 *number, int numberLength, __int64 bitPosition, int count)
{
	__int64 index = bitPosition >> 6;
	int shift = (int)bitPosition & 63;
	uint64 result = number[index] >> shift;
	shift = 64 - shift;
	if (shift < count && index + 1 < numberLength)
	{
		result |= number[index + 1] << shift;
	}
	result &= (1ULL << count) - 1ULL;
	return result;
}

static void WriteBits(uint64 *number, int numberLength, __int64 bitPosition, int count, uint64 value)
{
	__int64 index = bitPosition >> 6;
	int shift = (int)bitPosition & 63;
	uint64 mask = (1ULL << count) - 1ULL;
	value &= mask;
	number[index] = (number[index] & (~(mask << shift))) | (value << shift);
	shift = 64 - shift;
	if (shift < count && index + 1 < numberLength)
	{
		number[index + 1] = (number[index + 1] & (~(mask >> shift))) | (value >> shift);
	}
}

static void GroupDigits(uint64 *number, int numberLength, int bitsPerDigit, uint64 *result, int resultLength)
{
	__int64 numberOfDigits = minimumInt64((__int64)resultLength, ((__int64)numberLength * 64LL + bitsPerDigit - 1) / bitsPerDigit);
	__int64 position = numberOfDigits * bitsPerDigit;
	for (__int64 i = numberOfDigits << 1; (i -= 2) >= 0;)
	{
		position -= bitsPerDigit;
		result[i] = ReadBits(number, numberLength, position, bitsPerDigit);
		result[i + 1] = 0ULL;
	}
	if (resultLength > numberOfDigits)
	{
		memset(&result[numberOfDigits << 1], 0, (resultLength - numberOfDigits) << 4);
	}
}

static void UngroupDigits(uint64* number, int numberLength, int bitsPerDigit, uint64 *result, int resultLength)
{
	__int64 numberOfDigits = minimumInt64((__int64)numberLength, ((__int64)resultLength * 64LL + bitsPerDigit - 1) / bitsPerDigit);
	__int64 position = numberOfDigits * bitsPerDigit;
	for (__int64 i = numberOfDigits << 1; (i -= 2) >= 0;)
	{
		position -= bitsPerDigit;
		WriteBits(result, resultLength, position, bitsPerDigit, number[i]);
	}
}

static void shrmod(uint64* x, int shift)
{
	uint64 number[2] = { x[0], x[1] };
	while (shift >= 64)
	{
		uint64 carry = MultiplyDigitAndAdd(FFTPrimeModuloDigits, NegativeInverseLowP * number[0], number, 2);
		number[0] = number[1];
		number[1] = carry;
		shift -= 64;
	}
	if (shift != 0)
	{
		uint64 mask = (1ULL << shift) - 1ULL;
		uint64 carry = MultiplyDigitAndAdd(FFTPrimeModuloDigits, NegativeInverseLowP * number[0] & mask, number, 2);
		x[0] = (number[0] >> shift) | (number[1] << (64 - shift));
		x[1] = (number[1] >> shift) | (carry << (64 - shift));
		return;
	}
	x[0] = number[0];
	x[1] = number[1];
}

static void fourierMultiplication(uint64* input1, uint64* input2, uint64* result, int n)
{
	int fftBitsPerDigit, fftLog2Size;
	GetFFTParameters(n, fftBitsPerDigit, fftLog2Size);
	int fftN = 1 << fftLog2Size;

	uint64 *memory = new uint64[fftN * 4];
	uint64 *number1 = &memory[0];
	uint64 *number2 = &memory[fftN * 2];
	GroupDigits(input1, n, fftBitsPerDigit, number1, fftN);
	GroupDigits(input2, n, fftBitsPerDigit, number2, fftN);

	fft(number1, fftN);
	fft(number2, fftN);
	for (int i = fftN * 2; (i -= 2) >= 0;)
	{
		mulmod(&number1[i], &number2[i]);
		shrmod(&number1[i], fftLog2Size);
	}
	ifft(number1, fftN);

	PropagateCarries(number1, fftBitsPerDigit, fftN);
	UngroupDigits(number1, fftN, fftBitsPerDigit, result, n << 1);
	delete[] memory;
}

static void fourierSquare(uint64* input, uint64* result, int n)
{
	int fftBitsPerDigit, fftLog2Size;
	GetFFTParameters(n, fftBitsPerDigit, fftLog2Size);
	int fftN = 1 << fftLog2Size;

	uint64 *number1 = new uint64[fftN * 2];
	GroupDigits(input, n, fftBitsPerDigit, number1, fftN);

	fft(number1, fftN);
	for (int i = fftN * 2; (i -= 2) >= 0;)
	{
		mulmod(&number1[i], &number1[i]);
		shrmod(&number1[i], fftLog2Size);
	}
	ifft(number1, fftN);

	PropagateCarries(number1, fftBitsPerDigit, fftN);
	UngroupDigits(number1, fftN, fftBitsPerDigit, result, n << 1);
	delete[] number1;
}

static void fastestMultiplication(uint64* input1, uint64* input2, uint64* result, int n)
{
	if (n < FourierBreakPoint)
	{
		Karatsuba(input1, input2, result, n);
		return;
	}
	fourierMultiplication(input1, input2, result, n);
}

static void fastestSquare(uint64* input, uint64* result, int n)
{
	if (n < FourierBreakPoint)
	{
		KaratsubaSquare(input, result, n);
		return;
	}
	fourierSquare(input, result, n);
}

static uint64 carrylessMultipyAndXor(uint64* input, uint64 digit, uint64* result, int n)
{
	__m128i reg1 = _mm_set1_epi64x((__int64)digit);
	__m128i carry = _mm_setzero_si128();
	for (int i = 0; i < n; i++)
	{
		__m128i reg2 = _mm_loadl_epi64((__m128i*)(&input[i]));
		carry = _mm_xor_si128(carry, _mm_loadl_epi64((__m128i*)(&result[i])));
		carry = _mm_xor_si128(carry, _mm_clmulepi64_si128(reg1, reg2, 0x00));
		_mm_storel_epi64((__m128i*)(&result[i]), carry);
		carry = _mm_srli_si128(carry, 8);
	}
	return (uint64)_mm_extract_epi64(carry, 0);
}

static void carrylessMultiplication(uint64* input1, uint64* input2, uint64* result, int n)
{
	for (int i = n * 2; --i >= 0; )
	{
		result[i] = 0ULL;
	}
	int m = n & -2;
	for (int i = m; (i -= 2) >= 0; )
	{
		__m128i reg1 = _mm_loadu_si128((__m128i*)(&input1[i]));
		__m128i cache1 = _mm_xor_si128(reg1, _mm_srli_si128(reg1, 8));
		for (int j = m; (j -= 2) >= 0; )
		{
			__m128i resultLow = _mm_loadu_si128((__m128i*)(&result[i + j]));
			__m128i resultHigh = _mm_loadu_si128((__m128i*)(&result[i + j + 2]));
			__m128i reg2 = _mm_loadu_si128((__m128i*)(&input2[j]));
			__m128i cache2 = _mm_xor_si128(reg2, _mm_srli_si128(reg2, 8));
			__m128i multiplicationLow = _mm_clmulepi64_si128(reg1, reg2, 0x00);
			__m128i multiplicationHigh = _mm_clmulepi64_si128(reg1, reg2, 0x11);
			__m128i multiplicationMiddle = _mm_clmulepi64_si128(cache1, cache2, 0x00);
			multiplicationMiddle = _mm_xor_si128(_mm_xor_si128(multiplicationMiddle, multiplicationLow), multiplicationHigh);

			multiplicationLow = _mm_xor_si128(multiplicationLow, _mm_slli_si128(multiplicationMiddle, 8));
			multiplicationHigh = _mm_xor_si128(multiplicationHigh, _mm_srli_si128(multiplicationMiddle, 8));
			resultLow = _mm_xor_si128(resultLow, multiplicationLow);
			resultHigh = _mm_xor_si128(resultHigh, multiplicationHigh);

			_mm_storeu_si128((__m128i*)(&result[i + j]), resultLow);
			_mm_storeu_si128((__m128i*)(&result[i + j + 2]), resultHigh);
		}
	}
	if ((n & 1) != 0)
	{	//now: n = m + 1 
		result[m * 2] ^= carrylessMultipyAndXor(input2, input1[m], &result[m], m);
		result[m * 2 + 1] ^= carrylessMultipyAndXor(input1, input2[m], &result[m], n);
	}
}

static void	carrylessAdd(uint64* input1, uint64* input2, uint64* result, int n)
{
	if ((n & 1) != 0)
	{
		result[n - 1] = input1[n - 1] ^ input2[n - 1];
	}
	for (int i = n & -2; (i -= 2) >= 0; )
	{
		_mm_storeu_si128((__m128i*)(&result[i]), _mm_xor_si128(_mm_loadu_si128((__m128i*)(&input1[i])), _mm_loadu_si128((__m128i*)(&input2[i]))));
	}
}

static void carrylessKaratsuba(uint64* input1, uint64* input2, uint64* result, int n, uint64* temporaryBuffer)
{	//(ax+b)(cx+d) = ac xx + bd + x[ac+bd-(b-a)(d-c)]
	if (n < CarrylessKaratsubaBreakPoint)
	{
		carrylessMultiplication(input1, input2, result, n);
		return;
	}

	int hi_n = n >> 1;
	int lo_n = (n + 1) >> 1;
	carrylessKaratsuba(&input1[0], &input2[0], &result[0], lo_n, temporaryBuffer);				// b * d
	carrylessKaratsuba(&input1[lo_n], &input2[lo_n], &result[lo_n * 2], hi_n, temporaryBuffer);	// a * c

	carrylessAdd(&input1[0], &input1[lo_n], &temporaryBuffer[0], hi_n);		// b - a
	carrylessAdd(&input2[0], &input2[lo_n], &temporaryBuffer[lo_n], hi_n);	// d - c
	if (lo_n != hi_n)
	{
		temporaryBuffer[hi_n] = input1[hi_n];
		temporaryBuffer[lo_n + hi_n] = input2[hi_n];
	}
	carrylessKaratsuba(&temporaryBuffer[0], &temporaryBuffer[lo_n], &temporaryBuffer[lo_n * 2], lo_n, &temporaryBuffer[lo_n * 4]);

	carrylessAdd(&result[0], &temporaryBuffer[lo_n * 2], &temporaryBuffer[0], lo_n * 2);
	carrylessAdd(&temporaryBuffer[0], &result[lo_n * 2], &temporaryBuffer[0], hi_n * 2);
	carrylessAdd(&result[lo_n], &temporaryBuffer[0], &result[lo_n], lo_n * 2);
}

static void fastestCarrylessMultiplication(uint64* input1, uint64* input2, uint64* result, int n)
{
	if (n < CarrylessFourierBreakPoint)
	{
		int tempSize = GetKaratsubaCarrylessMultiplicationBufferSize(n);
		uint64* temporaryBuffer = new uint64[tempSize];
		carrylessKaratsuba(input1, input2, result, n, temporaryBuffer);
		delete[] temporaryBuffer;
	}
	else
	{
		int tempSize = GetKaratsubaCarrylessMultiplicationBufferSize(n);
		uint64* temporaryBuffer = new uint64[tempSize];
		carrylessKaratsuba(input1, input2, result, n, temporaryBuffer);
		delete[] temporaryBuffer;
	}
}

extern "C"
{
	X64_API void __cdecl FourierSquare(uint64* input, uint64* result, int n)
	{
		fourierSquare(input, result, n);
	}

	X64_API void __cdecl FourierMultiplication(uint64* input1, uint64* input2, uint64* result, int n)
	{
		fourierMultiplication(input1, input2, result, n);
	}

	X64_API void __cdecl FastestMultiplication(uint64* input1, uint64* input2, uint64* result, int n, bool isSigned = false)
	{
		fastestMultiplication(input1, input2, result, n);
		if (!isSigned)
		{
			return;
		}
		if ((input1[n - 1] & Bit63) != 0)
		{
			Subtract(&result[n], input2, &result[n], 0, n);
		}
		if ((input2[n - 1] & Bit63) != 0)
		{
			Subtract(&result[n], input1, &result[n], 0, n);
		}
	}

	X64_API uint64 __cdecl CarrylessMultipyAndXor(uint64* input, uint64 digit, uint64* result, int n)
	{
		return carrylessMultipyAndXor(input, digit, result, n);
	}

	X64_API void __cdecl CarrylessMultiplication(uint64* input1, uint64* input2, uint64* result, int n)
	{
		fastestCarrylessMultiplication(input1, input2, result, n);
	}

	X64_API void __cdecl FastestSquare(uint64* input, uint64* result, int n, bool isSigned = false)
	{
		fastestSquare(input, result, n);
		if (isSigned && (input[n - 1] & Bit63) != 0)
		{
			SubShl(&result[n], input, 1, n);
		}
	}

	///returns a = a * b mod 2^128 - 2^54 + 1
	X64_API void __cdecl mulmod(uint64* a, uint64* b)
	{	//2^128 == 2^54 - 1 mod 'P' ; P = 2^128 - 2^54 + 1
		//r3 r2 r1 r0 = a * b
		//   t1 t0
		//   s1 s0
		uint64 r1, r0 = _umul128(a[0], b[0], &r1);
		uint64 r3, r2 = _umul128(a[1], b[1], &r3);
		uint64 t1, t0 = _umul128(a[0], b[1], &t1);
		uint64 s1, s0 = _umul128(a[1], b[0], &s1);

		uint8 carry;
		carry = _addcarry_u64(0, r1, t0, &r1);
		carry = _addcarry_u64(carry, r2, t1, &r2);
		_addcarry_u64(carry, r3, 0ULL, &r3);

		carry = _addcarry_u64(0, r1, s0, &r1);
		carry = _addcarry_u64(carry, r2, s1, &r2);
		_addcarry_u64(carry, r3, 0ULL, &r3);

		//X == (r3,r2) << 54
		uint64 x2 = r3 >> 10;	//10 == 64 - 54
		uint64 x1 = (r3 << 54) | (r2 >> 10);
		uint64 x0 = r2 << 54;

		uint8 borrow;
		borrow = _subborrow_u64(0, x0, r2, &x0);
		borrow = _subborrow_u64(borrow, x1, r3, &x1);
		_subborrow_u64(borrow, x2, 0ULL, &x2);	//x[0,1,2] = (W << 54) - W

		carry = _addcarry_u64(0, r0, x0, &r0);
		carry = _addcarry_u64(carry, r1, x1, &r1);
		_addcarry_u64(carry, x2, 0ULL, &x2);	//[r0,r1,x2] = x[0,1,2] + r[0,1]

		x0 = x2 << 54;
		x1 = x2 >> 10;	//(x2 << 128) == (x2 << 54) - x2
		borrow = _subborrow_u64(0, x0, x2, &x0);
		_subborrow_u64(borrow, x1, 0ULL, &x1);	//x[0,1] = (x2 << 54) - x2

		carry = _addcarry_u64(0, r0, x0, &r0);
		carry = _addcarry_u64(carry, r1, x1, &r1);
		while (carry != 0)
		{
			carry = _addcarry_u64(0, r0, (1ULL << 54) - 1ULL, &r0);
			carry = _addcarry_u64(carry, r1, 0ULL, &r1);
		}

		a[0] = r0;
		a[1] = r1;
	}

	///returns a = a + b and b = a - b in parallel.
	X64_API void __cdecl butterfly(uint64* a, uint64* b)
	{
		uint64 r0 = a[0], r1 = a[1];
		uint64 s0 = b[0], s1 = b[1];

		uint64 y0, y1;
		uint8 borrow = _subborrow_u64(0, r0, s0, &y0);
		borrow = _subborrow_u64(borrow, r1, s1, &y1);
		while (borrow != 0)
		{
			borrow = _subborrow_u64(0, y0, (1ULL << 54) - 1ULL, &y0);
			borrow = _subborrow_u64(borrow, y1, 0ULL, &y1);
		}
		b[0] = y0;
		b[1] = y1;

		uint8 carry = _addcarry_u64(0, r0, s0, &r0);
		carry = _addcarry_u64(carry, r1, s1, &r1);
		while (carry != 0)
		{
			carry = _addcarry_u64(0, r0, (1ULL << 54) - 1ULL, &r0);
			carry = _addcarry_u64(carry, r1, 0ULL, &r1);
		}
		a[0] = r0;
		a[1] = r1;
	}

	///returns a = a + b modulo P = 2^128 - 2^54 + 1.
	X64_API void __cdecl addmod(uint64* a, uint64* b)
	{
		uint64 r0 = a[0], r1 = a[1];

		uint8 carry = _addcarry_u64(0, r0, b[0], &r0);
		carry = _addcarry_u64(carry, r1, b[1], &r1);
		while (carry != 0)
		{
			carry = _addcarry_u64(0, r0, (1ULL << 54) - 1ULL, &r0);
			carry = _addcarry_u64(carry, r1, 0ULL, &r1);
		}
		a[0] = r0;
		a[1] = r1;
	}

	///returns a = a - b and b = a - b in parallel.
	X64_API void __cdecl submod(uint64* a, uint64* b)
	{
		uint64 r0 = a[0], r1 = a[1];

		uint8 borrow = _subborrow_u64(0, r0, b[0], &r0);
		borrow = _subborrow_u64(borrow, r1, b[1], &r1);
		while (borrow != 0)
		{
			borrow = _subborrow_u64(0, r0, (1ULL << 54) - 1ULL, &r0);
			borrow = _subborrow_u64(borrow, r1, 0ULL, &r1);
		}
		a[0] = r0;
		a[1] = r1;
	}

	X64_API void __cdecl PowMod(uint64 *a, uint64* power)
	{
		if (power[1] == 0 && power[0] <= 1)
		{
			if (power[0] == 0)
			{
				a[0] = 1ULL;
				a[1] = 0ULL;
			}
			return;
		}
		uint64 savea[2] = { a[0], a[1] };
		uint64 halfPower[2] = { (power[0] >> 1) | (power[1] << 63), power[1] >> 1 };
		PowMod(a, halfPower);
		mulmod(a, a);
		if ((power[0] & 1ULL) != 0)
		{
			mulmod(a, savea);
		}
	}

	X64_API void __cdecl forwardFFTCore(uint64* number, int startIndex, int coreLength, uint64* wpower)
	{
		startIndex <<= 1;
		coreLength <<= 1;
		int middleIndex = startIndex + coreLength;
		for (int i = coreLength; (i -= 2) >= 0; )
		{
			mulmod(&number[middleIndex + i], wpower);
			butterfly(&number[startIndex + i], &number[middleIndex + i]);
		}
	}

	X64_API void __cdecl inverseFFTCore(uint64* number, int startIndex, int coreLength, uint64* wpower)
	{
		startIndex <<= 1;
		coreLength <<= 1;
		int middleIndex = startIndex + coreLength;
		for (int i = coreLength; (i -= 2) >= 0; )
		{
			butterfly(&number[startIndex + i], &number[middleIndex + i]);
			mulmod(&number[middleIndex + i], wpower);
		}
	}

	X64_API void __cdecl butterflyCore(uint64* number, int startIndex, int coreLength)
	{
		startIndex <<= 1;
		coreLength <<= 1;
		int middleIndex = startIndex + coreLength;
		for (int i = coreLength; (i -= 2) >= 0; )
		{
			butterfly(&number[startIndex + i], &number[middleIndex + i]);
		}
	}

	X64_API void __cdecl PropagateCarries(uint64* number, int bitsPerDigit, int numberLength)
	{
		uint64 carry0 = 0, carry1 = 0;
		uint64 mask = (1ULL << bitsPerDigit) - 1ULL;
		numberLength <<= 1;
		int reverseShift = 64 - bitsPerDigit;
		for (int i = 0; i < numberLength; i += 2)
		{
			remodP(&number[i]);
			_addcarry_u64(
				_addcarry_u64(0, carry0, number[i], &carry0),
				carry1, number[i + 1], &carry1);
			number[i] = carry0 & mask;
			number[i + 1] = 0ULL;
			carry0 = (carry0 >> bitsPerDigit) | (carry1 << reverseShift);
			carry1 >>= bitsPerDigit;
		}
	}

	X64_API void __cdecl OpMov(__m256i* destination, __m256i* source, int log2BitsPerOp, int n)
	{
		memmove(destination, source, n << (log2BitsPerOp >> 3));
	}

	X64_API void __cdecl OpNot(__m256i* destination, int log2BitsPerOp, int n)
	{
		__int64 bytes = (__int64)n << (log2BitsPerOp >> 3);
		__int64* sourceDestination = (__int64*)destination;
		for (__int64 i = bytes >> 3; --i >= 0; )
		{
			sourceDestination[i] = ~sourceDestination[i];
		}
		__int8 *byteSourceDestination = (__int8*)&sourceDestination[bytes >> 3];
		for (int i = bytes & 7; --i >= 0; )
		{
			byteSourceDestination[i] = ~byteSourceDestination[i];
		}
	}

#define REPEAT(count, log2Groups, operation, member, directOperation)	\
	for (int i = (count) >> (log2Groups); --i >= 0; )	{ _mm256_storeu_si256(&destination[i], operation(_mm256_loadu_si256(&source1[i]), _mm256_loadu_si256(&source2[i]))); }	\
	for (int j = (count) >> (log2Groups), i = (count) & ((1 << (log2Groups)) - 1); --i >= 0; ) { destination[j].##member[i] = source1[j].##member[i] ##directOperation source2[j].##member[i]; }
#define REPEAT256bit(count, operation)	\
	for (int i = (count); --i >= 0; )	{ _mm256_storeu_si256(&destination[i], operation(_mm256_loadu_si256(&source1[i]), _mm256_loadu_si256(&source2[i]))); }
#define uint64PointerTo(memory, index)	&(((uint64*)memory)[index])

	X64_API void __cdecl OpAdd(__m256i* source1, __m256i* source2, __m256i* destination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 3:	REPEAT(n, 5, _mm256_add_epi8, m256i_u8, +); break;
		case 4:	REPEAT(n, 4, _mm256_add_epi16, m256i_u16, +); break;
		case 5:	REPEAT(n, 3, _mm256_add_epi32, m256i_u32, +); break;
		case 6:	REPEAT(n, 2, _mm256_add_epi64, m256i_u64, +); break;
		default:
			int shift = log2BitsPerOp - 6;
			int k = 1 << shift;
			for (int i = n; --i >= 0; )
			{
				int index = i << shift;
				Add(uint64PointerTo(source1, index), uint64PointerTo(source2, index), uint64PointerTo(destination, index), 0, k);
			}
			break;
		}
	}

	X64_API void __cdecl OpSub(__m256i* source1, __m256i* source2, __m256i* destination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 3:	REPEAT(n, 5, _mm256_sub_epi8, m256i_u8, -); break;
		case 4:	REPEAT(n, 4, _mm256_sub_epi16, m256i_u16, -); break;
		case 5:	REPEAT(n, 3, _mm256_sub_epi32, m256i_u32, -); break;
		case 6:	REPEAT(n, 2, _mm256_sub_epi64, m256i_u64, -); break;
		default:
			int shift = log2BitsPerOp - 6;
			int k = 1 << shift;
			for (int i = n; --i >= 0; )
			{
				int index = i << shift;
				Subtract(uint64PointerTo(source1, index), uint64PointerTo(source2, index), uint64PointerTo(destination, index), 0, k);
			}
			break;
		}
	}

	X64_API void __cdecl OpAnd(__m256i* source1, __m256i* source2, __m256i* destination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 3:	REPEAT(n, 5, _mm256_and_si256, m256i_u8, &); break;
		case 4:	REPEAT(n, 4, _mm256_and_si256, m256i_u16, &); break;
		case 5:	REPEAT(n, 3, _mm256_and_si256, m256i_u32, &); break;
		case 6:	REPEAT(n, 2, _mm256_and_si256, m256i_u64, &); break;
		default:
			REPEAT256bit(n << (log2BitsPerOp - 6), _mm256_and_si256); break;
		}
	}

	X64_API void __cdecl OpOr(__m256i* source1, __m256i* source2, __m256i* destination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 3:	REPEAT(n, 5, _mm256_or_si256, m256i_u8, | ); break;
		case 4:	REPEAT(n, 4, _mm256_or_si256, m256i_u16, | ); break;
		case 5:	REPEAT(n, 3, _mm256_or_si256, m256i_u32, | ); break;
		case 6:	REPEAT(n, 2, _mm256_or_si256, m256i_u64, | ); break;
		default:
			REPEAT256bit(n << (log2BitsPerOp - 6), _mm256_or_si256); break;
			break;
		}
	}

	X64_API void __cdecl OpXor(__m256i* source1, __m256i* source2, __m256i* destination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 3:	REPEAT(n, 5, _mm256_xor_si256, m256i_u8, ^); break;
		case 4:	REPEAT(n, 4, _mm256_xor_si256, m256i_u16, ^); break;
		case 5:	REPEAT(n, 3, _mm256_xor_si256, m256i_u32, ^); break;
		case 6:	REPEAT(n, 2, _mm256_xor_si256, m256i_u64, ^); break;
		default:
			REPEAT256bit(n << (log2BitsPerOp - 6), _mm256_xor_si256); break;
			break;
		}
	}

	X64_API void __cdecl OpMul(__m256i* source1, __m256i* source2, __m256i* destination, int log2BitsPerOp, int n)
	{
		if (log2BitsPerOp >= 6)
		{
			int shift = log2BitsPerOp - 6;
			int k = 1 << shift;
			for (int i = n; --i >= 0; )
			{
				int index = i << shift;
				fastestMultiplication(uint64PointerTo(source1, index), uint64PointerTo(source2, index), uint64PointerTo(destination, index << 1), k);
			}
			return;
		}
	}

	X64_API void __cdecl OpCLMul(__m256i* source1, __m256i* source2, __m256i* destination, int log2BitsPerOp, int n)
	{
		if (log2BitsPerOp >= 6)
		{
			int shift = log2BitsPerOp - 6;
			int k = 1 << shift;
			for (int i = n; --i >= 0; )
			{
				int index = i << shift;
				carrylessMultiplication(uint64PointerTo(source1, index), uint64PointerTo(source2, index), uint64PointerTo(destination, index << 1), k);
			}
			return;
		}
	}

#define REPEAT_FP(count, itemsInRegister, loadop, storeop, operation, directOperation)	\
	for (int i = (count) & (-(itemsInRegister)); (i -= (itemsInRegister)) >= 0; )	{ storeop(&destination[i], operation(loadop(&source1[i]), loadop(&source2[i]))); }	\
	for (int j = (count) - 1, i = (count) & ((itemsInRegister) - 1); --i >= 0; j--) { destination[j] = source1[j] ##directOperation source2[j]; }

	X64_API void __cdecl OpFPAdd(__m256* inputSource1, __m256* inputSource2, __m256* inputDestination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 5:
		{
			float* source1 = (float*)inputSource1;
			float* source2 = (float*)inputSource2;
			float* destination = (float*)inputDestination;
			REPEAT_FP(n, 8, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_add_ps, +);
			break;
		}
		case 6:
		{
			double* source1 = (double*)inputSource1;
			double* source2 = (double*)inputSource2;
			double* destination = (double*)inputDestination;
			REPEAT_FP(n, 4, _mm256_loadu_pd, _mm256_storeu_pd, _mm256_add_pd, +);
			break;
		}
		default:
			break;
		}
	}

	X64_API void __cdecl OpFPSub(__m256* inputSource1, __m256* inputSource2, __m256* inputDestination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 5:
		{
			float* source1 = (float*)inputSource1;
			float* source2 = (float*)inputSource2;
			float* destination = (float*)inputDestination;
			REPEAT_FP(n, 8, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_sub_ps, -);
			break;
		}
		case 6:
		{
			double* source1 = (double*)inputSource1;
			double* source2 = (double*)inputSource2;
			double* destination = (double*)inputDestination;
			REPEAT_FP(n, 4, _mm256_loadu_pd, _mm256_storeu_pd, _mm256_sub_pd, -);
			break;
		}
		default:
			break;
		}
	}

	X64_API void __cdecl OpFPMul(__m256* inputSource1, __m256* inputSource2, __m256* inputDestination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 5:
		{
			float* source1 = (float*)inputSource1;
			float* source2 = (float*)inputSource2;
			float* destination = (float*)inputDestination;
			REPEAT_FP(n, 8, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_mul_ps, *);
			break;
		}
		case 6:
		{
			double* source1 = (double*)inputSource1;
			double* source2 = (double*)inputSource2;
			double* destination = (double*)inputDestination;
			REPEAT_FP(n, 4, _mm256_loadu_pd, _mm256_storeu_pd, _mm256_mul_pd, *);
			break;
		}
		default:
			break;
		}
	}

	X64_API void __cdecl OpFPDiv(__m256* inputSource1, __m256* inputSource2, __m256* inputDestination, int log2BitsPerOp, int n)
	{
		switch (log2BitsPerOp)
		{
		case 5:
		{
			float* source1 = (float*)inputSource1;
			float* source2 = (float*)inputSource2;
			float* destination = (float*)inputDestination;
			REPEAT_FP(n, 8, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_div_ps, / );
			break;
		}
		case 6:
		{
			double* source1 = (double*)inputSource1;
			double* source2 = (double*)inputSource2;
			double* destination = (double*)inputDestination;
			REPEAT_FP(n, 4, _mm256_loadu_pd, _mm256_storeu_pd, _mm256_div_pd, / );
			break;
		}
		default:
			break;
		}
	}

#undef REPEAT
#undef REPEAT_FP
#undef REPEAT256bit
#undef uint64PointerTo
};
