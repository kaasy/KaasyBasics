// AsmX64BitBasicOperations.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"

static uint8	inversion_table[128];
static bool 	inversionInitialized = false;

static uint64	invert(uint64 digit)
{	//return invert of X mod 2^64, note that X must be odd...
	if (!inversionInitialized)
	{
		for (int i = 128; --i >= 0; )
		{
			uint8 byteDigit = (uint8)(i * 2 + 1);
			uint8 inverse = byteDigit & 0x0F;
			if ((byteDigit & 7) == 3 || (byteDigit & 7) == 5)
			{
				inverse ^= 8;
			}
			inverse *= 2 - byteDigit * inverse;
			inversion_table[i] = inverse;
		}
		inversionInitialized = true;
	}
	uint64	x = inversion_table[(digit >> 1) & 0x7F];
	x *= 2 - digit * x;
	x *= 2 - digit * x;
	x *= 2 - digit * x;
	return	x;
}

static int maximum(int x, int y)
{
	return x >= y ? x : y;
}

static int minimum(int x, int y)
{
	return x <= y ? x : y;
}

static sint64 absoluteValue(sint64 x)
{
	return x >= 0 ? x : -x;
}

static int bsr(uint64 digit)
{
	unsigned long bit_position;
	_BitScanReverse64(&bit_position, digit);
	return bit_position;
}

static int bsf(uint64 digit)
{
	unsigned long bit_position;
	_BitScanForward64(&bit_position, digit);
	return bit_position;
}

//return position of the highest set bit of A or -1 if A[]= 0
static int high_bit(uint64* a, int n)
{
	for (int i = n; --i >= 0;)
	{
		if (a[i] != 0)
		{
			return  i * 64 + bsr(a[i]);
		}
	}
	return  -1;
}

//return position of the lowest set bit of A or n * 64 if a == 0.
static int low_bit(uint64* a, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (a[i] != 0)
		{
			return  i * 64 + bsf(a[i]);
		}
	}
	return n << 6;
}

int     int_digits(uint64* a, int n)
{	//return number of digits the number has...
	for (int i = n; --i >= 0;)
	{
		if (a[i] != 0)
		{
			return i + 1;
		}
	}
	return 0;
}

//use shift in [0..63]
static	uint64	add_shl(uint64* inputOutput, uint64 *input2, uint8 shift, int n)
{
	uint64 carry = 0;
	uint8 bitCarry = 0;
	for (int i = 0; i < n; i++)
	{
		uint64 nextValue = input2[i];
		carry = __shiftleft128(carry, nextValue, shift);
		bitCarry = _addcarry_u64(bitCarry, inputOutput[i], carry, &inputOutput[i]);
		carry = nextValue;
	}
	carry = __shiftleft128(carry, 0, shift) + bitCarry;
	return carry;
}

//use shift in [0..63]
static	uint64	sub_shl(uint64* inputOutput, uint64 *input2, uint8 shift, int n)
{
	uint64 borrow = 0;
	uint8 bitBorrow = 0;
	for (int i = 0; i < n; i++)
	{
		uint64 nextValue = input2[i];
		borrow = __shiftleft128(borrow, nextValue, shift);
		bitBorrow = _subborrow_u64(bitBorrow, inputOutput[i], borrow, &inputOutput[i]);
		borrow = nextValue;
	}
	borrow = __shiftleft128(borrow, 0, shift) + bitBorrow;
	return borrow;
}

static void makeZero(uint64* destination, int n)
{
	for (int i = n; --i >= 0;)
	{
		destination[i] = 0;
	}
}

static void let(uint64* destination, uint64 digit, int n)
{
	for (int i = n; --i >= 0;)
	{
		destination[i] = digit;
	}
}

//returns overflow part only for the higher digit if shift is fractional (lower 6 bits != 0)
static	uint64	int_shl(uint64* input, int _shift, int n)
{
	int blocks = _shift >> 6;
	if (blocks != 0)
	{
		for (int i = n - blocks; --i >= 0; )
		{
			input[i + blocks] = input[i];
		}
		makeZero(input, blocks);
	}
	uint8 shift = (uint8)(_shift & 63);
	uint64 carry = 0;
	if (shift != 0)
	{
		for (int i = blocks; i < n; i++)
		{
			uint64 nextValue = input[i];
			input[i] = __shiftleft128(carry, nextValue, shift);
			carry = nextValue;
		}
		carry = __shiftleft128(carry, 0, shift);
	}
	return carry;
}

//use shift >= 0, returns overflow part only for 1 digit/register value.
static	uint64	int_shr(uint64* input, int _shift, int n, bool isSigned = false)
{
	uint64 borrow = isSigned && n > 0 && (input[n - 1] & Bit63) != 0 ? MaxValueUint64 : 0ULL;
	int blocks = _shift >> 6;
	int end = n - blocks;
	if (blocks != 0)
	{
		for (int i = 0; i < end; i++)
		{
			input[i] = input[i + blocks];
		}
		let(&input[end], borrow, blocks);
	}
	uint8 shift = (uint8)(_shift & 63);
	if (shift != 0)
	{
		for (int i = end; --i >= 0;)
		{
			uint64 previousValue = input[i];
			input[i] = __shiftright128(previousValue, borrow, shift);
			borrow = previousValue;
		}
		borrow = __shiftright128(0, borrow, shift);
	}
	return borrow;
}

static	bool	isZero(uint64* number, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (number[i] != 0)
		{
			return false;
		}
	}
	return true;
}

static	void	let(uint64* destination, uint64 *source, int n)
{
	if (destination > source)
	{
		for (int i = n; --i >= 0;)
		{
			destination[i] = source[i];
		}
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			destination[i] = source[i];
		}
	}
}

unsigned char udiv128Data[] =
{
	//extern "C" digit64 udiv128(digit64 low, digit64 hi, digit64 divisor);
	//; Arguments
	//; RCX       Low Digit
	//; RDX       High Digit
	//; R8        Divisor
	//; RAX       Quotient upon return
	0x48, 0x89, 0xC8, // mov rax,rcx
	0x49, 0xF7, 0xF0, // div r8
	0xC3              // ret
};

uint64(__fastcall *udiv128)(uint64 numlo, uint64 numhi, uint64 divisor) = (uint64(__fastcall *)(uint64, uint64, uint64))(void*)udiv128Data;
static bool isFunctionProtected = false;

uint64	digit_divide(uint64 low, uint64 high, uint64 div)
{       //return [high:low] / div with overflow correction
	if (!isFunctionProtected)
	{
		DWORD dummy;
		VirtualProtect(udiv128Data, sizeof(udiv128Data), PAGE_EXECUTE_READWRITE, &dummy);
		isFunctionProtected = true;
	}

	if (high == div)
	{
		return (uint64)-1;
	}
	uint64 result = udiv128(low, high, div);
	return result;
}

static uint64	getSquareSum(uint64* ai, uint64* aj, uint64* ck, uint64 carry)
{	//operation: C[K]= sum(A[I--] * A[J++]), return carry
	uint64 sum2 = 0, sum1 = 0, sum0 = 0;
	uint8 bitCarry = 0;
	while (ai > aj)
	{
		uint64 hi, lo = _umul128(ai[0], aj[0], &hi);
		_addcarry_u64(_addcarry_u64(_addcarry_u64(0, sum0, lo, &sum0), sum1, hi, &sum1), sum2, 0, &sum2);
		ai--;
		aj++;
	}
	_addcarry_u64(_addcarry_u64(_addcarry_u64(0, sum0, sum0, &sum0), sum1, sum1, &sum1), sum2, sum2, &sum2);
	if (ai == aj)
	{
		uint64 hi, lo = _umul128(ai[0], aj[0], &hi);
		_addcarry_u64(_addcarry_u64(_addcarry_u64(0, sum0, lo, &sum0), sum1, hi, &sum1), sum2, 0, &sum2);
	}
	_addcarry_u64(_addcarry_u64(_addcarry_u64(0, ck[0], sum0, &ck[0]), carry, sum1, &ck[1]), sum2, 0, &carry);
	return carry;
}

static void negate(uint64* result, int n)
{
	uint8 carry = 1;
	while (--n >= 0)
	{
		carry = _addcarry_u64(carry, ~result[0], 0, result);
		result++;
	}
}

static int chooseStartPosition(uint64* power, int bitPosition, int cacheBits)
{	//return j <= bitPosition so bit(power, j) is set.
	int	j = bitPosition - cacheBits;
	if (j < 0)
	{
		j = 0;
	}
	while ((power[j >> 6] & (1ULL << (j & 63))) == 0)
	{
		j++;
	}
	return j;
}

static bool getBitAt(uint64* number, int bitIndex)
{
	return	(number[bitIndex >> 6] & (1ULL << (bitIndex & 63))) != 0;
}

static uint64 getBitsValue(uint64* number, int startPosition, int endPosition)
{
	startPosition = maximum(startPosition, 0);
	if (endPosition < startPosition)
	{
		return 0;
	}
	int startIndex = startPosition >> 6;
	uint64 result = number[startIndex] >> (startPosition & 63);
	int endIndex = endPosition >> 6;
	if (endIndex != startIndex)
	{
		result |= number[endIndex] << ((-startPosition) & 63);
	}
	int length = endPosition - startPosition + 1;
	if (length < 64)
	{
		result &= (1ULL << length) - 1;
	}
	return result;
}

extern "C"
{
	X64_API uint8 __cdecl AddDigit(uint64* result, uint64 digit, int n);
	X64_API uint8 __cdecl SubtractDigit(uint64* result, uint64 digit, int n);
	X64_API uint8 __cdecl Add(uint64* input1, uint64* input2, uint64* result, uint8 carry, int n);
	X64_API uint8 __cdecl Subtract(uint64* input1, uint64* input2, uint64* result, uint8 borrow, int n);
	X64_API void __cdecl Multiply(uint64* input1, uint64* input2, uint64* result, int n);
	X64_API void __cdecl Square(uint64* input, uint64* result, int n);
	X64_API uint64 __cdecl MultiplyDigitAndAdd(uint64* input, uint64 digit, uint64* result, int n);
	X64_API uint64 __cdecl MultiplyDigitAndSubtract(uint64* input, uint64 digit, uint64* result, int n);
}

static void shrModP(uint64* number, int m, int inputShift, uint64* mod, int n, uint64 inverseP0)
{
	bool isNegative = m > n && ((sint64)number[m - 1]) < 0;
	if (isNegative)
	{
		negate(number, m);
	}
	for (int shift = inputShift; shift > 0; shift -= 64)
	{	//number -= constant * P so that number[0] == 0 after the subtraction.
		uint64 mask = shift >= 64 ? 0ULL - 1ULL : (1ULL << shift) - 1ULL;
		uint64 borrow = MultiplyDigitAndSubtract(mod, number[0] * inverseP0 & mask, number, n);
		int_shr(number, shift >= 64 ? 64 : shift, m);
		if (SubtractDigit(&number[n - 1], shift >= 64 ? borrow : borrow << (64 - shift), 1) != 0)
		{
			while (Add(number, mod, number, 0, n) == 0);
		}
	}
	if (isNegative)
	{
		if (Subtract(mod, number, number, 0, n) != 0)
		{
			while (Add(number, mod, number, 0, n) == 0);
		}
	}
}

static void getSummation(uint64* x, uint64* _mx, uint64* result, int n)
{
	if (((sint64)_mx[1]) < 0)
	{
		uint64 mx[2] = { _mx[0], _mx[1] };
		negate(mx, 2);
		SubtractDigit(&result[n], MultiplyDigitAndSubtract(x, mx[0], &result[0], n), 2);
		if (mx[1] != 0)
		{
			SubtractDigit(&result[n + 1], MultiplyDigitAndSubtract(x, mx[1], &result[1], n), 1);
		}
	}
	else
	{
		AddDigit(&result[n], MultiplyDigitAndAdd(x, _mx[0], &result[0], n), 2);
		if (_mx[1] != 0)
		{
			AddDigit(&result[n + 1], MultiplyDigitAndAdd(x, _mx[1], &result[1], n), 1);
		}
	}
}

//operation input %= mod, fill quotient if quotient != NULL; (input[m], mod[n]) should be unsigned.
static void int_mod(uint64* input, int m, uint64* mod, int n, uint64* quotient = NULL)
{
	int i = int_digits(input, m),
		j = int_digits(mod, n);
	//zero(DIV, N* 2);
	if (i <= 0 || j <= 0 || j > i)
	{
		return;         //DIVISION BY 0 EXCEPTION OR input < mod
	}
	int shift = 63 - bsr(mod[j - 1]);
	int_shl(mod, shift, j);   //shifted B improves performance
	uint64
		hi_A = int_shl(input, shift, i),
		*fix = &hi_A;    //use one digit to fix div loop below
	if (hi_A == 0 && input[i - 1] < mod[j - 1])
	{
		fix = &input[--i];   //DIV_RESULT= 0 in this case
	}
	if (quotient != NULL)
	{
		makeZero(&quotient[i - j + 1], m - (i - j + 1));
	}
	for (int k = i - j + 1; --k >= 0;)
	{
		uint64 digit = digit_divide(input[k + j - 1], *fix, mod[j - 1]);
		(*fix) -= MultiplyDigitAndSubtract(mod, digit, &input[k], j);
		while (*fix != 0)
		{
			(*fix) += Add(&input[k], mod, &input[k], 0, j);
			digit--;
			//theoretically: loop maximum 2 times because mod[j - 1] has highest bit set.
		}
		if (quotient != NULL)
		{
			quotient[k] = digit;
		}
		fix = &input[k + j - 1];		//now fix digits of A (mak'em 0)
	}
	int_shr(mod, shift, j);
	int_shr(input, shift, j);			//A just got smaller
}

static void karatsuba(uint64* input1, uint64* input2, uint64* result, int n, uint64* temporaryBuffer)
{	//(ax+b)(cx+d) = ac xx + bd + x[ac+bd-(b-a)(d-c)]
	if (n < KaratsubaBreakPoint)
	{
		Multiply(input1, input2, result, n);
		return;
	}

	int hi_n = n >> 1;
	int lo_n = (n + 1) >> 1;
	karatsuba(&input1[0], &input2[0], &result[0], lo_n, temporaryBuffer);				// b * d
	karatsuba(&input1[lo_n], &input2[lo_n], &result[lo_n * 2], hi_n, temporaryBuffer);	// a * c

	uint8 carryA = Subtract(&input1[0], &input1[lo_n], &temporaryBuffer[0], 0, hi_n);		// b - a
	uint8 carryB = Subtract(&input2[0], &input2[lo_n], &temporaryBuffer[lo_n], 0, hi_n);	// d - c
	if (lo_n != hi_n)
	{
		carryA = Subtract(&input1[hi_n], zero, &temporaryBuffer[hi_n], carryA, 1);
		carryB = Subtract(&input2[hi_n], zero, &temporaryBuffer[lo_n + hi_n], carryB, 1);
	}
	bool minus = true;
	if (carryA != 0)
	{
		negate(&temporaryBuffer[0], lo_n);
		minus = !minus;
	}
	if (carryB != 0)
	{
		negate(&temporaryBuffer[lo_n], lo_n);
		minus = !minus;
	}
	karatsuba(&temporaryBuffer[0], &temporaryBuffer[lo_n], &temporaryBuffer[lo_n * 2], lo_n, &temporaryBuffer[lo_n * 4]);

	uint8 borrow1 = 0, carry1 = 0;
	if (minus)
	{
		borrow1 = Subtract(&result[0], &temporaryBuffer[lo_n * 2], &temporaryBuffer[0], 0, lo_n * 2);
	}
	else
	{
		carry1 = Add(&result[0], &temporaryBuffer[lo_n * 2], &temporaryBuffer[0], 0, lo_n * 2);
	}	//temp[0..] = bd - (b-a)(d-c)
	uint8 carry2 = Add(&temporaryBuffer[0], &result[lo_n * 2], &temporaryBuffer[0], 0, hi_n * 2);
	if (lo_n != hi_n)
	{
		carry2 = Add(&temporaryBuffer[hi_n * 2], zero, &temporaryBuffer[hi_n * 2], carry2, 2);
	}	//temp[0..] = ac + bd - (b-a)(d-c) = ad + bc

	uint8 carry3 = Add(&result[lo_n], &temporaryBuffer[0], &result[lo_n], 0, lo_n * 2);
	uint8 generalCarry = carry1 + carry2 + carry3 - borrow1;
	if (generalCarry != 0)
	{
		AddDigit(&result[lo_n * 3], (uint64)generalCarry, n * 2 - lo_n * 3);
	}
}

static void karatsubaSquare(uint64* input, uint64* result, int n, uint64* temporaryBuffer)
{	//(ax+b)(ax+b) = aa xx + bb + x[aa+bb-(b-a)(b-a)]
	if (n < KaratsubaBreakPoint)
	{
		Square(input, result, n);
		return;
	}

	int hi_n = n >> 1;
	int lo_n = (n + 1) >> 1;
	karatsubaSquare(&input[0], &result[0], lo_n, temporaryBuffer);				// b * b
	karatsubaSquare(&input[lo_n], &result[lo_n * 2], hi_n, temporaryBuffer);	// a * a

	uint8 borrow0 = Subtract(&input[0], &input[lo_n], &temporaryBuffer[0], 0, hi_n);		// b - a
	if (lo_n != hi_n)
	{
		borrow0 = Subtract(&input[hi_n], zero, &temporaryBuffer[hi_n], borrow0, 1);
	}
	if (borrow0 != 0)
	{
		negate(&temporaryBuffer[0], lo_n);
	}
	karatsubaSquare(&temporaryBuffer[0], &temporaryBuffer[lo_n], lo_n, &temporaryBuffer[lo_n * 3]);

	uint8 borrow = Subtract(&result[0], &temporaryBuffer[lo_n], &temporaryBuffer[0], 0, lo_n * 2);
	uint8 carry2 = Add(&temporaryBuffer[0], &result[lo_n * 2], &temporaryBuffer[0], 0, hi_n * 2);
	if (lo_n != hi_n)
	{
		carry2 = Add(&temporaryBuffer[hi_n * 2], zero, &temporaryBuffer[hi_n * 2], carry2, 2);
	}	//temp[0..] = aa + bb - (b-a)(b-a)

	uint8 carry3 = Add(&result[lo_n], &temporaryBuffer[0], &result[lo_n], 0, lo_n * 2);
	uint8 generalCarry = carry2 + carry3 - borrow;
	if (generalCarry != 0)
	{
		AddDigit(&result[lo_n * 3], (uint64)generalCarry, n * 2 - lo_n * 3);
	}
}

static void	montgomery_reduction(uint64* a, uint64* mod, uint64 inverse_of_mod_0, int n)
{   //hi_A= A* (BASE^^-N) mod C
	for (int i = 0; i < n; i++)
	{
		a[i] = MultiplyDigitAndSubtract(mod, a[i] * inverse_of_mod_0, &a[i], n);
	}
	uint8    overflow = Subtract(a + n, a, a, 0, n);
	if (overflow != 0)
	{
		while (Add(a, mod, a, 0, n) == 0);
	}
}

#define	double_shift_right(lo, hi)	lo = ((lo) >> 1) | ((hi) << 63); hi = (uint64)(((sint64)hi) >> 1)
#define	subtract_with_borrow(lo1, hi1, lo2, hi2)	_subborrow_u64(_subborrow_u64(0, lo1, lo2, &lo1), hi1, hi2, &hi1);
#define	swap_numbers(first, second, temp)	uint64 temp = first; first = second; second = temp

static void estimateBinaryShift(
	uint64 hix, uint64 hiy, uint64 lox, uint64 loy,
	uint64 *resulta, uint64 *resultb, uint64 *resultc, uint64 *resultd)
{
	int count = 64;
	//x <- (x * a + y * b) >> 64
	//y <- (x * c + y * d) >> 64
	uint64 a1 = 1, b1 = 0, c1 = 0, d1 = 1;
	uint64 a0 = 0, b0 = 0, c0 = 0, d0 = 0;
	goto start;
	while (count > 0)
	{
		if (hix < hiy)
		{
			swap_numbers(hix, hiy, auxiliary1);
			swap_numbers(lox, loy, auxiliary2);
			swap_numbers(a0, c0, auxiliary3);
			swap_numbers(b0, d0, auxiliary4);
			swap_numbers(a1, c1, auxiliary5);
			swap_numbers(b1, d1, auxiliary6);
		}
		lox -= loy;
		hix -= hiy;
		subtract_with_borrow(a0, a1, c0, c1);
		subtract_with_borrow(b0, b1, d0, d1);
	start:
		if (lox == 0)
		{
			break;
		}
		while ((lox & 1) == 0 && count > 0)
		{
			count--;
			lox >>= 1;
			hix >>= 1;
			double_shift_right(a0, a1);
			double_shift_right(b0, b1);
		}
	}

	resulta[0] = a0; resulta[1] = a1;
	resultb[0] = b0; resultb[1] = b1;
	resultc[0] = c0; resultc[1] = c1;
	resultd[0] = d0; resultd[1] = d1;
}

#undef	swap_numbers
#undef	double_shift_right
#undef	subtract_with_borrow

unsigned char udivmod128Data[] =
{
	//extern "C" digit64 udivmod128(uint64 *lowAndQuotient, uint64 hi, uint64 divisor);
	//; Arguments
	//; RCX       *lowAndQuotient
	//; RDX       High Digit
	//; R8        Divisor
	//; RAX       Reminder upon return
	0x48, 0x8B, 0x01, // mov rax, [rcx]
	0x49, 0xF7, 0xF0, // div r8
	0x48, 0x89, 0x01, // mov [rcx], rax
	0x48, 0x89, 0xD0, // mov rax, rdx
	0xC3              // ret
};

uint64(__fastcall *udivmod128)(uint64 *lowAndQuotient, uint64 numhi, uint64 divisor) = (uint64(__fastcall *)(uint64*, uint64, uint64))(void*)udivmod128Data;
static bool isUdivmod128FunctionProtected = false;

extern "C"
{
	X64_API void __cdecl SetKaratsubaBreakPoint(int karatsubaBreakPoint)
	{
		KaratsubaBreakPoint = karatsubaBreakPoint;
	}

	X64_API void __cdecl SetFourierBreakPoint(int fourierBreakPoint)
	{
		FourierBreakPoint = fourierBreakPoint;
	}

	X64_API	int GetKaratsubaMultiplicationBufferSize(int n)
	{
		int tempSize = 0;
		for (int m = n; m >= KaratsubaBreakPoint;)
		{
			m = (m + 1) >> 1;
			tempSize += m * 4;
		}
		return tempSize;
	}

	X64_API	int GetKaratsubaCarrylessMultiplicationBufferSize(int n)
	{
		int tempSize = 0;
		for (int m = n; m >= CarrylessKaratsubaBreakPoint;)
		{
			m = (m + 1) >> 1;
			tempSize += m * 4;
		}
		return tempSize;
	}

	X64_API	int GetKaratsubaSquareBufferSize(int n)
	{
		int tempSize = 0;
		for (int m = n; m >= KaratsubaBreakPoint;)
		{
			m = (m + 1) >> 1;
			tempSize += m * 3;
		}
		return tempSize;
	}

	X64_API void __cdecl Negate(uint64* result, int n)
	{
		negate(result, n);
	}

	X64_API uint64 __cdecl IntShr(uint64* result, int shift, int n, bool isSigned = false)
	{
		return int_shr(result, shift, n, isSigned);
	}

	X64_API uint64 __cdecl IntShl(uint64* result, int shift, int n)
	{
		return int_shl(result, shift, n);
	}

	X64_API uint64 __cdecl AddShl(uint64* inputOutput, uint64 *input2, uint8 shift, int n)
	{
		return add_shl(inputOutput, input2, shift, n);
	}

	X64_API uint64 __cdecl SubShl(uint64* inputOutput, uint64 *input2, uint8 shift, int n)
	{
		return sub_shl(inputOutput, input2, shift, n);
	}

	X64_API uint8 __cdecl AddDigit(uint64* result, uint64 digit, int n)
	{
		uint8 carry = _addcarry_u64(0, result[0], digit, result);
		n--;
		result++;
		while (--n >= 0 && carry != 0)
		{
			carry = _addcarry_u64(carry, result[0], 0, result);
			result++;
		}
		return carry;
	}

	X64_API void __cdecl MultiplyDigit(uint64* input, uint64 digit, uint64* resultDigitsPlus1, int n, bool isSigned = false)
	{
		makeZero(resultDigitsPlus1, n);
		resultDigitsPlus1[n] = MultiplyDigitAndAdd(input, digit, resultDigitsPlus1, n);
		if (!isSigned)
		{
			return;
		}
		if ((input[n - 1] & Bit63) != 0)
		{
			SubtractDigit(&resultDigitsPlus1[n], digit, 1);
		}
		if ((digit & Bit63) != 0)
		{
			Subtract(&resultDigitsPlus1[1], input, &resultDigitsPlus1[1], 0, n);
		}
	}

	X64_API uint64 __cdecl DivideDigit(uint64* result, uint64 digit, int n, bool isSigned = false)
	{
		if (!isUdivmod128FunctionProtected)
		{
			DWORD dummy;
			VirtualProtect(udivmod128Data, sizeof(udivmod128Data), PAGE_EXECUTE_READWRITE, &dummy);
			isUdivmod128FunctionProtected = true;
		}

		bool negative = false;
		if (isSigned)
		{
			if ((digit & Bit63) != 0)
			{
				digit = (uint64)(-(__int64)digit);
				negative = true;
			}
			if (n > 0 && (result[n - 1] & Bit63) != 0)
			{
				negate(result, n);
				negative = !negative;
			}
		}
		uint64 mod = 0;
		for (int i = n; --i >= 0; )
		{
			mod = udivmod128(&result[i], mod, digit);
		}
		if (negative)
		{
			mod = (uint64)(-(__int64)mod);
			negate(result, n);
		}
		return mod;
	}

	X64_API uint8 __cdecl SubtractDigit(uint64* result, uint64 digit, int n)
	{
		uint8 borrow = _subborrow_u64(0, result[0], digit, result);
		n--;
		result++;
		while (--n >= 0 && borrow != 0)
		{
			borrow = _subborrow_u64(borrow, result[0], 0, result);
			result++;
		}
		return borrow;
	}

	X64_API uint8 __cdecl Add(uint64* input1, uint64* input2, uint64* result, uint8 carry, int n)
	{
		for (int i = 0; i < n; i++)
		{
			carry = _addcarry_u64(carry, input1[i], input2[i], &result[i]);
		}
		return carry;
	}

	X64_API uint8 __cdecl Subtract(uint64* input1, uint64* input2, uint64* result, uint8 borrow, int n)
	{
		for (int i = 0; i < n; i++)
		{
			borrow = _subborrow_u64(borrow, input1[i], input2[i], &result[i]);
		}
		return borrow;
	}

	///result += input * digit
	X64_API uint64 __cdecl MultiplyDigitAndAdd(uint64* input, uint64 digit, uint64* result, int n)
	{
		uint64 carry = 0;
		for (int i = 0; i < n; i++)
		{
			uint64 hi, lo = _umul128(input[i], digit, &hi);
			_addcarry_u64(_addcarry_u64(0, lo, carry, &lo), hi, 0, &hi);
			_addcarry_u64(_addcarry_u64(0, result[i], lo, &result[i]), hi, 0, &carry);
		}
		return carry;
	}

	///result += input * digit
	X64_API uint64 __cdecl MultiplyDigitAndSubtract(uint64* input, uint64 digit, uint64* result, int n)
	{
		uint64 borrow = 0;
		for (int i = 0; i < n; i++)
		{
			uint64 hi, lo = _umul128(input[i], digit, &hi);
			_addcarry_u64(_addcarry_u64(0, lo, borrow, &lo), hi, 0, &hi);
			_addcarry_u64(_subborrow_u64(0, result[i], lo, &result[i]), hi, 0, &borrow);
		}
		return borrow;
	}

	//input 1 and 2 of size n; result of size n * 2 qwords.
	X64_API void __cdecl Multiply(uint64* input1, uint64* input2, uint64* result, int n)
	{
		for (int i = n; --i >= 0;)
		{
			result[i] = 0;
		}
		for (int i = 0; i < n; i++)
		{
			result[n + i] = MultiplyDigitAndAdd(input1, input2[i], &result[i], n);
		}
	}

	X64_API void __cdecl Square(uint64* input, uint64* result, int n)
	{
		uint64    carry;
		result[0] = carry = 0;
		for (int i = 0; i < n; i++)
		{
			carry = getSquareSum(&input[i], input, result++, carry);
		}
		for (int i = 1; i < n; i++)
		{
			carry = getSquareSum(&input[n - 1], &input[i], result++, carry);
		}
	}

	//input 1 and 2 of size n; result of size n * 2 qwords.
	X64_API void __cdecl Karatsuba(uint64* input1, uint64* input2, uint64* result, int n, uint64* temporaryBuffer = NULL)
	{
		bool mustDisposeTemporaryBuffer = false;
		if (temporaryBuffer == NULL)
		{
			int tempSize = GetKaratsubaMultiplicationBufferSize(n);
			temporaryBuffer = new uint64[tempSize];
			mustDisposeTemporaryBuffer = true;
		}
		karatsuba(input1, input2, result, n, temporaryBuffer);
		if (mustDisposeTemporaryBuffer)
		{
			delete[] temporaryBuffer;
		}
	}

	//input 1 and 2 of size n; result of size n * 2 qwords.
	X64_API void __cdecl KaratsubaSquare(uint64* input, uint64* result, int n, uint64* temporaryBuffer = NULL)
	{
		bool mustDisposeTemporaryBuffer = false;
		if (temporaryBuffer == NULL)
		{
			int tempSize = GetKaratsubaSquareBufferSize(n);
			temporaryBuffer = new uint64[tempSize];
			mustDisposeTemporaryBuffer = true;
		}
		karatsubaSquare(input, result, n, temporaryBuffer);
		if (mustDisposeTemporaryBuffer)
		{
			delete[] temporaryBuffer;
		}
	}

	static void montgomeryMultiplyAndReduce(uint64* number1, uint64* number2, uint64* mod, uint64 inv_mod_0, uint64* temporary2x, uint64* destination, int n, uint64* temporaryBuffer)
	{
		Karatsuba(number1, number2, temporary2x, n, temporaryBuffer);
		montgomery_reduction(temporary2x, mod, inv_mod_0, n);
		let(destination, temporary2x, n);
	}

	static void montgomerySquareAndReduce(uint64* sourceAndDestination, uint64* mod, uint64 inv_mod_0, uint64* temporary2x, int n, uint64* temporaryBuffer)
	{
		KaratsubaSquare(sourceAndDestination, temporary2x, n, temporaryBuffer);
		montgomery_reduction(temporary2x, mod, inv_mod_0, n);
		let(sourceAndDestination, temporary2x, n);
	}

	//operation: base = base^power modulo 'mod' using montgomery reduction with exp_cache
	//RESTRICTION: 'mod' must have lowest and highest bits set (for RSA there are)
	X64_API void __cdecl MontgomeryExpMod(uint64* base, uint64* power, uint64* mod, int n, int cacheBits)
	{
		int hb = high_bit(power, n);
		if (hb < 0)
		{
			makeZero(base, n);
			base[0] = 1;
			return;
		}

		int temporaryKaratsubaSize = GetKaratsubaMultiplicationBufferSize(n);

		int cacheN = 1 << cacheBits;
		uint64* temporaryBuffer = new uint64[n * 4 + temporaryKaratsubaSize + cacheN * n];
		uint64
			*r1 = &temporaryBuffer[0],
			*r2 = &temporaryBuffer[n * 2],
			*temporaryKaratsubaBuffer = &temporaryBuffer[n * 4],
			inv_mod_0 = invert(mod[0]);
		uint64* cache = &temporaryBuffer[n * 4 + temporaryKaratsubaSize];

		makeZero(r1, n);
		let(&r1[n], base, n);
		int_mod(r1, n * 2, mod, n);
		let(&cache[0], r1, n);
		montgomerySquareAndReduce(r1, mod, inv_mod_0, r2, n, temporaryKaratsubaBuffer);

		for (int entry = 1; entry < cacheN; entry++)
		{
			montgomeryMultiplyAndReduce(&cache[(entry - 1) * n], r1, mod, inv_mod_0, r2, &cache[entry * n], n, temporaryKaratsubaBuffer);
		}

		int endPosition = hb;
		int startPosition = chooseStartPosition(power, hb, cacheBits);
		uint64 *source = &cache[getBitsValue(power, startPosition + 1, endPosition) * n];
		let(r1, source, n);
		for (int i = startPosition; --i >= 0;)
		{
			if ((power[i >> 6] & (1ULL << (i & 63))) != 0)
			{	//bit(power, i) is set
				startPosition = chooseStartPosition(power, i, cacheBits);
				for (int count = i - startPosition + 1; --count >= 0;)
				{
					montgomerySquareAndReduce(r1, mod, inv_mod_0, r2, n, temporaryKaratsubaBuffer);
					//r1 = r1^2 % mod
				}
				montgomeryMultiplyAndReduce(r1, &cache[getBitsValue(power, startPosition + 1, i) * n], mod, inv_mod_0, r2, r1, n, temporaryKaratsubaBuffer);
				//r1= r1 * cache % mod
				i = startPosition;
				continue;
			}
			montgomerySquareAndReduce(r1, mod, inv_mod_0, r2, n, temporaryKaratsubaBuffer);
			//r1 = r1^2 % mod
		}

		makeZero(&r1[n], n);
		montgomery_reduction(r1, mod, inv_mod_0, n);
		let(base, r1, n);

		delete[] temporaryBuffer;
	}

	static void multiplyAndReduce(uint64* number1, uint64* number2, uint64* mod, uint64* temporary2x, uint64* destination, int n, uint64* temporaryBuffer)
	{
		Karatsuba(number1, number2, temporary2x, n, temporaryBuffer);
		int_mod(temporary2x, n * 2, mod, n);
		let(destination, temporary2x, n);
	}

	static void squareAndReduce(uint64* sourceAndDestination, uint64* mod, uint64* temporary2x, int n, uint64* temporaryBuffer)
	{
		KaratsubaSquare(sourceAndDestination, temporary2x, n, temporaryBuffer);
		int_mod(temporary2x, n * 2, mod, n);
		let(sourceAndDestination, temporary2x, n);
	}

	X64_API void __cdecl ExpMod(uint64* base, uint64* power, uint64* mod, int n, int cacheBits)
	{
		int     hb = high_bit(power, n);
		if (hb < 0)
		{
			makeZero(base, n);
			base[0] = 1;
			return;
		}

		int temporaryKaratsubaSize = GetKaratsubaMultiplicationBufferSize(n);

		int cacheN = 1 << cacheBits;
		uint64* temporaryBuffer = new uint64[n * 4 + temporaryKaratsubaSize + cacheN * n];
		uint64
			*r1 = &temporaryBuffer[0],
			*r2 = &temporaryBuffer[n * 2],
			*temporaryKaratsubaBuffer = &temporaryBuffer[n * 4],
			*cache = &temporaryBuffer[n * 4 + temporaryKaratsubaSize];

		let(&cache[0], base, n);
		let(r1, base, n);
		squareAndReduce(r1, mod, r2, n, temporaryKaratsubaBuffer);	//r1 = r1^2 % mod
		for (int entry = 1; entry < cacheN; entry++)
		{
			multiplyAndReduce(&cache[(entry - 1) * n], r1, mod, r2, &cache[entry * n], n, temporaryKaratsubaBuffer);
		}

		int endPosition = hb;
		int startPosition = chooseStartPosition(power, hb, cacheBits);
		uint64 *source = &cache[getBitsValue(power, startPosition + 1, endPosition) * n];
		let(r1, source, n);
		for (int i = startPosition; --i >= 0;)
		{
			if ((power[i >> 6] & (1ULL << (i & 63))) != 0)
			{	//bit(power, i) is set
				startPosition = chooseStartPosition(power, i, cacheBits);
				for (int count = i - startPosition + 1; --count >= 0;)
				{
					squareAndReduce(r1, mod, r2, n, temporaryKaratsubaBuffer);
					//r1 = r1^2 % mod
				}
				multiplyAndReduce(&cache[getBitsValue(power, startPosition + 1, i) * n], r1, mod, r2, r1, n, temporaryKaratsubaBuffer);
				//r1= r1 * cache % mod
				i = startPosition;
				continue;
			}
			squareAndReduce(r1, mod, r2, n, temporaryKaratsubaBuffer);
			//r1 = r1^2 % mod
		}

		let(base, r1, n);
		delete[] temporaryBuffer;
	}

	//return a = a / X 'modulo' P; use P = odd number.
	X64_API void __cdecl DivideModP(uint64* _A, uint64* _X, uint64* P, int n)
	{
		if (isZero(_X, n) || isZero(P, n))
		{
			makeZero(_A, n);
			return;
		}
		uint64
			*temporaryBuffer = new uint64[n * 7 + 8],
			*x = &temporaryBuffer[0],
			*y = &temporaryBuffer[n],
			*c = &temporaryBuffer[n * 2],
			*t1 = &temporaryBuffer[n * 3],
			*t2 = &temporaryBuffer[n * 4 + 2],
			*t3 = &temporaryBuffer[n * 5 + 4],
			*t4 = &temporaryBuffer[n * 6 + 6],
			*a = _A,
			inverseP0 = invert(P[0]);
		//rule: a/_A*_X + b*P = x
		//rule: c/_A*_X + d*P = y
		let(x, _X, n);
		let(y, P, n);
		makeZero(c, n);

		{
			int shift = low_bit(x, n);
			int_shr(x, shift, n);
			shrModP(a, n, shift, P, n, inverseP0);
		}
		{
			int shift = low_bit(y, n);
			int_shr(y, shift, n);
			//shrModP(c, n, shift, P, n, inverseP0);
		}

		bool isZero1 = false, isZero2 = false;
		while (!isZero1 && !isZero2)
		{
			int hi = maximum(high_bit(x, n), high_bit(y, n));
			uint64 ra[2], rb[2], rc[2], rd[2];
			uint64 hix = getBitsValue(x, hi - 63, hi);
			uint64 hiy = getBitsValue(y, hi - 63, hi);
			estimateBinaryShift(hix, hiy, x[0], y[0], ra, rb, rc, rd);
			//x <- (x * a + y * b) >> shr
			//y <- (x * c + y * d) >> shr
			makeZero(t1, n + 2);
			makeZero(t2, n + 2);
			makeZero(t3, n + 2);
			makeZero(t4, n + 2);
			getSummation(x, ra, t1, n);
			getSummation(y, rb, t1, n);
			getSummation(x, rc, t2, n);
			getSummation(y, rd, t2, n);
			getSummation(a, ra, t3, n);
			getSummation(c, rb, t3, n);
			getSummation(a, rc, t4, n);
			getSummation(c, rd, t4, n);

			if (((sint64)t1[n + 1]) < 0)
			{
				negate(t1, n + 2);
				negate(t3, n + 2);
			}
			if (((sint64)t2[n + 1]) < 0)
			{
				negate(t2, n + 2);
				negate(t4, n + 2);
			}

			isZero1 = isZero(t1, n + 2);
			isZero2 = isZero(t2, n + 2);

			if (!isZero1)
			{
				int shift = low_bit(t1, n + 2);
				int_shr(t1, shift, n + 2);
				shrModP(t3, n + 2, shift, P, n, inverseP0);
			}
			let(x, t1, n);
			let(a, t3, n);

			if (!isZero2)
			{
				int shift = low_bit(t2, n + 2);
				int_shr(t2, shift, n + 2);
				shrModP(t4, n + 2, shift, P, n, inverseP0);
			}
			let(y, t2, n);
			let(c, t4, n);
		}

		let(_A, isZero1 ? c : a, n);
		delete[] temporaryBuffer;
	}

#define	swap_pointers(first, second, temp)	uint64* temp = first; first = second; second = temp
	//return a = a / X 'modulo' P; use P = odd number.
	X64_API void __cdecl DivideModPSlow(uint64* _A, uint64* _X, uint64* P, int n)
	{
		uint64
			*temporaryBuffer = new uint64[n * 4],
			*x = &temporaryBuffer[0],
			*y = &temporaryBuffer[n],
			*c = &temporaryBuffer[n * 2],
			*t = &temporaryBuffer[n * 3],
			*a = _A,
			inverseP0 = invert(P[0]);
		//rule: a/_A*_X + b*P = x
		//rule: c/_A*_X + d*P = y
		let(x, _X, n);
		let(y, P, n);
		makeZero(c, n);
		goto start;

		while (true)
		{
			if (Subtract(x, y, t, 0, n))
			{
				negate(t, n);
				swap_pointers(x, y, auxiliary1);
				swap_pointers(a, c, auxiliary2);
			}
			{	//let(x, t, n);
				swap_pointers(x, t, auxiliary3);
			}
			if (Subtract(a, c, a, 0, n))
			{
				Add(a, P, a, 0, n);     //A = A - C 'modulo' P
			}
		start:
			if (int_digits(x, n) <= 0)
			{
				break;
			}
			//when x is even it is corrected in the while below...
			int shift = low_bit(x, n);
			int_shr(x, shift, n);
			shrModP(a, n, shift, P, n, inverseP0);
		}
		if (_A != c)
		{
			let(_A, c, n);
		}
		delete[] temporaryBuffer;
	}

#undef	swap_pointers

	X64_API void __cdecl GetDivMod(uint64* input2x, uint64* mod, int n, uint64* quotient2x = NULL, bool isSigned = false)
	{
		if (!isSigned)
		{
			int_mod(input2x, n * 2, mod, n, quotient2x);
			return;
		}

		bool negativeInput = false;
		if (n > 0 && (input2x[n * 2 - 1] & Bit63) != 0)
		{
			negate(input2x, n * 2);
			negativeInput = true;
		}
		bool negativeMod = false;
		if (n > 0 && (mod[n - 1] & Bit63) != 0)
		{
			negate(mod, n);
			negativeMod = true;
		}
		int_mod(input2x, n * 2, mod, n, quotient2x);
		if (negativeMod)
		{
			negate(mod, n);
		}
		if (negativeInput != negativeMod)
		{
			if (quotient2x != NULL)
			{
				negate(quotient2x, n * 2);
			}
			negate(input2x, n * 2);
		}
	}

	//prime = 2 ^ (n * 64) - primeLowPart
	X64_API void __cdecl GetModFriendlyPrime(uint64* input2x, uint64 primeLowPart, int n)
	{
		uint64 carry = MultiplyDigitAndAdd(&input2x[n], primeLowPart, input2x, n);
		uint64 digit = MultiplyDigitAndAdd(&carry, primeLowPart, input2x, 1);
		uint8 overflow = AddDigit(&input2x[1], digit, n - 1);
		while (overflow != 0)
		{
			overflow = AddDigit(input2x, primeLowPart, n);
		}
	}

	X64_API void __cdecl RemodFriendlyPrime(uint64* input1x, uint64 primeLowPart, int n)
	{
		uint64 exactLowMaximum = 0ULL - primeLowPart;
		if (input1x[0] < exactLowMaximum)
		{
			return;	//already < P.
		}
		for (int i = n; --i > 0; )
		{
			if (input1x[i] != MaxValueUint64)
			{
				return;	//already < P.
			}
		}
		for (int i = n; --i > 0; )
		{
			input1x[i] = 0;
		}
		input1x[0] -= exactLowMaximum;
	}

	X64_API void __cdecl GetModFriendlyPrime521(uint64* input2x, int n)
	{
		uint64 carry = add_shl(input2x, &input2x[n], 55, n);
		//carry has to be added at bit position: 55.
		uint64 adder[2] = { (input2x[n - 1] >> 9) | (carry << 55), carry >> 9 };
		input2x[n - 1] &= 0x01FFULL;
		uint8 bitCarry = Add(input2x, adder, input2x, 0, 2);
		if (bitCarry != 0)
		{
			AddDigit(&input2x[2], 1ULL, n - 2);
		}
	}

	X64_API void __cdecl RemodFriendlyPrime521(uint64* input1x, int n)
	{
		while (input1x[n - 1] >= 0x0200ULL)
		{
			uint64 digit = input1x[n - 1] >> 9;
			input1x[n - 1] &= 0x01FFULL;
			AddDigit(input1x, digit, n);
		}
		if (input1x[n - 1] != 0x01FFULL)
		{
			return;
		}
		for (int i = n - 1; --i >= 0; )
		{
			if (input1x[i] != MaxValueUint64)
			{
				return;
			}
		}
		for (int i = n; --i >= 0; )
		{
			input1x[i] = 0ULL;
		}
	}

	X64_API void __cdecl AddModP(uint64* input1, uint64* input2, uint64* destination, uint64 LowP, int n)
	{	//2^256 == LowP modulo 'P' because P = 2^256 - LowP
		uint8 carry = Add(input1, input2, destination, 0, n);
		while (carry != 0)
		{
			carry = AddDigit(destination, LowP, n);
		}
	}

	X64_API void __cdecl SubtractModP(uint64* input1, uint64* input2, uint64* destination, uint64 LowP, int n)
	{
		uint8 borrow = Subtract(input1, input2, destination, 0, n);
		while (borrow != 0)
		{
			borrow = SubtractDigit(destination, LowP, n);
		}
	}

	X64_API void __cdecl AddScaledModP(uint64* destination, uint64* input, uint64 multiplier, uint64 LowP, int n)
	{
		uint64 carry = MultiplyDigitAndAdd(input, multiplier, destination, n);
		carry = MultiplyDigitAndAdd(&LowP, carry, destination, 1);
		if (carry != 0)
		{
			carry = AddDigit(&destination[1], carry, n - 1);
		}
		while (carry != 0)
		{
			carry = AddDigit(destination, LowP, n);
		}
	}

	X64_API void __cdecl SubScaledModP(uint64* destination, uint64* input, uint64 multiplier, uint64 LowP, int n)
	{
		uint64 borrow = MultiplyDigitAndSubtract(input, multiplier, destination, n);
		borrow = MultiplyDigitAndSubtract(&LowP, borrow, destination, 1);
		if (borrow != 0)
		{
			borrow = SubtractDigit(&destination[1], borrow, n - 1);
		}
		while (borrow != 0)
		{
			borrow = SubtractDigit(destination, LowP, n);
		}
	}

	X64_API void __cdecl AddModP521(uint64* input1, uint64* input2, uint64* destination, int n)
	{	//2^256 == LowP modulo 'P' because P = 2^256 - LowP
		uint8 carry = Add(input1, input2, destination, 0, n);
		if (carry != 0)
		{
			uint64 digit = (1ULL << (64 - 9)) | (destination[n - 1] >> 9);
			destination[n - 1] &= 511ULL;
			AddDigit(destination, digit, n);
		}
	}

	X64_API void __cdecl SubtractModP521(uint64* input1, uint64* input2, uint64* destination, int n)
	{
		uint8 borrow = Subtract(input1, input2, destination, 0, n);
		while (borrow != 0)
		{
			uint64 digit = (1ULL << (64 - 9)) - (destination[n - 1] >> 9);
			destination[n - 1] &= 511ULL;
			borrow = SubtractDigit(destination, digit, n);
		}
	}

	X64_API void __cdecl AddScaledModP521(uint64* destination, uint64* input, uint64 multiplier, int n)
	{
		uint64 carry = MultiplyDigitAndAdd(input, multiplier, destination, n);
		if (carry == 0)
		{
			return;
		}
		if (carry < 512ULL)
		{
			carry = (destination[n - 1] >> 9) | (carry << (64 - 9));
			destination[n - 1] &= 511ULL;
			AddDigit(destination, carry, n);
			return;
		}
		uint64 specialCarry[2] = { (destination[n - 1] >> 9) | (carry << (64 - 9)), carry >> 9 };
		destination[n - 1] &= 511ULL;
		if (Add(destination, specialCarry, destination, 0, 2) != 0)
		{
			AddDigit(&destination[2], 1ULL, n - 2);
		}
	}

	X64_API void __cdecl SubScaledModP521(uint64* destination, uint64* input, uint64 multiplier, int n)
	{
		uint64 borrow = MultiplyDigitAndSubtract(input, multiplier, destination, n);
		if (borrow == 0)
		{
			return;
		}
		uint8 bitBorrow;
		if (borrow < 512ULL)
		{
			borrow = (borrow << (64 - 9)) - (destination[n - 1] >> 9);
			destination[n - 1] &= 511ULL;
			bitBorrow = SubtractDigit(destination, borrow, n);
		}
		else
		{
			uint64 specialCarry[2] = { (borrow << (64 - 9)) - (destination[n - 1] >> 9), borrow >> 9 };
			destination[n - 1] &= 511ULL;
			if (Subtract(destination, specialCarry, destination, 0, 2) == 0)
			{
				bitBorrow = 0;
			}
			else
			{
				bitBorrow = SubtractDigit(&destination[2], 1ULL, n - 2);
			}
		}
		while (bitBorrow != 0)
		{
			uint64 digit = (1ULL << (64 - 9)) - (destination[n - 1] >> 9);
			destination[n - 1] &= 511ULL;
			bitBorrow = SubtractDigit(destination, digit, n);
		}
	}
};
