#include "BigInt.h"

BigInt::BigInt() :nInt(0), nBit(0){
}
BigInt::BigInt(int _nBit) : nBit(_nBit){
	nInt = (int)ceil(1.0*_nBit / INTSZ);
	ints.resize(nInt);
}
BigInt::BigInt(const BigInt& bint){
	nInt = bint.nInt;
	nBit = bint.nBit;
	ints = bint.ints;
}
BigInt::~BigInt(){}

BigInt BigInt::operator& (const BigInt &bint){
	BigInt newBInt(*this);
	for (int i = 0; i < nInt; ++i){
		newBInt.ints[i] &= bint.ints[i];
	}
	return newBInt;
}
BigInt BigInt::operator^ (const BigInt &bint){
	BigInt newBInt(*this);
	for (int i = 0; i < nInt; ++i){
		newBInt.ints[i] ^= bint.ints[i];
	}
	return newBInt;
}
BigInt BigInt::operator| (const BigInt &bint){
	BigInt newBInt(*this);
	for (int i = 0; i < nInt; ++i){
		newBInt.ints[i] |= bint.ints[i];
	}
	return newBInt;
}
void BigInt::operator|= (const BigInt &bint){
	for (int i = 0; i < nInt; ++i){
		ints[i] |= bint.ints[i];
	}
}
void BigInt::operator-= (const BigInt &bint){
	for (int i = 0; i < nInt; ++i){
		ints[i] = ints[i] ^ (ints[i] & bint.ints[i]);
	}
}
BigInt & BigInt::operator=(const BigInt &bint){
	nBit = bint.nBit;
	nInt = bint.nInt;
	ints = bint.ints;
	return *this;
}
bool BigInt::operator==(const BigInt &bint) const {
	return (nBit == bint.nBit) && (ints == bint.ints);
}
