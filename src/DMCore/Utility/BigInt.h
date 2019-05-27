#ifndef _BIGINT_H_
#define _BIGINT_H_

#include <stdint.h>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define USEINT64
#ifdef USEINT64

#ifndef INT
#define INT uint64_t
#endif
#define INTSZ 64
#define RMASK 63 
#define DIVIER 6

#else

#ifndef INT
#define INT uint32_t
#endif
#define INTSZ 32
#define RMASK 31 
#define DIVIER 5

#endif

using namespace std;

class BigInt{

	inline void myitoa(INT num, char str[]);
	inline void reverse(char str[], int length);
public:
	vector<INT> ints;
	int nInt;
	int nBit;

	BigInt();
	BigInt(int _nBit);
	BigInt(const BigInt &bint);
	~BigInt();

	inline void set(int pos);
	inline void set(vector<int> & poss);
	inline void unset(int pos);
	inline bool get(int pos);
	inline vector<int> getOnesInd () const;
	inline int getOnesNum();
	inline void clear();
	inline bool isZero();
	inline bool isAllOne();
	inline bool shareBit(const BigInt &bint);
	inline int popOne();
	inline int find1st(int bit);
	inline void print();
	inline void flip();
	inline void resize(int _nBit);
	BigInt operator& (const BigInt &bint);
	BigInt operator^ (const BigInt &bint);
	BigInt operator| (const BigInt &bint);
	void operator|= (const BigInt &bint);
	void operator-= (const BigInt &bint);
	BigInt & operator= (const BigInt &bint);
	bool operator == (const BigInt &bint) const;

	inline void pack(int nbit_k, int nbit_S, int nbit_g,
		int n, vector<int> &Ss,
		int k, vector<int> &gs);
	inline void unpack(int nbit_k, int nbit_S, int nbit_g,
		int n, vector<int> &Ss,
		int &k, vector<int> &gs);
	inline void pack(int nbit_G, vector<int> &Gs);
	inline void unpack(int nbit_G, int n, vector<int> &Gs);
	inline void setINT(INT val, int nbit, int& inti, int& nExBit);
	inline void getINT(int &val, int nbit, INT mask, int& inti, int& nExBit);

	inline int size(){
		return nBit;
	}
};

class myHash{
public:
	inline size_t operator() (const BigInt &x) const {

		INT nHash = 0;
		for (INT i = 0; i < x.nInt; i++) {
			nHash += x.ints[i];
			nHash += (nHash << 10);
			nHash ^= (nHash >> 6);
		}
		nHash += (nHash << 3);
		nHash ^= (nHash >> 11);
		nHash += (nHash << 15);
		return (size_t)nHash;
	}
};

class myEq{
public:
	inline bool operator() (const BigInt &x, const BigInt &y) const {
		return (x == y);
	}
};

typedef unordered_map<BigInt, int, myHash, myEq> bitIntHash;

inline void BigInt::resize(int _nBit){
	nBit = _nBit;
	nInt = (int)ceil(1.0*_nBit / INTSZ);
	ints.resize(nInt);

	int move = (nBit - INTSZ * (nInt - 1));
	if(move < INTSZ){
		ints[nInt - 1] &= ((INT)1 << (nBit - INTSZ * (nInt - 1))) - 1;
	}
}

inline void BigInt::set(int pos){
	ints[pos >> DIVIER] |= ((INT)1 << (pos & RMASK));
}
inline void BigInt::set(vector<int> & poss){
	for (unsigned int i = 0; i < poss.size(); ++i){
		set(poss[i]);
	}
}
inline void BigInt::unset(int pos){
	ints[pos >> DIVIER] &= (~((INT)1 << (pos & RMASK)));
}
inline bool BigInt::get(int pos){
	return (ints[pos >> DIVIER] & ((INT)1 << (pos & RMASK))) != 0;
}
inline void BigInt::flip(){
	for (int i = 0; i < nInt; ++i){
		ints[i] = ~ints[i];
	}
	int move = (nBit - INTSZ * (nInt - 1));
	if(move < INTSZ){
		ints[nInt - 1] &= ((INT)1 << (nBit - INTSZ * (nInt - 1))) - 1;
	}
}
inline void BigInt::clear(){
	fill(ints.begin(), ints.end(), (INT)0);
}
inline bool BigInt::isZero(){
	for (int i = 0; i < nInt; ++i){
		if (ints[i] != 0)
			return false;
	}
	return true;
}
inline bool BigInt::isAllOne(){
	BigInt tmp = *this;
	tmp.flip();
	return tmp.isZero();
}
inline bool BigInt::shareBit(const BigInt &bint){
	for (int i = 0; i < nInt; ++i){
		if (ints[i] & bint.ints[i])
			return true;
	}
	return false;
}
inline int BigInt::find1st(int bit){
	BigInt tmp = *this;
	if (bit == 0) {
		tmp.flip();
	}
	return tmp.popOne();
}
inline int BigInt::popOne(){
	int index = -1;
	for (int i = 0; i < nInt; ++i){
		if (ints[i] != 0){
			index = (int)(log(ints[i] & ~(ints[i] - 1))/log(2));
			ints[i] &= ~((INT)1 << index);
			index += i * INTSZ;
			break;
		}
	}
	return index;
}
inline vector<int> BigInt::getOnesInd() const{
	vector<int> res;
	BigInt tmp = *this;
	int ind;
	while ((ind = tmp.popOne()) != -1){
		res.push_back(ind);
	}
	return res;
}
inline int BigInt::getOnesNum(){
	int res = 0;
	BigInt tmp = *this;
	int ind;
	while ((ind = tmp.popOne()) != -1){
		res++;
	}
	return res;
}
inline void BigInt::print(){
	cout << "\nnBit = " << nBit << "  nInt = " << nInt << endl;
	char buffer[INTSZ+1];
#ifdef USEINT64
	cout << "|| 321.987654321.987654321.987654321.987654321.987654321.987654321." << endl;
#else
	cout << "|| 1.987654321.987654321.987654321." << endl;
#endif
	for (int i = 0; i < nInt; ++i){
		myitoa(ints[i], buffer);
		cout << i << ": " << buffer << endl;
	}
	cout << endl;
}
inline void BigInt::reverse(char str[], int length) {
	int start = 0;
	int end = length - 1;
	while (start < end) {
		std::swap(*(str + start), *(str + end));
		start++;
		end--;
	}
}
inline void BigInt::myitoa(INT num, char str[]) {
	int ind = 0;
	for (ind = 0; ind<INTSZ; ++ind){
		str[ind] = '.';
	}
	str[ind] = '\0';

	ind = 0;
	while (num != 0) {
		INT rem = num % 2;
		str[ind++] = (rem > 0) ? '1' : '.';
		num = num >> 1;
	}
	reverse(str, INTSZ);
}
inline void BigInt::pack(int nbit_k, int nbit_S, int nbit_g, int n, vector<int> &Ss, int k, vector<int> &gs){
	ints[0] |= (INT)k;

	int nExBit = nbit_k;
	int inti = 0;
	for (int si : Ss){
		setINT((INT)si, nbit_S, inti, nExBit);
	}

	for (int gi : gs){
		setINT((INT)gi, nbit_g, inti, nExBit);
	}
}
inline void BigInt::pack(int nbit_G, vector<int> &Gs){
	int nExBit = 0;
	int inti = 0;
	for (int gi : Gs){
		setINT((INT)gi, nbit_G, inti, nExBit);
	}
}
inline void BigInt::setINT(INT val, int nbit, int &inti, int &nExBit){
	int nCurBit = nExBit + nbit;
	if (nCurBit <= INTSZ){
		ints[inti] |= (val << nExBit);
		nExBit = nCurBit;
		if (nCurBit == INTSZ){
			nExBit = 0;
			inti++;
		}
		else if (nCurBit > INTSZ){
			nExBit = nCurBit % INTSZ;
			inti += nCurBit / INTSZ;
		}
	} else {
		int nRight = INTSZ - nExBit;
		INT mask_r = ((INT)1 << nRight) - 1;
		ints[inti] |= (val & mask_r) << (nExBit);
		inti++;
		ints[inti] |= (val >> nRight);
		nExBit = nCurBit - INTSZ;
	}
}
inline void BigInt::getINT(int &val, int nbit, INT mask, int& inti, int& nExBit){
	int nCurBit = nExBit + nbit;
	val = 0;
	if (nCurBit <= INTSZ){
		val = (int)((ints[inti] >> nExBit) & mask);
		nExBit = nCurBit;
		if (nCurBit == INTSZ){
			nExBit = 0;
			inti++;
		}
		else if (nCurBit > INTSZ){
			nExBit = nCurBit % INTSZ;
			inti += nCurBit / INTSZ;
		}
	}
	else {
		int nRight = INTSZ - nExBit;
		int nLeft = nCurBit - INTSZ;
		INT mask_r = ((INT)1 << nRight) - 1;
		INT mask_l = ((INT)1 << nLeft) - 1;
		val |= (int)((ints[inti] >> nExBit) & mask_r);
		inti++;
		val |= (int)((ints[inti] & mask_l) << nRight);
		nExBit = nLeft;
	}
}
inline void BigInt::unpack(int nbit_k, int nbit_S, int nbit_g, int n, vector<int> &Ss, int &k, vector<int> &gs){

	INT mask_k = ((INT)1 << nbit_k) - 1;
	k = (int)(ints[0] & mask_k);

	int nExBit = nbit_k;
	int inti = 0;
	Ss.resize(n);
	INT mask_s = ((INT)1 << nbit_S) - 1;
	for (int si = 0; si < n; ++si){
		getINT(Ss[si], nbit_S, mask_s, inti, nExBit);
	}
	gs.resize(k);
	INT mask_g = ((INT)1 << nbit_g) - 1;
	for (int gi = 0; gi < k; ++gi){
		getINT(gs[gi], nbit_g, mask_g, inti, nExBit);
	}
}

inline void BigInt::unpack(int nbit_G, int n, vector<int> &Gs){
	int nExBit = 0;
	int inti = 0;
	Gs.resize(n);
	INT mask_s = ((INT)1 << nbit_G) - 1;
	for (int gi = 0; gi < n; ++gi){
		getINT(Gs[gi], nbit_G, mask_s, inti, nExBit);
	}
}



#endif //_BIGINT_H_
