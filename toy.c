#include "toy.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


// This function fills a buffer with small random values
// The values are either 0, 1, or -1 (represented as q-1)
// The distribution is approximately binomial
// Inputs: buf - a pointer to a buffer of type short
//         n - the size of the buffer
static void toy_fill_small(short *buf, int n)
{
#if 1
	for (int k=0;k<n;++k)
	{
		short val=rand();
		val=(val>>1&1)-(val&1);
		if(val<0)
			{
				val+=TK_Q;
			}
		buf[k]=val;
	}

#else 
	memset(buf, 0, n*sizeof(short));
#endif
}

// This function performs a naive polynomial multiplication modulo q
// It assumes that the polynomials have degree N-1 and coefficients in Z_q
// It uses a simple schoolbook method with O(N^2) complexity
// Inputs: dst - a pointer to a buffer of type short for the result
//         a - a pointer to a buffer of type short for the first polynomial
//         b - a pointer to a buffer of type short for the second polynomial
//         add - a flag indicating whether to add or overwrite the result

static void toy_polmul_naive(short *dst, const short *a, const short *b, int add)
{

	dst[0]=((dst[0]&-add)+a[0]*b[0]+NEG(a[3])*b[1]+NEG(a[2] )*b[2]+NEG(a[1])*b[3])%TK_Q;
	
	dst[1]=((dst[1]&-add)+a[1]*b[0]+a[0]*b[1]+NEG(a[3])*b[2]+NEG(a[2])*b[3])%TK_Q;
	
	dst[2]=((dst[2]&-add)+a[2]*b[0]+a[1]*b[1]+a[0]*b[2]+NEG(a[3])*b[3])%TK_Q;
	
	dst[3]=((dst[3]&-add)+a[3]*b[0]+a[2]*b[1]+a[1]*b[2]+a[0]*b[3])%TK_Q;
}

// This function performs a matrix-vector multiplication modulo q
// It assumes that the matrix has size K x K and the vector has size K
// It uses the toy_polmul_naive function for each entry of the result
// Inputs: dst - a pointer to a buffer of type short for the result
//         mat - a pointer to a buffer of type short for the matrix
//         vec - a pointer to a buffer of type short for the vector

static void toy_mulmv(short *dst, const short *mat, const short *vec)
{

memset(dst, 0, TK_K*TK_N*sizeof(short));

for(int kv=0, idx=0;kv<TK_K*TK_N; kv+=TK_N)

{
	for(int k=0;k<TK_K*TK_N;k+=TK_N, idx+=TK_N)
	toy_polmul_naive(dst+kv, mat+idx, vec+k, 1);
}
}

// This function performs a transposed matrix-vector multiplication modulo q
// It assumes that the matrix has size K x K and the vector has size K
// It uses the toy_polmul_naive function for each entry of the result
// Inputs: dst - a pointer to a buffer of type short for the result
//         mat - a pointer to a buffer of type short for the matrix
//         vec - a pointer to a buffer of type short for the vector

static void toy_mulmTv(short *dst, const short *mat, const short *vec)
{
    memset(dst, 0, TK_K*TK_N*sizeof(short));
    
    for(int kv=0; kv<TK_K*TK_N; kv+=TK_N)
    {
        for(int k=0; k<TK_K*TK_N; k+=TK_N)
        {
            toy_polmul_naive(dst+kv, mat+TK_K*k+kv, vec+k, 1);
        }
    }
}


// This function performs a dot product of two vectors modulo q
// It assumes that the vectors have size K and each entry is a/ polynomial of degree N-1
// It uses the toy_polmul_naive function for each entry of the result
// Inputs: dst - a pointer to a buffer of type short for the result
//         v1 - a pointer to a buffer of type short for the first vector
//         v2 - a pointer to a buffer of type short for the second vector

static void toy_dot(short *dst, const short *v1, const short *v2)
{
    memset(dst, 0, TK_N*sizeof(short));

    for(int k=0; k<TK_K*TK_N; k+=TK_N)
    {
        toy_polmul_naive(dst, v1+k, v2+k, 1);
    }
}

// This function performs a vector addition or subtraction modulo q
// It assumes that the vectors have the same size and each entry is an integer in Z_q
// Inputs: dst - a pointer to a buffer of type short for the result
//         v1 - a pointer to a buffer of type short for the first vector
//         v2 - a pointer to a buffer of type short for the second vector
//         count - the size of the vectors
//         v2_neg - a flag indicating whether to subtract or add v2

static void toy_add(short *dst, const short *v1, const short *v2,
                    int count, int v2_neg)
{
    for(int k=0; k<count; ++k)
    {
        short val = v2[k];
        
        if(v2_neg)
            val = NEG(val);
        dst[k] = (v1[k] + val) % TK_Q;
    }
}

// This function generates the public and private keys for the toy NTRU scheme
// It uses the toy_fill_small and toy_mulmv functions to generate random values
// Inputs: A - a pointer to a buffer of type short for the public matrix
//         t - a pointer to a buffer of type short for the public vector
//         s - a pointer to a buffer of type short for the private vector

void toy_gen(short *A, short *t, short *s)
{
    short e[TK_K*TK_N];
    
    // Fill A with random values modulo q
    for(int k = 0; k < TK_K*TK_K*TK_N; ++k)
    {
        A[k] = rand() % TK_Q;
    }
    
    // Fill s and e with small random values
    toy_fill_small(s, TK_K*TK_N);
    toy_fill_small(e, TK_K*TK_N);
    
    // Compute t = A * s + e
    toy_mulmv(t, A, s);  
    toy_add(t, t, e, TK_K*TK_N, 0);  
}

// This function encrypts a plaintext message using the toy NTRU scheme
// It uses the toy_fill_small, toy_mulmTv, toy_dot, and toy_add functions to generate random values and perform the encryption
// Inputs: A - a pointer to a buffer of type short for the public matrix
//         t - a pointer to a buffer of type short for the public vector
//         plain - an integer representing the plaintext message as a binary vector of size N
//         u - a pointer to a buffer of type short for the first part of the ciphertext
//         v - a pointer to a buffer of type short for the second part of the ciphertext

void toy_enc(const short *A, const short *t, int plain, short *u, short *v)
{
    short r[TK_K*TK_N], e1[TK_K*TK_N], e2[TK_N];
    
    // Fill r, e1, and e2 with small random values
    toy_fill_small(r, TK_K*TK_N);
    toy_fill_small(e1, TK_K*TK_N);
    toy_fill_small(e2, TK_N);

    // Compute u = A^T * r + e1
    toy_mulmTv(u, A, r);
	toy_add(u,u,e1,TK_K*TK_N,0);
  
    // Compute v = t * r + e2 + m * (q/2)
    toy_dot(v,t,r);
    toy_add(v,v,e2,TK_N,0); //noise vectors for more security
    
    // Add the plaintext message to v by shifting it by q/2
    for(int k=0; k<TK_N; ++k)
        v[k] = (v[k] + ((TK_Q>>1) & -(plain>>(TK_N-1-k)&1))) % TK_Q;
}

// This function decrypts a ciphertext using the toy NTRU scheme
// It uses the toy_dot and toy_add functions to perform the decryption
// Inputs: s - a pointer to a buffer of type short for the private vector
//         u - a pointer to a buffer of type short for the first part of the ciphertext
//         v - a pointer to a buffer of type short for the second part of the ciphertext

int toy_dec(const short *s, const short *u, const short *v) 
{
    short p[TK_N], plain;
    
    // Compute p = v - s * u
    toy_dot(p, s, u);
    toy_add(p, v, p, TK_N, 1);

    plain = 0;
    
    for (int k = 0; k < TK_N; ++k) 
    {
        int val = p[k];
        
        // Adjust the value if it is negative
        if (val > TK_Q / 2) 
        {
            val -= TK_Q;
            
        }
        printf("%5d", val);
        
        // Determine the bit by comparing the absolute value to q/4
        int bit = abs(val) > TK_Q / 4;
      

        // Set the bit in the plaintext message
        plain |= bit <<(TK_N-1-k);
    }

    return plain;
}

// This function tests the toy NTRU scheme by encrypting and decrypting some messages
// It uses the toy_gen, toy_enc, and toy_dec functions
// Inputs: none

void toy_test()
{
#if 1

		short A[TK_K*TK_K*TK_N], s[TK_K*TK_N], e[TK_K*TK_N],
			r[TK_K*TK_N], e1[TK_K*TK_N], e2[TK_N];

// Fill A with random values modulo q
for(int k=0;k<TK_K*TK_K*TK_N;++k)
		A[k]=rand()%TK_Q;
//memset(A, 9, TK_K*TK_K*TK_N*sizeof(short));

// Fill s, e, r, e1, and e2 with small random values
toy_fill_small(s, TK_K*TK_N); 
toy_fill_small(e, TK_K*TK_N);
toy_fill_small(r, TK_K*TK_N);
toy_fill_small(e1, TK_K*TK_N);
toy_fill_small(e2, TK_N);

short temp[TK_K*TK_N];

// Compute temp = A * s + e
toy_mulmv(temp, A, s);
toy_add(temp, temp, e, TK_K*TK_N, 0);

short acc[TK_N];

// Compute acc = temp * r + e2
toy_dot(acc, temp, r);
toy_add(acc, acc, e2, TK_N, 0);

// Compute temp = A^T * r + e1
toy_mulmTv(temp, A, r);
toy_add(temp, temp, e1, TK_K*TK_N, 0);

short temp2[TK_N];

// Compute temp2 = s * temp
toy_dot(temp2, s, temp);

// Compute acc = acc - temp2
toy_add(acc, acc, temp2, TK_N, 1);

// Check if the coefficients of acc are less than q/4
for(int k=0;k<TK_N;++k)//coeffs must be < q/4

{
int val=acc[k]; 
if (val>TK_Q/2)
		val-=TK_Q;
if(abs(val)>TK_Q/4)
printf("ERROR %5d\n", val);
}

#endif
}

// This function prints a binary representation of an integer
// It assumes that the integer has at most nbits bits
// Inputs: x - an integer to be printed
//         nbits - the number of bits to be printed

