# CSE451s:Computer and Network Security Assignment 2 - Toy Kyber Cryptosystem


## Introduction

The Kyber cryptosystem is a lattice-based public-key encryption scheme that is designed to resist quantum attacks. It is based on the hardness of the learning with errors (LWE) problem over polynomial rings. The main idea is to encrypt a message by adding a small error to a linear combination of public polynomials, and decrypt it by using a secret polynomial to remove the error.

The toy-version of the Kyber cryptosystem is a simplified version that uses smaller parameters and simpler operations. It is intended for educational purposes only and should not be used for any security-critical applications.

## Parameters

The toy-version of the Kyber cryptosystem uses the following parameters:

- q = 97, a prime number that is the modulus of the polynomial ring
- N = 4, a power of two that is the degree of the polynomials
- K = 3, the dimension of the matrix and the vectors
- q/2 = 48, the midpoint of the modulus that is used to encode and decode binary messages

## Key Generation

The key generation algorithm generates a public key and a private key as follows:




- The public key consists of a matrix A of size K x K and a vector t of size K, where each entry is a polynomial of degree N-1 with coefficients in Z_q
- The private key consists of a vector s of size K, where each entry is a polynomial of degree N-1 with coefficients in {-1, 0, 1}
- The matrix A is filled with random polynomials modulo q
- The vector s and a vector e of size K are filled with small random polynomials with coefficients in {-1, 0, 1}
- The vector t is computed as t = A * s + e, where * denotes the matrix-vector multiplication modulo q

## Encryption

The encryption algorithm encrypts a plaintext message m, which is an integer representing a binary vector of size N, as follows:

- The ciphertext consists of a pair (u, v), where u and v are vectors of size K and N respectively, and each entry is a polynomial of degree N-1 with coefficients in Z_q
- A vector r of size K is filled with small random polynomials with coefficients in {-1, 0, 1}
- Two vectors e1 and e2 of size K and N respectively are filled with small random polynomials with coefficients in {-1, 0, 1}
- The vector u is computed as u = A^T * r + e1, where A^T denotes the transpose of the matrix A and * denotes the matrix-vector multiplication modulo q
- The vector v is computed as v = t * r + e2 + m * (q/2), where * denotes the dot product of two vectors modulo q and m * (q/2) denotes the addition of the message to v by shifting it by q/2

## Decryption

The decryption algorithm decrypts a ciphertext (u, v) as follows:

- The plaintext message m is recovered as an integer representing a binary vector of size N
- A vector p of size N is computed as p = v - s * u, where * denotes the dot product of two vectors modulo q
- The vector p is rounded to the nearest multiple of q/2 and then reduced modulo 2 to obtain the message m

## Compilation

You can use the gcc compiler to compile  program as follows:

```bash
gcc -o main toy.c main.c

## Execution
Run your program with the following command:

`./main`
