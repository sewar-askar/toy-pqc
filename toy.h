
//header file for the toy Kyber cryptosystem
// It contains the declarations of the constants and the functions

#ifndef TOY_H
#define TOY_H

// Define the constants for the parameters
#define TK_K 3 // The dimension of the matrix A and the vectors s, r, e1
#define TK_N 4 // The degree of the polynomials
#define TK_Q 97 // The modulus of the ring
#define NEG(X) (TK_Q-(X)) // The negation of a value modulo q

// Declare the functions for the toy Kyber scheme
void toy_gen(short *A, short *t, short *s); // Key generation
void toy_enc(const short *A, const short *t, int plain, short *u, short *v); // Encryption
int toy_dec(const short *s, const short *u, const short *v); // Decryption



#endif