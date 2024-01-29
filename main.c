
#include "toy.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// This function prints a binary representation of an integer
// It assumes that the integer has at most nbits bits
// Inputs: x - an integer to be printed
//         nbits - the number of bits to be printed

static void print_bin(int x, int nbits)
{
for(int k=0;k<nbits;++k)
    printf("%d", x>>(nbits-1-k)&1);
}

int main(int argc, char **argv) {
    
short A[TK_K*TK_K*TK_N], t[TK_K*TK_N];
short s[TK_K*TK_N];

// Generate the keys
toy_gen(A, t, s); 

// Encrypt and decrypt some messages
for(int msg=0;msg<16;++msg)
{

    short u[TK_K*TK_N], v[TK_N];//ciphertext
    toy_enc(A, t, msg, u, v);
    short m2=toy_dec(s, u, v);
    printf(" %2d %2d ", msg, m2);
    print_bin(msg, TK_N);
    printf("  ");
    print_bin(m2, TK_N);
    printf("  ");
    print_bin(msg^m2, TK_N);
    printf("\n");

    }

    printf("Done.\n");
    
    return 0;
}
