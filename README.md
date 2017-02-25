# Elliptic curve study code 

This is my program to investigate elliptic curves.
It implements the elliptic curve group operation as described in
Understanding Cryptography by Christof Paar and Jan Pelzl.
It implements the double and add method for group multiplication.
By default, it uses the secp256k1 curve, but any elliptic curve can be
used by specifying the parameters in the command line arguments.

It uses the GMP library to deal with large numbers.

To compile: 

    g++ -O2 -std=c++11 -o ec ec.cpp -lgmpxx -lgmp

Usage: 

    ./ec [OPTIONS] <private key>
	    
    The elliptic curve is y^2 = x^3 + ax + b
        and generating point G modulo p.

    All numbers should be specified in hex, including the private key.

    The options are:
        -v      verbose print
        -a #    specify a in hex
        -b #    specify b in hex
        -p #    specify modulus in hex
        -n #    specify cycle length in hex
        -G #    specify generating point in hex
    <private key> can be specified in the last argument or interactively

Notes:

    Format for generating point G (X,Y) is:
    Uncompressed: 04 + X + Y
    Compressed: 02 + X if Y is even
                03 + X if Y is odd
    Where + is string concatenation and X and Y are in Hex


