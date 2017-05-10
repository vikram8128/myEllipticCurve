//
//  ec.cpp
//
//  Created by Vikram Singh on 17/02/2017.
//  Copyright © 2017 Vikram Singh.
//
//  This is my program to investigate elliptic curves.
//  It implements the elliptic curve group operation as described in
//  Understanding Cryptography by Christof Paar and Jan Pelzl.
//  It implements the double and add method for group multiplication.
//  By default, it uses the secp256k1 curve, but any elliptic curve can be
//  used by specifying the parameters in the command line arguments.
//  
//  To compile: g++ -O2 -std=c++11 -o ec ec.cpp -lgmpxx -lgmp
//
//  Notes:
//  Public key (X,Y) format is:
//   Uncompressed: 04 + X + Y
//   Compressed: 02 + X if Y is even
//		 03 + X if Y is odd
//   Where + is string concatenation and X and Y are in Hex
//  

#include <iostream>
#include <cstring>
#include <iomanip>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmpxx.h>
using namespace std;

char secp256k1_a_str[] = "0";
char secp256k1_b_str[] = "7";
char secp256k1_p_str[] = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
char secp256k1_n_str[] = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";
char secp256k1_G_str[] = "0279BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";

bool verbose = false;

class Point {
public:
  mpz_class x,y;
  bool y_inf;
  Point() {
    y_inf = false;
  };
  Point(mpz_class x, mpz_class y) {
    this->x = x;
    this->y = y;
    y_inf = false;
  };

};

class Curve {
public:
  mpz_class a;
  mpz_class b;
  mpz_class p;
  mpz_class n;
  Point G;
 
  Curve() {

  }
};

mpz_class setHexMpz(const char* str) {
  char tempStr[strlen(str) + 3];
  tempStr[0] = '0';
  tempStr[1] = 'x';
  strcpy(tempStr+2, str);
  return mpz_class(tempStr);
}

void setGxGy(char* G_str, mpz_class& Gx, mpz_class& Gy, mpz_class a, mpz_class b, mpz_class p) {
  if ((G_str[0] != '0') || 
      ((G_str[1] != '2') && (G_str[1] != '3') && (G_str[1] != '4'))) {
    cout << "Error: G is invalid" << endl;
    exit(5);
  }
  long len = strlen(G_str);
  char s[len+1];
  if (G_str[1] == '4') {
    strncpy(s, G_str + 2, (len-2)/2);
    s[(len-2)/2] = '\0';
    Gx = setHexMpz(s);
    strncpy(s, G_str + 2 + (len-2)/2, (len-2)/2); 
    s[(len-2)/2] = '\0';
    Gy = setHexMpz(s);
    return;
  }
  
  // check if p = 3 (mod 4)
  mpz_t res;
  mpz_class temp("4");
  mpz_init(res);
  mpz_mod(res, p.get_mpz_t(), temp.get_mpz_t());
  if (cmp(mpz_class(res),3)) {
    cout << "Error: p must be 3 (mod 4)" << endl;
    exit(5);
  }

  Gx = setHexMpz(G_str + 2);

  bool oddy = G_str[1] == '3';
  
  mpz_class r;
  r = Gx*Gx*Gx + a*Gx + b;
  mpz_mod(res, r.get_mpz_t(), p.get_mpz_t());
  r = mpz_class(res);

  // Take the square root of r
  temp = (p+1)/4;
  mpz_powm(res, r.get_mpz_t(), temp.get_mpz_t(), p.get_mpz_t());
  r = mpz_class(res);

  temp = mpz_class("2");
  mpz_mod(res, r.get_mpz_t(), temp.get_mpz_t());
  if ((cmp(mpz_class(res), 1) && oddy) || 
      (cmp(mpz_class(res), 0) && !oddy)) {
    r = p - r;
  }
  Gy = r;
}


void ecPlus(Curve& c, Point A, Point B, Point& res) {
  mpz_class s, temp;
  mpz_t r;
  mpz_init(r);
  if(A.y_inf) {
    res.x = B.x;
    res.y = B.y;
    res.y_inf = B.y_inf;
    return;
  }
  if(B.y_inf) {
    res.x = A.x;
    res.y = A.y;
    res.y_inf = A.y_inf;
    return;
  }
  
  if ((B.x == A.x) && (B.y == A.y)) {
    temp = A.y * 2;
    mpz_invert(r, temp.get_mpz_t(), c.p.get_mpz_t());
    temp = mpz_class(r) * (A.x*A.x*3 + c.a);
    mpz_mod(r, temp.get_mpz_t(), c.p.get_mpz_t());
    s = mpz_class(r); 
  } else { 
    temp = B.x - A.x;
    if (cmp(temp,0) == 0) {
      res.y_inf = true;
      return;
    }
    mpz_mod(r, temp.get_mpz_t(), c.p.get_mpz_t());
    temp = mpz_class(r);
    mpz_invert(r, temp.get_mpz_t(), c.p.get_mpz_t());
    temp = mpz_class(r) * (B.y - A.y);
    mpz_mod(r, temp.get_mpz_t(), c.p.get_mpz_t());
    s = mpz_class(r);
  }
  temp = s*s - A.x - B.x;
  mpz_mod(r, temp.get_mpz_t(), c.p.get_mpz_t());
  res.x = mpz_class(r);
  temp = s*(A.x - res.x) - A.y;
  mpz_mod(r, temp.get_mpz_t(), c.p.get_mpz_t());
  res.y = mpz_class(r);
  res.y_inf = false;
}

void ecMult(Curve& c, Point& pubK, mpz_class priK) {
  long bitCount = mpz_sizeinbase(priK.get_mpz_t(), 2);
  mpz_t r;
  mpz_init(r);
  mpz_class two(2);
  mpz_pow_ui(r, two.get_mpz_t(), bitCount-1);
  mpz_class bitI(r);

  pubK.x = c.G.x;
  pubK.y = c.G.y;
  bitI = bitI / 2;

  while (cmp(bitI, 0)) {
    ecPlus(c, pubK, pubK, pubK);
    mpz_and(r, bitI.get_mpz_t(), priK.get_mpz_t());
    if (cmp(mpz_class(r),0)) {
      ecPlus(c, pubK, c.G, pubK);
    }
    bitI = bitI / 2;
  }
}

void ec(char* a_str, char* b_str, char* p_str, char* n_str, char* G_str, const char* priK_str) {
  Curve c;
  c.a = setHexMpz(a_str);
  c.b = setHexMpz(b_str);
  c.p = setHexMpz(p_str);
  c.n = setHexMpz(n_str);
  mpz_class priK = setHexMpz(priK_str);
  
  setGxGy(G_str, c.G.x, c.G.y, c.a, c.b, c.p);

  cout << hex;
  if (verbose) {
    cout << "a = " << c.a << " (" << mpz_sizeinbase(c.a.get_mpz_t(), 2) << " bits)" << endl;
    cout << "b = " << c.b << " (" << mpz_sizeinbase(c.b.get_mpz_t(), 2) << " bits)" << endl;
    cout << "p = " << c.p << " (" << mpz_sizeinbase(c.p.get_mpz_t(), 2) << " bits)" << endl;
    //cout << "n = " << c.n << " (" << mpz_sizeinbase(c.n.get_mpz_t(), 2) << " bits)" << endl;
    cout << "Gx = " << c.G.x << " (" << mpz_sizeinbase(c.G.x.get_mpz_t(), 2) << " bits)" << endl;
    cout << "Gy = " << c.G.y << " (" << mpz_sizeinbase(c.G.y.get_mpz_t(), 2) << " bits)" << endl;
    cout << "-------" << endl;
    verbose = false;
  }

  cout << "private key: " << endl << "  " << priK << " (" << mpz_sizeinbase(priK.get_mpz_t(), 2) << " bits)" << endl;

  Point pubK;
  ecMult(c, pubK, priK);

  cout << "public key: " << endl;
  if (pubK.y_inf) {
    cout << " *Point at Infinity*" << endl;
  } else {
    cout << " x = " << setfill('0') << setw(64) << pubK.x << endl;
    cout << " y = " << setfill('0') << setw(64) << pubK.y << endl;
  }
}

void printhelp(char* argv[]) {
  cout << "usage: " << argv[0] << " [OPTIONS] <private key>" << endl << endl;
  cout << "The elliptic curve is y^2 = x^3 + ax + b \n\tand generating point G modulo p." << endl;
  cout << "\tAll numbers should be specified in hex, including the private key." << endl;
  cout << "The options are:" << endl;
  cout << "\t-v\tverbose print" << endl;
  cout << "\t-a #\tspecify a in hex" << endl;
  cout << "\t-b #\tspecify b in hex" << endl;
  cout << "\t-p #\tspecify modulus in hex" << endl;
  cout << "\t-n #\tspecify cycle length in hex" << endl;
  cout << "\t-G #\tspecify generating point in hex" << endl;
  cout << "<private key> can be specified in the last argument or in stdin" << endl;
  cout << endl;
}

int main (int argc, char **argv) {
  char *a_str = secp256k1_a_str;
  char *b_str = secp256k1_b_str;
  char *p_str = secp256k1_p_str;
  char *n_str = secp256k1_n_str;
  char *G_str = secp256k1_G_str;

  char c;

  while ((c = getopt (argc, argv, "va:b:p:n:G:")) != -1) {
    switch (c) {
      case 'v':
        verbose = true;
        break;
      case 'a':
        a_str = optarg;
        break;
      case 'b':
        b_str = optarg;
        break;
      case 'p':
        p_str = optarg;
        break;
      case 'n':
        n_str = optarg;
        break;
      case 'G':
        G_str = optarg;
        break;
      case '?':
        if ((optopt == 'a') || (optopt == 'b') || (optopt == 'p') || 
            (optopt == 'n') || (optopt == 'G') || (optopt == 'v')) 
          printf("Option -%c requires an argument.\n", optopt);
        else 
          printf("Unknown option `-%c'.\n", optopt);
      default:
        printhelp(argv);
        exit(6);
    }
  }

  if (optind == argc) {
    string priK_str;
    printf("Enter private key in hex or ^D to quit:\n");
    while(getline(cin, priK_str)) {
      ec(a_str, b_str, p_str, n_str, G_str, priK_str.c_str());
      printf("\nEnter private key in hex or ^D to quit:\n");
    }
    printf("done\n");
    exit(0);
  } else if (optind + 1 != argc) {
    printf("Error: Too many arguments:\n");
    for (int i = optind; i < argc; i++)
      printf ("  unknown argument %s\n", argv[i]);
    printhelp(argv); 
    exit(8);
  }

  ec(a_str, b_str, p_str, n_str, G_str, argv[optind]);
    
  return 0;
}


