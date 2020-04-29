#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mod2sparse.h"

int syndrome
( mod2sparse *H,	/* Parity check matrix */
  char *er,		/* Guess for codeword */
  char *pchk		/* Place to store parity checks */
)
{
  mod2sparse_mulvec (H, er, pchk);

  return 0;

}