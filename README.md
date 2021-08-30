# BP+OSD: A decoder for quantum LDPC codes 
A C library implementing belief propagation with ordered statistics post-processing for decoding sparse quantum LDPC codes as described in [arXiv:2005.07016](https://arxiv.org/abs/2005.07016).

## Installation from PyPi (recommended method)

Installtion from [PyPi](https://pypi.org/project/ldpc/) requires Python>=3.6.
To install via pip, run:

```
pip install ldpc
```

## Installation (from source)

Installation from source requires Python>=3.6 and a local C compiler (eg. 'gcc' in Linux or 'clang' in Windows). The LDPC package can then be installed by running:

```
git clone https://github.com/quantumgizmos/ldpc.git
cd ldpc
pip install -e ldpc
```

## Dependencies
This package makes use of the `mod2sparse` data structure from Radford Neal's [Software for Low Density Parity Check Codes](https://www.cs.toronto.edu/~radford/ftp/LDPC-2012-02-11/index.html) C package.



# Basic usage

## Constructing CSS codes

The `bposd.css.css_code` class can be used to create a CSS code from two classical codes. As an example, we can create a [[7,4,3]] Steane code from the classical Hamming code


```python
from ldpc.codes import hamming_code
from bposd.css import css_code
h=hamming_code(3) #Hamming code parity check matrix
steane_code=css_code(hx=h,hz=h) #create Steane code where both hx and hz are Hamming codes
print("Hx")
print(steane_code.hx)
print("Hz")
print(steane_code.hz)
```

    Hx
    [[0 0 0 1 1 1 1]
     [0 1 1 0 0 1 1]
     [1 0 1 0 1 0 1]]
    Hz
    [[0 0 0 1 1 1 1]
     [0 1 1 0 0 1 1]
     [1 0 1 0 1 0 1]]


The `bposd.css.css_code` class automatically computes the logical operators of the code.


```python
print("Lx Logical")
print(steane_code.lx)
print("Lz Logical")
print(steane_code.lz)
```

    Lx Logical
    [[1 1 1 0 0 0 0]]
    Lz Logical
    [[1 1 1 0 0 0 0]]


## BP+OSD Decoding

BP+OSD decoding is useful for codes that do not perform well under standard-BP. To use the BP+OSD decoder, we first call the `bposd.bposd_decoder` class:


```python
import numpy as np
from bposd import bposd_decoder

bpd=bposd_decoder(
    steane_code.hx,#the parity check matrix
    error_rate=0.1,
    max_iter=steane_code.N #the maximum number of iterations for BP)
)

for i in range(steane_code.N):
    error=np.zeros(steane_code.N).astype(int)
    error[i]=1
    syndrome=steane_code.hz@error%2
    bpd.decode(syndrome)
    print(bpd.osdw_decoding)
```

    [0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0]



```python

```


```python

```
