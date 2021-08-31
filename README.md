# BP+OSD: A decoder for quantum LDPC codes 
A C library implementing belief propagation with ordered statistics post-processing for decoding sparse quantum LDPC codes as described in [arXiv:2005.07016](https://arxiv.org/abs/2005.07016).

## Installation from PyPi (recommended method)

Installtion from [PyPi](https://pypi.org/project/bposd/) requires Python>=3.6.
To install via pip, run:

```
pip install bposd
```

## Installation (from source)

Installation from source requires Python>=3.6 and a local C compiler (eg. 'gcc' in Linux or 'clang' in Windows). The LDPC package can then be installed by running:

```
git clone https://github.com/quantumgizmos/bposd.git
cd bposd
pip install -e bposd
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


Not all combinations of the `hx` and `hz` matrices will produce a valid CSS code. Use the `bposd.css.css_code.test` function to check whether the code is valid. For example, we can easily check that the Steane code passes all the CSS code tests:


```python
steane_code.test()
```

    <Unnamed CSS code>, (3,4)-[[7,1,nan]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (3,4)-[[7,1,nan]]





    True



As an example of a code that isn't valid, consider the case when `hx` and `hz` are repetition codes:


```python
from ldpc.codes import rep_code

hx=hz=rep_code(7)
qcode=css_code(hx,hz)
qcode.test()
```

    <Unnamed CSS code>, (2,2)-[[7,-5,nan]]
     -Block dimensions incorrect
     -PCMs commute hz@hx.T==0: Fail
     -PCMs commute hx@hz.T==0: Fail
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anitcommute: Fail





    False



## Hypergraph product codes

The hypergraph product can be used to construct a valid CSS code from any pair of classical seed codes. To use the the hypergraph product, call the `bposd.hgp.hgp` function. Below is an example of how the distance-3 surface code can be constructed by taking the hypergraph product of two distance-3 repetition codes.


```python
from bposd.hgp import hgp
h=rep_code(3)
surface_code=hgp(h1=h,h2=h,compute_distance=True) #nb. set compute_distance=False for larger codes
surface_code.test()
```

    <Unnamed CSS code>, (2,4)-[[13,1,3]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (2,4)-[[13,1,3]]





    True



## BP+OSD Decoding

BP+OSD decoding is useful for codes that do not perform well under standard-BP. To use the BP+OSD decoder, we first call the `bposd.bposd_decoder` class:


```python
import numpy as np
from bposd import bposd_decoder

bpd=bposd_decoder(
    surface_code.hz,#the parity check matrix
    error_rate=0.05,
    max_iter=surface_code.N, #the maximum number of iterations for BP)
    bp_method="ps",
    osd_method="osd_cs", #the OSD method
    osd_order=-1)
```

We can then decode by passing a syndrome to the `bposd.bposd_decoder.decode` method:


```python
error=np.zeros(surface_code.N).astype(int)
error[[5,12]]=1
syndrome=surface_code.hz@error %2
bpd.decode(syndrome)

print("Error")
print(error)
print("BP+OSD Decoding")
print(bpd.bp_decoding)
print(bpd.osd0_decoding)
print(bpd.osdw_decoding)
residual_error=(bpd.osdw_decoding+error) %2
a=(surface_code.lz@residual_error%2).any()
if a: a="Yes"
else: a="No"
print(f"Logical Error: {a}\n")

```

    Error
    [0 0 0 0 0 1 0 0 0 0 0 0 1]
    BP+OSD Decoding
    [0 0 0 0 0 0 0 0 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 1 0 0 0 0]
    Logical Error: No
    



```python



def random_error(N,error_rate):
    error=np.zeros(N).astype(int)

    for i in range(N):
        if np.random.random()<error_rate: error[i]=1
        else: error[i]=0

    return error

error_rate=0.10

bpd=bposd_decoder(
    surface_code.hz,#the parity check matrix
    error_rate=0.05,
    max_iter=3, #the maximum number of iterations for BP)
    bp_method="ps",
    osd_method="osd_e", #the OSD method
    osd_order=7)

for i in range(3):
    error=random_error(N=surface_code.N,error_rate=error_rate)
    syndrome=surface_code.hz@error%2
    bpd.decode(syndrome)
    print("Error")
    print(error)
    print("BP+OSD Decoding")
    print(bpd.osdw_decoding)
    print(bpd.converge)
    a=(surface_code.lz@((bpd.osdw_decoding+error) %2)%2).any()
    if a: a="Yes"
    else: a="No"
    print(f"Logical Error: {a}\n")
```

    Error
    [0 0 0 0 0 0 0 0 0 0 0 0 0]
    BP+OSD Decoding
    [0 0 0 0 0 0 0 0 0 0 0 0 0]
    1
    Logical Error: No
    
    Error
    [0 0 0 0 0 0 0 0 0 0 0 0 0]
    BP+OSD Decoding
    [0 0 0 0 0 0 0 0 0 0 0 0 0]
    1
    Logical Error: No
    
    Error
    [0 0 0 0 0 0 0 0 0 0 0 0 0]
    BP+OSD Decoding
    [0 0 0 0 0 0 0 0 0 0 0 0 0]
    1
    Logical Error: No
    



```python
surface_code.hz
```




    array([[1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
           [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0],
           [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1],
           [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1]])




```python

```
