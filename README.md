
# BP+OSD: A decoder for quantum LDPC codes 
A Python library implementing belief propagation with ordered statistics post-processing for decoding sparse quantum LDPC codes as described in [arXiv:2005.07016](https://arxiv.org/abs/2005.07016). Note, this library has recently been completly rewritten using Python and Cython. The bulk of the code now resides in the [LDPC](https://github.com/quantumgizmos/ldpc) repository. The original C++ version can be found in the `cpp_version` branch of this repository.

## Installation from PyPi (recommended method)

Installation from [PyPi](https://pypi.org/project/bposd/) requires Python>=3.6.
To install via pip, run:

```
pip install -U bposd
```

## Documentation
This package buids upon the [LDPC](https://github.com/quantumgizmos/ldpc) python package. The documentation for LDPC can be found [here](https://roffe.eu/software/ldpc/index.html).

## Attribution
If you use this software in your research, please cite the following research paper:

```
@article{roffe_decoding_2020,
   title={Decoding across the quantum low-density parity-check code landscape},
   volume={2},
   ISSN={2643-1564},
   url={http://dx.doi.org/10.1103/PhysRevResearch.2.043423},
   DOI={10.1103/physrevresearch.2.043423},
   number={4},
   journal={Physical Review Research},
   publisher={American Physical Society (APS)},
   author={Roffe, Joschka and White, David R. and Burton, Simon and Campbell, Earl},
   year={2020},
   month={Dec}
}
```

Please also cite the LDPC software package:

```
@software{Roffe_LDPC_Python_tools_2022,
author = {Roffe, Joschka},
title = {{LDPC: Python tools for low density parity check codes}},
url = {https://pypi.org/project/ldpc/},
year = {2022}
}
```

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
from ldpc.codes import rep_code
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
from ldpc import bposd_decoder

bpd=bposd_decoder(
    surface_code.hz,#the parity check matrix
    error_rate=0.05,
    channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
    max_iter=surface_code.N, #the maximum number of iterations for BP)
    bp_method="ms",
    ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
    osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
    osd_order=7 #the osd search depth
    )
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
print(bpd.osdw_decoding)
#Decoding is successful if the residual error commutes with the logical operators
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
    Logical Error: No
    

