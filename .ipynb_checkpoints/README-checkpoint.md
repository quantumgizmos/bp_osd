# BP+OSD: A decoder for sparse quantum codes
A python library implementing belief propagation with ordered statistics post-processing for decoding sparse quantum codes as described in [arXiv:2005.07016](https://arxiv.org/abs/2005.07016). 

## Installation from PyPi

Installtion from PyPi requires Python>=3.6.
To install via pip, run:

```
pip install bposd
```

## Installation (from source)

Installation from sources requires Python>=3.6 and a local C compiler (eg. 'gcc' in Linux or 'clang' in Windows). Once these requirements have been met, navigate to the repository root and install using pip:

```
git clone https://github.com/quantumgizmos/bp_osd.git
cd bp_osd
pip install -e .
```

## Demo Scripts

To check whether the package is working run one of the demo scripts in the examples folder.

## Software
This library makes use of the following software:
- [Sofware for low density parity check codes](https://github.com/radfordneal/LDPC-codes), Radford M. Neal

## license

MIT License

Copyright (c) 2020 Joschka Roffe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
