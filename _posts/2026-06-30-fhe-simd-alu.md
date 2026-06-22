---
layout: post
title: >
  FHE for SIMD Arithmetic Logic Units with Amortized O(1) Bootstrapping per Ciphertext
date: 2026-06-22 00:00:00-0000
description: >
  TL;DR: We propose a new CKKS-compatible encoding framework that uses multiple complex slots to represent one number in a trianglular form, enabling both arithmetic and Boolean operations for a vector of, for example, 64-bit integers with amortized O(1) bootstrapping when O(n) ciphertexts are batched. In particular, for integer-arithmetic-only workloads, each refresh requires only two bootstrapping operations, regardless of the bit-width. The prototype is available at <https://github.com/tsinghua-ideal/fhe-simd-alu>. 
author: Hongren Zheng
tags: 
categories: 
related_posts: false
toc:
  sidebar: right
giscus_comments: true
---


- Written by [Hongren Zheng](https://hongrenzhe.ng) (Tsinghua University)
- Based on [https://ia.cr/2026/233](https://eprint.iacr.org/2026/233)

## Arithmetic and Boolean Operations 

In software, we usually need both arithmetic and Boolean operations. One example is computing a value and then using it to condition execution. The computation itself mostly involves arithmetic operations, while the logical step often involves Boolean operations.

```c
long func(long a, long b, long c) {
    long res = 42 * a + 100 * b - c; // Arithmetic Add/Mult
    if (res < 1024) { // Comparison
        // ...
    }
}
```

In this article, we focus on computer integers, also known as "machine words". They are among the most basic types in programming.

For arithmetic operations, considering only addition, subtraction and multiplication (side note: division can be understood as a mixture of arithmetic and Boolean operations), $n$-bit machine words can be understood as elements in the ring $\\mathbb{Z}_{2^{n}}$. We refer to this representation as the "whole" paradigm.

On the other hand, Boolean operations require the *bit representation* of machine words. For example, the "less than" comparison operation for integers $a$ and $b$ requires the most significant bit (MSB) of $(a - b)$. We refer to this representation as the "bit-vector" paradigm.

In plaintext computation, it is relatively cheap to convert between the two paradigms. However, in homomorphic encryption, the situation is different. For the "whole" paradigm, we have the following limitations:

  + BGV/BFV: They typically support $\mathbb{Z}\_{p}$ arithmetic with prime $p$ (directly supporting $\mathbb{Z}\_{2^n}$ leads to limited SIMD capability). Applications can only use $\mathbb{Z}\_{p}$ to simulate $\mathbb{Z}\_{2^n}$ semantics. However, unless we know the exact range of all operands, we can only use worst-case range analysis and $p$ is often large.
  + CKKS: It is possible to directly embed 64-bit integers in CKKS complex slots. However, after one multiplication, we obtain 128-bit integers in complex slots, and we need more-than-128-bit-precision bootstrapping to keep the lower bits (resulting message) and truncate the higher bits (overflows) to allow further computation.
  + "Whole" paradigm are hard to convert to the "bit vector" paradigm: they either need polynomial evaluation to extract all bits, or use a series of bootstrappings.

In this regard, the "bit-vector" paradigm is more natural for existing solutions, such as TFHE and RadixCKKS. However, these schemes then need to manage carries and overflows *iteratively*. A typical TFHE  multiplication circuit uses $O(n^2)$ programmable bootstrapping operations. RadixCKKS requires $O(\log (n/4))$ CKKS bootstrapping operations after each multiplication.

<div align="center" markdown="1">

|Scheme| Arithmetic | Boolean | 
|:---:|:---:|:---:|
|BGV/BFV with $$\mathbb{Z}_{2^n}$$ | No SIMD | Hard |
|BGV/BFV with $$\mathbb{Z}_{p}$$ | Good but range issue | Hard |
|CKKS | Precision issue | Hard |
| TFHE | Expensive carry logic | Natural |
| RadixCKKS | $O(\log (n/4))$ BTS | $O(1)$ BTS |
| This work | "Leveled" with $O(1)$ BTS | Amortized $O(1)$ BTS |

</div>

<div class="caption">
    Table 1. Comparison of schemes
</div>

## Triangle Encoding

We depart from existing paradigms and introduce a middle ground called "triangle encoding". It addresses the *iterative bootstrapping* issue by representing messages, carries, and overflows in a form that allows efficient extraction and cleanup.

Recall that REFHE (Eurocrypt'26) used a ring isomorphism from the "whole" paradigm to a polynomial quotient ring:

$$
\mathbb{Z}_{2^n} \cong \mathbb{Z}[X]/\langle X^n-X+2, X-2 \rangle
$$

The ring $\mathbb{Z}[X]$ contains all polynomials $f(X)=\sum_{k=0}^{d}f_k X^k$ with integer coefficients of any polynomial degree $d$. The quotient structure means the polynomial addition and multiplication are carried out modulo both $$X^n-X+2$$ and $$X-2$$. Concretely, this can be understood as replacing every occurrence of $X^n$ with $X-2$ and then replacing every occurrence of $X-2$ with 0.

More specifically, each element $$m\in \mathbb{Z}_{2^n}$$ has a corresponding *flattened* polynomial:

$$
[m]_t = b_0 + b_1 X + \cdots + b_{n-1} X^{n-1}
$$

The $i$-th coefficient is the $i$-th bit of $m$, since $$m=\sum_{k=0}^{n-1}b_k 2^k=\sum_{k=0}^{n-1}b_k X^k$$.
Here we use the notation $t=X-2$. Every polynomial addition and multiplication of flattened polynomials modulo both $$X^n-X+2$$ and $$X-2$$ corresponds to addition and multiplication in $$\mathbb{Z}_{2^n}$$, respectively.

However, maintaining both quotient relations homomorphically is difficult. Our strategy is to simulate the two quotient relations step by step. We first work in the polynomial ring $\mathbb{Z}[X]/\langle X^n-X+2 \rangle$ which enforces the first relation. To simulate the second modulo structure, we then separate any polynomial $f(X)$ in the polynomial ring $\mathbb{Z}[X]/\langle X^n-X+2 \rangle$ into two parts:

$$
f(X) = [m]_t + (X-2) \cdot I(X)
$$

Here, $[m]_t$ is the remainder and $I(X)$ is the quotient. Equivalently, $(X-2) I(X)$ lives in the $\langle X-2\rangle$ polynomial ideal in $$\mathbb{Z}[X]/\langle X^n-X+2 \rangle$$, and quotienting by this ideal gives us the flattened polynomial.

Our key observation is that we can turn the polynomial ideal into an integer ideal by multiplying by $t^{-1}$ in an appropriate ring. We do this in the the ring $\mathbb{R}[X]/\langle X^n-X+2\rangle$, where we have the inverse of $t$. A direct computation in this ring gives:

$$
I(X) + \frac{[m]_t}{t}= \sum_{k=0}^{n-1}I_kX^k + \left(\textcolor{gray}{\frac{m}{2^n}} - \sum_{k=0}^{n-1}\frac{1}{2^{k+1}}{\underset{0\leq i\leq k}{[m]}}X^{k}\right)
$$

where $\underset{0\leq i\leq k}{[m]} = \sum_{i=0}^{k}b_i2^i$ is the lower $k+1$ bits of $m$. After scaling by $2^{k+1}$, these prefix sums form a triangular pattern.

We now take homomorphic encryption concerns into account. In RLWE, homomorphic encryption will introduce noise. To actually simulate the $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$ arithmetic with noise, we follow the idea of the CKKS scaling factor. The triangle encoding has the following form, where $e$ denotes the noise:

$$
\Delta \left( I + \frac{[m]_t}{t}\right) + e
$$

Figure 1 illustrates the encoding. The extra term $m/2^n$ in the constant coefficient is not shown. Here, each column represents one coefficient. The first column is the constant coefficient, the second column is the coefficient for $X^1$, and so on. The column is displayed from low bits to high bits.

<div class="row mt-3">
    <div class="col-sm-12 mt-3 mt-md-0 mx-auto d-block">
        {% include figure.liquid loading="eager" path="assets/img/blog/2606_Hongren/triangle.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 1. Triangle Encoding
</div>

This representation has the following useful properties:
  * (Encoding Correctness) The noise $e$ can contaminate part of the triangle, as long as the "upper side" (in purple) of the triangle is not contaminated.
  * (Arithmetic Refreshing) If we can efficiently reset the overflow $I$ and the noise $e$, we allow further arithmetic operations.
  * (Arithmetic-To-Boolean Conversion) If we extract the "upper side" of the triangle, we get a Boolean representation of the number.

## SIMD Homomorphism Chain

Before discussing these management operations, we first explain how to obtain SIMD, since it is crucial for application efficiency. We also want to use the standard RLWE ring $$\mathcal{R}_Q=\mathbb{Z}[X]/\langle X^N+1, Q\rangle$$ (e.g., $$N=2^{16}$$) that is widely studied and accelerated.

The idea is to use a *chain* of homomorphisms, as shown in Figure 2. We first connect the RLWE ring to complex slots, which is standard CKKS. Then, we group the complex slots, with each group containing $n/2$ complex slots. The crucial observation is that there is an isomorphism from $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$ to $$\mathbb{C}^{n/2}$$, because $X^n-X+2$ for power-of-two $n$ (e.g., $n=64$) has $n/2$ pairs of conjugate roots with no repeated roots. We then use the triangle encoding to embed $$\mathbb{Z}_{2^n}$$ into the polynomial ring structure. Thus, we obtain $N/n$ $$\mathbb{Z}_{2^n}$$ slots from the usual RLWE ring. In a typical parameterization where $N=2^{16}$ and $n=64$, we obtain 1024 slots.
<div class="row mt-3">
    <div class="col-sm-12 mt-3 mt-md-0 mx-auto d-block">
        {% include figure.liquid loading="eager" path="assets/img/blog/2606_Hongren/simd.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 2. SIMD Homomorphism Chain
</div>

To clarify, the homomorphisms in this chain are evaluated in plaintext during encoding/decoding. Encryption and homomorphic computation always take place in the ring $$\mathbb{Z}[X]/\langle X^N+1, Q\rangle$$.

## O(1) Bootstrapping for Arithmetic Refreshing

The goal of arithmetic refreshing is to reset the overflow $I$ and noise $e$ so that arithmetic computations can continue. We first focus on the overflow $I$, which has been made into a regular structure by our triangle encoding.

On the other hand, the homomorphism chain makes the encoded polynomial appear "irregular" in $$\mathbb{Z}[X]/\langle X^N+1, Q\rangle$$. In CKKS, there are linear transformations that move data between coefficients and slots. Here, we use linear transformations to move coefficients from $$(\mathbb{R}[X]/\langle X^n-X+2\rangle)^{N/n}$$ to $$\mathcal{R}_Q=\mathbb{Z}[X]/\langle X^N+1, Q\rangle$$.

After applying linear transformations, the ciphertext has the $$\Delta (I + [m]_t/t) + e$$ coefficient structure in the RLWE ring $$\mathcal{R}$$. We can then reset the overflow $I$ using the "Truncate and ModRaise" operation. The "Truncate" operation is inspired by "Modular Reduction in CKKS" which multiplies by $$Q_{\ell}/\Delta$$ to remove the higher bits, and then modulus switching back to $$q_0$$. The ModRaise operation then introduces another $I'$ while raising the modulus back to $$Q_{L}$$. By using a sparse key encapsulation technique, the norm of the new $I'$ can be very small.

<div class="row mt-3">
    <div class="col-sm-12 mt-3 mt-md-0 mx-auto d-block">
        {% include figure.liquid loading="eager" path="assets/img/blog/2606_Hongren/truncate.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 3. The "Truncate and ModRaise" Operation
</div>

After the "Truncate and ModRaise" operation, we again use linear transformations to move $$\Delta(I' + [m]_t/t) + e$$ back to the original representation. The "linear transformations + Truncate and modRaise + linear transformations" structure has a cost similar to one CKKS bootstrapping, dominated by the linear transformations.

To manage the noise $e$, we follow a similar idea to META-BTS. We first multiply by $t$ in the original ring to obtain $$\Delta (tI + [m]_t) + te$$, so that $$\Delta (tI + [m]_t)$$ can be truncated away using "Truncate and ModRaise" as $(tI + [m]_t)$ has integer coefficients. We then bootstrap $te$, recover $e$, and subtract it from the input ciphertext.

Therefore, we only need two CKKS bootstrapping operations for arithmetic-to-arithmetic refreshing. We use $O(1)$ notation as the number does not change with the bit-width $n$, unlike previous works. We further remark that refreshing is needed only when $I$ or $e$ becomes large enough, so several multiplications can be performed before one refreshing ("leveled"), and both bootstrapping operations can be executed independently.

## Amortized O(1) Bootstrapping for Arithmetic-To-Boolean Conversion


To support Boolean operations, we need a representation where the bits are placed in the complex slots

$$
\Delta (b_0 , b_1, \cdots, b_{n-1}) \in \mathbb{C}^{n}
$$

In this representation, bitwise AND/OR can be simulated as $(ab)$ and $(a + b - ab)$ over complex slots. Bit shifting can be simulated using CKKS rotation.

On the other hand, this representation resembles the "upper side" (purple line) in Figure 1. Therefore, the central task of Arithmetic-to-Boolean conversion is to extract the bits of the triangle efficiently.
We omit some details here and focus on the central idea.

The observation is that, because the representation is a triangle, we can extract $b_0$ using one bootstrapping (we need bootstrapping due to the overflow $I$) and subtract the $b_0$ from other coefficients (in concrete execution, we remove some of them and treat other lower bits as noise). Then, $b_1$ becomes exposed to us and we can continue the extraction. Therefore, after $O(n)$ bootstrapping operations, we obtain all the bits.

<div class="row mt-3">
    <div class="col-sm-12 mt-3 mt-md-0 mx-auto d-block">
        {% include figure.liquid loading="eager" path="assets/img/blog/2606_Hongren/bit-remove.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 4. Extract $b_0$ then remove $b_0$ in other coefficients
</div>

Then we observe that, in the bootstrapping for each bit $b_i$, we actually use one real slot out of $n/2$ complex slots occupied by the triangle. This allows us to combine the $b_i$'s of $O(n)$ ciphertexts into one input ciphertext for one bootstrapping, and we can use $O(n)$ bootstrapping to process $O(n)$ ciphertexts, resulting in amortized O(1) bootstrapping per ciphertext.

<div class="row mt-3">
    <div class="col-sm-12 mt-3 mt-md-0 mx-auto d-block">
        {% include figure.liquid loading="eager" path="assets/img/blog/2606_Hongren/batch.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 5. Batch the bootstrapping of $b_0$ in $O(n)$ ciphertexts
</div>

In our concrete execution, the bit extraction is carried out in blocks (4 bits per bootstrapping). We need $n/4$ bootstrapping to process $n/4$ ciphertexts, therefore amortized to 1 bootstrapping per ciphertext.

## Summary

In summary, we introduce triangle encoding, a new encoding for machine words. In arithmetic mode, each refresh uses two bootstrapping-like operations, independent of the word size. When we want to mix in Boolean operations, we can achieve amortized $O(1)$ bootstrapping when there are a sufficient number of ciphertexts.

Our scheme is most suitable for applications with heavy arithmetic and lightweight Boolean operations. The recommended word sizes for applications are $n=32$ and $n=64$, as these are common word sizes in modern computer. We also demonstrated the use of $n=256$ for high-precision arithmetic in the paper. In principle, we can support arbitrary bit-width $n$, but then the linear transformation cost becomes dominant, unlike in RadixCKKS where a DFT-like structure can be used.

This post focuses on the amortized conversion strategy presented in the paper. More direct arithmetic-to-Boolean conversion methods, without relying on the same amortization regime, will be discussed separately.