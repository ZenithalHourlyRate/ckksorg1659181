---
layout: post
title: >
  FHE for SIMD Arithmetic Logic Units with Amortized O(1) Bootstrapping per Ciphertext
date: 2026-06-22 00:00:00-0000
description: >
  TL;DR: We propose a new CKKS-compatible encoding framework that supports both arithmetic and Boolean operations for a vector of, for example, 64-bit integers. The key idea is using multiple complex slots to represent one integer, with a special ring isomorphism to maintain the desired integer arithmetic. For arithmetic-only workloads, each refreshing requires only two bootstrapping operations for one ciphertext. For Boolean operations, the arithmetic-to-Boolean conversion can batch $O(n)$ ciphertexts, resulting in amortized $O(1)$ bootstrapping per ciphertext. The prototype is available at <https://github.com/tsinghua-ideal/fhe-simd-alu>. 
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

In this regard, the "bit-vector" paradigm is more natural for existing solutions, such as TFHE and RadixCKKS (Eurocrypt'26) [6]. However, these schemes then need to manage carries and overflows *iteratively*. A typical TFHE multiplication circuit uses $O(n^2)$ programmable bootstrapping operations. RadixCKKS requires $O(\log (n/4))$ CKKS bootstrapping operations after each multiplication.

<div align="center" markdown="1">

|Scheme| Arithmetic | Boolean | 
|:---:|:---:|:---:|
|BGV/BFV with $$\mathbb{Z}_{2^n}$$ | No SIMD | Hard |
|BGV/BFV with $$\mathbb{Z}_{p}$$ | Good but range issue | Hard |
|CKKS | Precision issue | Hard |
| TFHE | Expensive carry logic | Natural |
| RadixCKKS [6] | $O(\log (n/4))$ BTS | $O(1)$ BTS |
| REFHE [5] | "Leveled" with $O(n)$ PBS | $O(n)$ PBS |
| This work | "Leveled" with $O(1)$ BTS | Amortized $O(1)$ BTS |

</div>

<div class="caption">
    Table 1. Comparison of schemes
</div>

## Triangle Encoding

We depart from existing paradigms and introduce a middle ground called "triangle encoding". It addresses the *iterative bootstrapping* issue by representing messages, carries, and overflows in a form that allows efficient extraction and cleanup.

This work follows the line of "polynomial plaintext modulus" technique. This technique dates all the way back to 2000 from Hoffstein and Silverman in the NTRU literature [1]. In the FHE literature, this idea appears as early as CLPX (CT-RSA'18) [2], and recent works like GBFV (Eurocrypt'25, S&P'25) [3,4] and REFHE (Eurocrypt'26) [5] quickly gain gravity in the FHE community. Their bootstrapping methods are constantly evolving: from leveled HE, to FHE with BFV-style bootstrapping, to FHE with TFHE-style bootstrapping. Technically, we tweak the flattened polynomial further into a triangular format, making it friendly to CKKS-style bootstrapping and obtaining (amortized) O(1) bootstrapping.

Let's directly jump into REFHE (Eurocrypt'26). It used a ring isomorphism from the "whole" paradigm to a polynomial quotient ring:

$$
\mathbb{Z}_{2^n} \cong \mathbb{Z}[X]/\langle X^n-X+2, X-2 \rangle
$$

The ring $\mathbb{Z}[X]$ contains all polynomials $f(X)=\sum_{k=0}^{d}f_k X^k$ with integer coefficients of any polynomial degree $d$. The quotient structure means the polynomial addition and multiplication are carried out modulo both $$X^n-X+2$$ and $$X-2$$. Concretely, this can be understood as replacing every occurrence of $X^n$ with $X-2$ and then replacing every occurrence of $X-2$ with 0.

To instantiate the ring isomorphism, there is a corresponding  "flatten operation" that expands a scalar message into a polynomial while preserving its arithmetic properties. It can be understood as modulo by a polynomial, therefore it uses the notation $$[\cdot]_t$$. Here we use $$t=X-2$$. Formally, for each element $$m\in \mathbb{Z}_{2^n}$$, there is a corresponding *flattened* polynomial in the ring $$\mathbb{Z}[X]/\langle X^n-X+2 \rangle$$:

$$
[m]_t = b_0 + b_1 X + \cdots + b_{n-1} X^{n-1}
$$

The $k$-th coefficient $$b_k\in \{0, 1\}$$ is the $k$-th bit of $m$. Another way to understand it is $$m=\sum_{k=0}^{n-1}b_k 2^k=\sum_{k=0}^{n-1}b_k X^k$$.
Note that as we modulo by $X-2$, which means we can replace every occurence of $X$ with 2, the modulo-by-$(X^n-2+2)$ structure can be understood as a modulo-by-$2^n$ operation. Therefore, every polynomial addition and multiplication of flattened polynomials modulo both $$X^n-X+2$$ and $$X-2$$ correspond to addition and multiplication in $$\mathbb{Z}_{2^n}$$, respectively.

However, maintaining the flattened structure homomorphically is difficult, as it involves two quotient relations. The usual strategy in existing works is to use the RLWE ring to induce the first relation, while using "plaintext modulus" to induce the second relation. However, this makes the arithmetic tightly coupled with the RLWE ring, like those in CLPX and GBFV. REFHE chose to change the RLWE ring, but it then became inconvenient for implementation.

Our strategy is to simulate the two quotient relations step by step. Assume now that we can obtain $\mathbb{R}[X]/\langle X^n-X+2 \rangle$ arithmetic which enforces the first relation. To simulate the second modulo structure, we then separate any polynomial $f(X)$ in the polynomial ring $\mathbb{Z}[X]/\langle X^n-X+2 \rangle$ into two parts:

$$
f(X) = [m]_t + (X-2) \cdot I(X)
$$

Here, $[m]_t$ is the remainder and $I(X)$ is the quotient, a polynomial with integer coefficients. Equivalently, $(X-2) I(X)$ lives in the $\langle X-2\rangle$ polynomial ideal in $$\mathbb{Z}[X]/\langle X^n-X+2 \rangle$$, and quotienting by this ideal gives us the flattened polynomial.

However, the polynomial ideal is hard to deal with using homomorphic operations. REFHE needs to manage it iteratively. Our key observation is that we can turn the polynomial ideal into an integer ideal by multiplying by $t^{-1}$, so that the $I(X)$ can be managed in a coefficient-wise manner. We have the inverse of $t$ in the the ring $\mathbb{R}[X]/\langle X^n-X+2\rangle$. A direct computation in this ring gives:

$$
I(X) + \frac{[m]_t}{t}= \sum_{k=0}^{n-1}I_kX^k + \left(\textcolor{gray}{\frac{m}{2^n}} - \sum_{k=0}^{n-1}\frac{1}{2^{k+1}}{\underset{0\leq i\leq k}{[m]}}X^{k}\right)
$$

where $\underset{0\leq i\leq k}{[m]} = \sum_{i=0}^{k}b_i2^i$ is the lower $k+1$ bits of $m$. After scaling by $2^{k+1}$, these prefix sums form a triangular pattern.

We now take homomorphic encryption concerns into account. In RLWE, homomorphic encryption will introduce noise. To actually simulate the $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$ arithmetic with noise, we follow the idea of the CKKS scaling factor. The triangle encoding has the following form, where $e$ denotes the noise:

$$
\Delta \left( I + \frac{[m]_t}{t}\right) + e
$$

Figure 1 illustrates the encoding. Here, each column represents one polynomial coefficient. The first column is the constant coefficient, the second column is the coefficient for $X^1$, and so on. The column is displayed from low bits to high bits. The extra term $m/2^n$ in the constant coefficient is not shown.

<div class="row mt-3">
    <div class="col-sm-12 mt-3 mt-md-0 mx-auto d-block">
        {% include figure.liquid loading="eager" path="assets/img/blog/2606_Hongren/triangle.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 1. Triangle Encoding
</div>

This representation has the following useful properties:
  * (Exact Correctness) The noise $e$ can contaminate part of the triangle, as long as the "top row" (in purple) of the triangle is not contaminated. As long as $$\lVert te \rVert < \Delta /2$$, we can recover the integer *exactly*.
  * (Arithmetic Refreshing) If we can efficiently reset the noise $e$ (and the overflow $I$), we allow further arithmetic operations.
  * (Arithmetic-To-Boolean Conversion) If we extract the "top row" of the triangle, we get a Boolean representation of the number.

## SIMD Homomorphism Chain

Before discussing these management operations, we first explain how to obtain the $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$ arithmetic. We also explain how to obtain SIMD, which is crucial for application efficiency. 

Our starting point is the celebrated RLWE ring $$\mathcal{R}_Q=\mathbb{Z}[X]/\langle X^N+1, Q\rangle$$ (e.g., $$N=2^{16}$$) that is widely studied and accelerated. It can simulate complex slots using CKKS.

The key observation is that we can use some complex slots to simulate the $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$ arithmetic because of a ring isomorphism. Here we assume $n$ is a power-of-two, then $X^n-X+2$ has $n/2$ pairs of conjugate roots in the complex plane with no repeated roots.

$$
\mathbb{C}^{n/2} \cong \mathbb{R}[X]/\langle X^n-X+2\rangle
$$

The whole picture is Figure 2, where we obtain a *chain* of homomorphisms. We first connect the RLWE ring to complex slots. Then, we group the complex slots, with each group containing $n/2$ complex slots, and connect the complex slots to the polynomial ring $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$. We then use the triangle encoding to embed $$\mathbb{Z}_{2^n}$$ into the polynomial ring structure. In the end, we obtain $N/n$ number of $$\mathbb{Z}_{2^n}$$ slots from the usual RLWE ring. In a typical parameterization where $N=2^{16}$ and $n=64$, we obtain 1024 slots.
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

The goal of arithmetic refreshing is to reset the noise $e$ (and the overflow $I$) so that arithmetic computations can continue. 

The reason for managing the overflow $I$ is that, it contributes to the noise growth in multiplication. Suppose we have two triangles in the form of $$\Delta I + e$$ and $$\Delta I' + e'$$ (just two terms to simplify the discussion), after multiplication and rescaling, we can find that the new noise contains cross terms like $$Ie'$$ and $$I'e$$, and the new overflow contains $II'$ term. Therefore, smaller overflow $I$ will lead to smaller noise growth.

In the triangle encoding section, we have made the overflow $I$ into an expression that can be managed in a coefficient-wise manner in the ring $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$. On the other hand, the encoding process through the homomorphism chain makes the encoded polynomial appear "irregular" in the ring $$\mathbb{Z}[X]/\langle X^N+1, Q\rangle$$. Note that in the RLWE ring, there is a modulo-by-Q structure that *natively* modulo  every coeffcient, and we want to exploit that. 

The idea is to move the "regular" coefficients from one polynomial ring to another. In CKKS, there are linear transformations that move data between coefficients and slots. Here, we use linear transformations to move polynomial coefficients from $$(\mathbb{R}[X]/\langle X^n-X+2\rangle)^{N/n}$$ to $$\mathcal{R}_Q=\mathbb{Z}[X]/\langle X^N+1, Q\rangle$$. More technically, we first move from the coefficients from the first polynomial rings to $$C^{N/2}$$ using linear transformation that costs $$O(\sqrt{n})$$ key-switching, then use CKKS-style SlotsToCoeffs to move them into the RLWE ring.

After applying linear transformations, the ciphertext encrypts the $$\Delta (I + [m]_t/t) + e$$ coefficient structure in the RLWE ring $$\mathcal{R}$$. We can then reset the overflow $I$ using the "Truncate and ModRaise" operation. The "Truncate" operation is inspired by [7] which multiplies by $$Q_{\ell}/\Delta$$ to remove the higher bits, and then modulus switching back to $$q_0$$. The ModRaise operation then introduces another $I'$ while raising the modulus back to $$Q_{L}$$. By using a sparse key encapsulation technique, the norm of the coefficients of the new $I'$ can be very small.

<div class="row mt-3">
    <div class="col-sm-12 mt-3 mt-md-0 mx-auto d-block">
        {% include figure.liquid loading="eager" path="assets/img/blog/2606_Hongren/truncate.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 3. The "Truncate and ModRaise" Operation
</div>

After the "Truncate and ModRaise" operation, we again use linear transformations to move $$\Delta(I' + [m]_t/t) + e$$ back to the original representation. The "linear transformations + Truncate and modRaise + linear transformations" structure has a cost similar to one CKKS bootstrapping, dominated by the linear transformations.

To manage the noise $e$, we follow the idea of bootstrapping the noise, like in META-BTS [8]. We first multiply by $t$ in the original ring to obtain $$\Delta (tI + [m]_t) + te$$, so that $$\Delta (tI + [m]_t)$$ can be truncated away using "Truncate and ModRaise" as $(tI + [m]_t)$ has integer coefficients. We then bootstrap $te$, recover $e$, and subtract it from the input ciphertext.

Therefore, we only need two CKKS bootstrapping operations for arithmetic-to-arithmetic refreshing. We use $O(1)$ notation as the CKKS bootstrapping number does not change with the bit-width $n$, unlike previous works. We see CoeffsToSlots and SlotsToCoeffs as the major cost while seeing the $O(\sqrt{n})$ cost for linear transformation as minor because $n$ is moderate. We further remark that refreshing is needed only when $I$ or $e$ becomes large enough, so several multiplications can be performed before one refreshing ("leveled"), and both bootstrapping operations can be executed independently.

As for the concrete number of multiplications that can be performed, it depends on the noise and overflow size of the operands, and the scaling factor. In our parameter setting, we estimate 2 to 3 multiplications could be performed before a refreshing, depending on the whether the input operands are plaintext, fresh encryption or refreshed ciphertexts. 

## Amortized O(1) Bootstrapping for Arithmetic-To-Boolean Conversion


To support Boolean operations, we need a representation where the bits are placed in the complex slots

$$
\Delta (b_0 , b_1, \cdots, b_{n-1}) \in \mathbb{C}^{n}
$$

In this representation, bitwise AND/OR can be simulated as $(ab)$ and $(a + b - ab)$ over complex slots. Bit shifting can be simulated using CKKS rotation.

On the other hand, this representation resembles the "top row" (purple line) in Figure 1. Therefore, the central task of Arithmetic-to-Boolean conversion is to extract the bits of the triangle efficiently.
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

In summary, we introduce the triangle encoding, a new encoding for machine words. In arithmetic mode, each refresh uses two bootstrapping-like operations, independent of the word size. When we want to mix in Boolean operations, we can achieve amortized $O(1)$ bootstrapping when there are a sufficient number of ciphertexts.

Our scheme is most suitable for applications with heavy arithmetic and lightweight Boolean operations. The recommended word sizes for applications are $n=32$ and $n=64$, as these are common word sizes in modern computer. We also demonstrated the use of $n=256$ for high-precision arithmetic in the paper. In principle, we can support arbitrary bit-width $n$, but then the linear transformation cost between $$\mathbb{C}^{n/2}$$ and $$\mathbb{R}[X]/\langle X^n-X+2\rangle$$ grows in the order of $$O(\sqrt{n})$$, unlike in RadixCKKS where a DFT-like structure can be used.

This post focuses on the amortized conversion strategy presented in the paper. More direct arithmetic-to-Boolean conversion methods, without relying on the same amortization regime, will be discussed separately.

<br />
## References
<div style="margin-top: 1.5em;"></div>

[1] Jeffrey Hoffstein and Joseph Silverman, "Optimizations for NTRU". Proc. the Conf. on Public Key Cryptography and Computational Number Theory, Warsaw. pp. 77–88 (2000)

[2] Hao Chen, Kim Laine, Rachel Player, and Yuhou Xia, "High-Precision Arithmetic in Homomorphic Encryption." CT-RSA'18

[3] Robin Geelen and Frederik Vercauteren, "Fully Homomorphic Encryption for Cyclotomic Prime Moduli." EUROCRYPT'25

[4] Hyunho Cha, Intak Hwang, Seonhong Min, Jinyeong Seo, and Yongsoo Song, "MatriGear: Accelerating Authenticated Matrix Triple Generation with Scalable Prime Fields via Optimized HE Packing." IEEE S&P'25

[5] Zvika Brakerski, Offir Friedman, Daniel Golan, Alon Gurny, Dolev Mutzari, and Ohad Sheinfeld, "REFHE: Fully Homomorphic ALU", EUROCRYPT'26

[6] Gyeongwon Cha, Dongjin Park, and Joon-Woo Lee, "Improved Radix-based Approximate Homomorphic Encryption for Large Integers via Lightweight Bootstrapped Digit Carry", EUROCRYPT'26

[7] Jaehyung Kim and Taeyeong Noh, "Modular Reduction in CKKS" CIC'25

[8] Youngjin Bae, Jung Hee Cheon, Wonhee Cho, Jaehyung Kim, and Taekyung Kim, "META-BTS: Bootstrapping Precision Beyond the Limit" CCS'22
