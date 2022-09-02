# special.js
A collection of special math functions written in JavaScript with full complex support! I'm looking for contributors to help test and more functionality to the library so feel free to make an issue or puul request with you suggestions!

## Usage
Grab the `special.js` or `special.min.js` from the repo and add it to your project to get started

# Documentation
Documentation for the library lives here. All the functions are in a namespace and can be used using the `special.` prefix

# Complex Numbers
## Complex(a: Number, b: Number)
A constructor used to create a complex number. Can be used without the `new` keyword using `special.complex(a, b)`

## Complex.add(a: Complex)
Returns the sum of a complex number

## Complex.sub(a: Complex)
Returns the difference of a complex number

## Complex.mult(a: Complex)
Returns the product of a complex number

## Complex.div(a: Complex)
Returns the quotient of a complex number

## Complex.get_conjugate()
Returns the conjugate of a complex number

## Complex.magnitude()
Returns the magnitude (absolute value) of a complex number


## add(a: Complex, b: Complex)
Returns the sum of 2 complex numbers

## sub(a: Complex, b: Complex)
Returns the difference of 2 complex numbers

## mult(a: Complex, b: Complex)
Returns the product of 2 complex numbers

## div(a: Complex, b: Complex)
Returns the quotient of 2 complex numbers

## exp(a: Complex)
Returns the exponential function of a complex number

## sqrt(a: Complex, b: Complex)
Returns the square root of a complex number

## polar(a: Complex)
Returns the polar form of a cartesian complex number

## cartesian(a: Complex)
Returns the cartesian form of a polar complex number

## pow(a: Complex, b: Complex)
Returns the power of a complex number to another

## abs(a: Complex)
Returns the absolute value of a complex number

## arg(a: Complex)
Returns the argument of a complex number

## log(a: Complex)
Returns the natural logarithm of a complex number

## sin(a: Complex)
Returns the `sin` of a complex number

## cos(a: Complex)
Returns the `cos` of a complex number

## tan(a: Complex)
Returns the `tan` of a complex number

## sinh(a: Complex)
Returns the `sinh` of a complex number

## cosh(a: Complex)
Returns the `cosh` of a complex number

## tanh(a: Complex)
Returns the `tanh` of a complex number

## asin(a: Complex)
Returns the `arcsin` of a complex number

## acos(a: Complex)
Returns the `arccos` of a complex number

## atan(a: Complex)
Returns the `arctan` of a complex number

# Special Functions

## gcd(a: Integer, b: Integer)
Returns the `gcd` of a 2 integers

## lcm(a: Integer, b: Integer)
Returns the `lcm` of a 2 integers

## isPrime(a: Integer)
Returns a boolean value determining the primality of an integer

## mobius(a: Integer)
Returns the mobius `mu` function of an integer

## nCr(a: Complex, b: Complex)
Returns `nCr` of 2 numbers making use of the gamma function to extend the definition to all of the complex plane

## gamma(a: Complex)
Returns the gamma function of a complex number defined by $\displaystyle \Gamma(z)=\int_0^\infty e^{-t}t^{z-1}dt$ making use of the Lanczos Approximation. Precision is fixed.

## bernoulli(n: Complex)
Returns the nth Bernoulli number making use of the equation $B_n=-n\zeta(1-n)$ where $\zeta(s)$ is the Riemann Zeta Function.

## factorial(a: Integer)
Returns the factorial of a positive integer

## riemann_zeta(s: Complex, prec: Number)
Returns the analytic continuation of the Riemann Zeta Function for the entire complex plane aside from $s=1$. Precision can be set and is `1e-3` at default. Computation makes use of the function $\displaystyle d_k = n\sum_{j=k}^n \frac{(n+j-1)!4^j}{(n-j)!(2j)!}$ where $n$ is chosen arbitrarily making use of the precision and uses $\displaystyle\zeta(s) = \frac{1}{d_0(1-2^{1-s})}\sum_{k=1}^n \frac{(-1)^{k-1}d_k}{k^s}$ for $\Re(s) > 1$ and uses the functional equation $\displaystyle \zeta(s) = 2^s\pi^{s-1}\sin\bigg(\frac{\pi s}{2}\bigg)\Gamma(1-s)\zeta(1-s)$ where $\Gamma(s)$ is the gamma function to extend it to the entire complex plane.

## beta(a: Complex, b: Complex)
Returns the beta function of 2 complex numbers defined by $\displaystyle B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}$

## lower_incomplete_gamma(s: Complex, z: Complex, prec: Number)
Returns the lower incomplete gamma function defined by $\displaystyle \int_0^z t^{s-1}e^{-t}dt$ to a specified precision  with default `1e-3`

## incomplete_gamma(s: Complex, z: Complex, prec: Number)
Returns the incomplete gamma function defined by $\displaystyle \int_x^\infty t^{s-1}e^{-t}dt$ using $\displaystyle \gamma(s,z)+\Gamma(s,z) = \Gamma(s)$ where $\gamma(s,z)$ is the lower incomplete gamma function computed to a specified precision with default `1e-3`

## Li(s: Complex, prec: Number)
Returns the Logarithmic Integral function defined by $\displaystyle \int_0^x \frac{dt}{\ln(x)}$ using a sum given by Ramanujan defined by $\displaystyle Li(x)= \gamma + \ln\ln x + \sqrt{x}\sum_{n=1}^\infty \frac{(-1)^{n-1}(\ln x)^n}{n! \ 2^{n-1}} \sum_{k=0}^{\lfloor \frac{n-1}{2}\rfloor} \frac{1}{2k+1}$ where $\gamma$ is the Euler-Mascheroni constant $0.57221...$. Precision is set to `1e-3` by default

## Ei(s: Complex, prec: Number)
Returns the Exponential Integral function defined by $\displaystyle \int_{-\infty}^x \frac{e^t}{t}dt$ using a sum by Ramanujan defined by $\displaystyle Ei(x)= \gamma + \ln x + \exp{\frac{x}{2}}\sum_{n=1}^\infty \frac{(-1)^{n-1}(x)^n}{n! \ 2^{n-1}} \sum_{k=0}^{\lfloor \frac{n-1}{2}\rfloor} \frac{1}{2k+1}$ where $\gamma$ is the Euler-Mascheroni constant $0.57221...$. Precision is set to `1e-3` by default

## Si(s: Complex, prec: Number)
Returns the Sine Integral function defined by $\displaystyle \int_{0}^x \frac{\sin t}{t}dt$ using a sum by Ramanujan defined by $\displaystyle Si(x)=\sum_{n=0}^\infty \frac{(-1)^n x^{2n+1}}{(2n+1)(2n+1)!}$ Precision is set to `1e-3` by default

## Ci(s: Complex, prec: Number)
Returns the Cosine Integral function defined by $\displaystyle \int_{0}^x \frac{\cos t}{t}dt$ using a sum by Ramanujan defined by $\displaystyle Si(x)=\gamma + \ln x + \sum_{n=0}^\infty \frac{(-1)^n x^{2n}}{(2n)(2n)!}$ where $\gamma$ is the Euler-Mascheroni constant $0.57221...$. Precision is set to `1e-3` by default

## erf(s: Complex, prec: Number)
Returns the Error function defined by $\displaystyle \frac{2}{\sqrt \pi}\int_{0}^x e^{-t^2}dt$ using and indentity defined by $\displaystyle \text{erf}(x) = \pi ^ {-\frac{1}{2}} \gamma\bigg(\frac{1}{2}, x^2\bigg)$ where $\gamma$ is the lower incomplete gamma function. Precision is set to `1e-3` by default

## LambertW(x: Complex, prec: Number)
Returns the Lambert W function defined by begin the inverse of $xe^x$ using the Halley's Iteration method $\displaystyle w_{j+1} = w_j - \frac{w_je^{w_j}-z}{e^{w_j}(w_j+1)-\frac{(w_j+2)(w_je^{w_j}-z)}{2w_j+2}}$ with a default precision of `1e-3`

## AGM(a: Complex, b: Complex, prec: Number)
Returns the Arithmetic-Geometric Mean of 2 numbers using the iteration formula $\displaystyle a_{n+1} = \frac{1}{2}(a_n + b_n)$ and $\displaystyle b_{n+1} = \sqrt{a_nb_n}$ until the difference between $a_n$ and $b_n$ is less than the precision. Default precision is `1e-3`




