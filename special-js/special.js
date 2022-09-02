
const special = {
    //constants
    gammaP: [676.5203681218851
        ,-1259.1392167224028
        ,771.32342877765313
        ,-176.61502916214059
        ,12.507343278686905
        ,-0.13857109526572012
        ,9.9843695780195716e-6
        ,1.5056327351493116e-7
    ],
    //first 1000 zeros of the zeta function
    //Complex functions
    complex: function(real, imaginary){
        return new special.Complex(real, imaginary);
    },
    add: function(a, b){
        //add two complex numbers
        if (a instanceof special.Complex == false){
            a = new special.Complex(a, 0);
        }
        if (b instanceof special.Complex == false){
            b = new special.Complex(b, 0);
        }
        return new special.Complex(a.real + b.real, a.imaginary + b.imaginary);
    },
    sub: function(a, b){
        //subtract two complex numbers
        if (a instanceof special.Complex == false){
            a = new special.Complex(a);
        }
        if (b instanceof special.Complex == false){
            b = new special.Complex(b);
        }
        return new special.Complex(a.real - b.real, a.imaginary - b.imaginary);
    },
    mult: function(a, b){
        //multiply two complex numbers
        if (a instanceof special.Complex == false){
            a = new special.Complex(a);
        }
        if (b instanceof special.Complex == false){
            b = new special.Complex(b);
        }
        return new special.Complex(a.real*b.real - a.imaginary*b.imaginary, a.real*b.imaginary + a.imaginary*b.real);
    },
    div: function(a, b){
        //divide two complex numbers
        if (a instanceof special.Complex == false){
            a = new special.Complex(a);
        }
        if (b instanceof special.Complex == false){
            b = new special.Complex(b);
        }
        return special.mult(a, b.get_conjugate()).div(b.magnitude()**2);
    },
    exp: function(a){
        //exponential of a complex number
        if (a instanceof special.Complex == false){
            return Math.exp(a)
        }
        return new special.Complex(Math.exp(a.real)*Math.cos(a.imaginary), Math.exp(a.real)*Math.sin(a.imaginary));
    },
    sqrt: function(a){
        //square root of a complex number
        if (a instanceof special.Complex == false && a >= 0){
            return Math.sqrt(a)
        }
        if (a instanceof special.Complex == false){
            a = special.complex(a, 0)
        }
        
        if (a.imaginary == 0){
            if (a.real < 0){
                return special.complex(0, Math.sqrt(-a.real))
            }
            return special.complex(Math.sqrt(a.real), 0)
        }
        return new special.Complex(Math.sqrt((a.real + a.magnitude())/2), a.imaginary/Math.abs(a.imaginary)*Math.sqrt((-a.real + a.magnitude())/2));
    },
    polar: function(a){
        //polar form of a complex number
        if (a instanceof special.Complex == false){
            return new special.Complex(a, 0)
        }
        return new special.Complex(a.magnitude(), Math.atan2(a.imaginary, a.real));
    },
    cartesian: function(a){
        //cartesian form of a complex number
        if (a instanceof special.Complex == false){
            return new special.Complex(a, 0)
        }
        return new special.Complex(a.real*Math.cos(a.imaginary), a.real*Math.sin(a.imaginary));
    },
    pow: function(a_, b_){
        //power of a complex number
        if (a_ instanceof special.Complex == false){
            a_ = special.complex(a_)
        }
        if (b_ instanceof special.Complex == false){
            b_ = special.complex(b_)
        }
        let a = a_.real
        let b = a_.imaginary
        let c = b_.real
        let d = b_.imaginary

        let f = (a**2 + b**2) ** (c/2) * Math.exp(-d * special.arg(a_))
        let g = Math.cos(c * special.arg(a_) + d * Math.log(a**2 + b**2) / 2)
        let h = Math.sin(c * special.arg(a_) + d * Math.log(a**2 + b**2) / 2)
        return special.mult(f, special.complex(g, h))
    },
    abs: function(a){
        //absolute value of a complex number
        if (a instanceof special.Complex == false){
            return Math.abs(a)
        }
        return a.magnitude();
    },
    arg: function(a){
        //argument of a complex number
        if (a instanceof special.Complex == false){
            a = new special.Complex(a);
        }
        return Math.atan2(a.imaginary, a.real);
    },
    log: function(a){
        //natural log of a complex number
        if (a instanceof special.Complex == false){
            return Math.log(a)
        }
        return new special.Complex(Math.log(a.magnitude()), special.arg(a));
    },
    //Trigonometric functions
    sin: function(a){
        //sine of a complex number
        if (a instanceof special.Complex == false){
            return Math.sin(a)
        }
        return new special.Complex(Math.sin(a.real)*Math.cosh(a.imaginary), Math.cos(a.real)*Math.sinh(a.imaginary));
    },
    cos: function(a){
        //cosine of a complex number
        if (a instanceof special.Complex == false){
            return Math.cos(a)
        }
        return new special.Complex(Math.cos(a.real)*Math.cosh(a.imaginary), -Math.sin(a.real)*Math.sinh(a.imaginary));
    },
    tan: function(a){
        //tangent of a complex number
        if (a instanceof special.Complex == false){
            return Math.tan(a)
        }
        return special.sin(a).div(special.cos(a));
    },
    sinh: function(a){
        //hyperbolic sine of a complex number
        if (a instanceof special.Complex == false){
            return Math.sinh(a)
        }
        return new special.Complex(Math.sinh(a.real)*Math.cos(a.imaginary), Math.cosh(a.real)*Math.sin(a.imaginary));
    },
    cosh: function(a){
        //hyperbolic cosine of a complex number
        if (a instanceof special.Complex == false){
            return Math.cosh(a)
        }
        return new special.Complex(Math.cosh(a.real)*Math.cos(a.imaginary), Math.sinh(a.real)*Math.sin(a.imaginary));
    },
    tanh: function(a){
        //hyperbolic tangent of a complex number
        if (a instanceof special.Complex == false){
            return Math.tanh(a)
        }
        return special.sinh(a).div(special.cosh(a));
    },
    asin: function(a){
        //arc sine of a complex number
        if (a instanceof special.Complex == false){
            return Math.asin(a)
        }
        let sq = special.sqrt(new special.Complex(1, 0).sub(a.mult(a)));
        return new special.Complex(0, -1).mult(special.log( special.mult(special.i, a).add(sq)  ))
    },
    acos: function(a){
        //arc cosine of a complex number
        if (a instanceof special.Complex == false){
            return Math.acos(a)
        }
        let sq = special.sqrt(new special.Complex(1, 0).sub(a.mult(a)));
        return new special.Complex(0, 1).mult(special.log( special.mult(special.i, a).add(sq)  )).add(Math.PI/2)

    },
    atan: function(a){
        //arc tangent of a complex number
        if (a instanceof special.Complex == false){
            return Math.atan(a)
        }
        let f = special.mult(special.i, special.log(special.sub(1, special.mult(special.i, a)))).div(2)
        let s = special.mult(special.i, special.log(special.add(1, special.mult(special.i, a)))).div(2)
        return special.sub(f, s)
    },
    //benchmarking functons
    benchmark: function(func, args){
        //benchmark a function
        let start = performance.now()
        console.log("function output:", func(...args))
        console.log("Time taken: " + (performance.now() - start) + "ms")
    },
    //Special functions (the good stuff)
    gcd: function(a, b){
        //greatest common divisor
        let x = a
        let y = b
        while (y != 0){
            let temp = y
            y = x % y
            x = temp
        }
        return x
    },
    lcm: function(a, b){
        //least common multiple
        return a*b/special.gcd(a, b)
    },
    isPrime: function(n){
        //check if a number is prime
        if (n == 2){
            return true
        }
        if (n % 2 == 0){
            return false
        }
        for (let i = 3; i <= Math.sqrt(n); i+=2){
            if (n % i == 0){
                return false
            }
        }
        return true
    },
    mobius: function(n){
        //mobius function
        if (n == 1){
            return 1
          }
          if (n == 2){
            return -1
          }
          let p = 0
        
          if (n % 2 == 0){
            n = Math.floor(n/2)
            p++
        
            if (n%2==0){
              return 0
            }
          }
          for (let i=3;i<Math.floor(Math.sqrt(n))+1;i++){
            if (n%i == 0){
              n = Math.floor(n/i)
              p++
              if (n%i==0){
                return 0
              }
            }
            i+=2
          }
          if (p%2 == 0){
            return -1
          }
          else{
            return 1
          }
    },
    nCr: function(n, r){
        //n choose r
        if (n instanceof special.Complex || r instanceof special.Complex || Math.floor(n) != n || Math.floor(r) != r){
            if (n instanceof special.Complex == false){
                n = special.complex(n, 0)
            }
            if (r instanceof special.Complex == false){
                r = special.complex(r, 0)
            }
            if (special.abs(n.sub(r).add(1)) == 0){
                return special.complex(0,0)
            }
            let o = special.gamma(n.add(1)).div(special.gamma(r.add(1)).mult(special.gamma(n.sub(r).add(1))))
            return o
        }
        return special.factorial(n) / (special.factorial(r) * special.factorial(n-r))
    },
    gamma: function(a){
        //gamma function
        if (a instanceof special.Complex == false){
            a = new special.Complex(a, 0)
        }
        if (a.real < 0.5){
            let t = special.complex(1-a.real, -a.imaginary)
            let r = special.complex(Math.PI*a.real, Math.PI*a.imaginary)
            let y = special.complex(Math.PI, 0)
            return special.div(y, special.sin(r)).div(special.gamma(t))
        }
        let n = special.complex(a.real-1, a.imaginary)
        let x = special.complex(1, 0)
        for (let i=0;i<special.gammaP.length;++i){
            x = x.add(special.complex(special.gammaP[i], 0).div(n.add(i+1)))
        }
        let t = special.complex(n.real + special.gammaP.length - 0.5, n.imaginary)
        return x.mult(special.complex(Math.sqrt(2*Math.PI), 0)).mult(special.pow(t, special.add(n, 0.5))).mult(special.exp(special.complex(-t.real, -t.imaginary)))
    },
    bernoulli: function(n, prec){
        return special.mult(n, -1).mult(special.riemann_zeta(special.sub(1, n), prec))
    },
    factorial: function(n){
        //factorial
        if (n instanceof special.Complex || Math.floor(n) != n){
            let o = special.gamma(n+1)
            return o.imaginary == 0 ? o.real : o
        }
        let ans = 1
        for (let i=1;i<=n;i++){
            ans *= i
        }
        return ans
    },
    riemann_zeta: function(s, prec){
        if (s instanceof special.Complex == false){
            s = new special.Complex(s, 0)
        }
        if (s.real == 0 && s.imaginary == 0){
            return -0.5
        }
        prec = prec == undefined ? 1e-3 : prec
        let dec = -Math.round(Math.log(prec*0.1)/Math.log(10))
        let n = Math.round(1.3*dec + 0.9*Math.abs(s.imaginary))
        n = Math.min(n, 60)
        
        function d(k){
            let S = 0
            for (let j=k;j<=n;j++){
                S += (special.factorial(n+j-1) * (4**j) / (special.factorial(n-j)*special.factorial(2*j)  ))
            }
            return n*S
        }
        function f(s){
            let c = special.div(1, 
                special.mult(d(0), special.sub(1, special.pow(2, special.sub(1, s)))))
            let S = special.complex(0, 0)
            for (let k=1;k<=n;k++){
                S = S.add(
                    special.div((-1)**(k-1) * d(k), special.pow(k, s))
                )
            }
            return special.mult(c, S)
        }
        if (s.real > 1 || s.real == 0){
            console.log("Running "+ n +" iterations")
            return f(s)
        }
        else{
            let c = special.mult(special.pow(2, s), special.pow(Math.PI, special.sub(s, 1))).mult(special.sin(special.mult(Math.PI/2, s))).mult(special.gamma(special.sub(1, s)))
            return special.mult(c, special.riemann_zeta(special.sub(1, s)))
        }
    },
    prime_counting: function(x){
        //number of primes less than x
        let ans = 0
        for (let i=2;i<x;i++){
            if (special.isPrime(i)){
                ans++
            }
        }
        return ans
    },
    beta: function(a, b){
        return special.div(special.mult(special.gamma(a), special.gamma(b)), special.gamma(special.add(a, b)))
    },
    lower_incomplete_gamma: function(s, z, prec){
        if (s instanceof special.Complex == false){
            s = new special.Complex(s, 0)
        }
        if (z instanceof special.Complex == false){
            z = new special.Complex(z, 0)
        }
        prec = prec == undefined ? 1e-3 : prec
        if (z.real == 0 && z.imaginary == 0){
            return special.gamma(s)
        }

        let S = special.complex(0, 0)
        for (let k=0;k<=100;k++){
            let old = special.complex(S.real, S.imaginary)

            let re = (-1)**k / special.factorial(k)
            let im = special.div(special.pow(z, special.add(s, k)), special.add(s, k))
            S = S.add(special.mult(re, im))
            if (special.abs(special.sub(S, old)) < prec){
                break
            }
        }
        return S
    },
    incomplete_gamma: function(s, z, prec){
        if (s instanceof special.Complex == false){
            s = new special.Complex(s, 0)
        }
        if (z instanceof special.Complex == false){
            z = new special.Complex(z, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        return special.sub(special.gamma(s), special.lower_incomplete_gamma(s, z, prec))
    },
    Li: function(x, prec){
        if (x instanceof special.Complex == false){
            x = special.complex(x, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        let S = special.complex(0, 0)
        for (let n=1;n<=100;n++){
            let old = special.complex(S.real, S.imaginary)
            let S2 = 0
            for (let k=0;k<=Math.floor((n-1)/2);k++){
                S2 += 1/(2*k+1)
            }
            let re = (-1) ** (n-1) / (special.factorial(n)* 2**(n-1))
            let im = special.pow(special.log(x), n)

            S = S.add(special.mult(re, im).mult(S2))
            if (special.abs(special.sub(S, old)) < prec){
                break
            }
        }
        return special.add(0.577215664901532, special.log(special.log(x))).add(special.mult(special.sqrt(x), S))
    },
    Ei: function(x, prec){
        if (x instanceof special.Complex == false){
            x = special.complex(x, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        let S = special.complex(0, 0)
        for (let n=1;n<=100;n++){
            let old = special.complex(S.real, S.imaginary)
            let S2 = 0
            for (let k=0;k<=Math.floor((n-1)/2);k++){
                S2 += 1/(2*k+1)
            }
            let re = (-1) ** (n-1) / (special.factorial(n)* 2**(n-1))
            let im = special.pow(x, n)

            S = S.add(special.mult(re, im).mult(S2))
            if (special.abs(special.sub(S, old)) < prec){
                break
            }
        }
        return special.add(0.577215664901532, special.log(x)).add(special.mult(special.exp(special.div(x, 2)), S))
    },
    Si: function(x, prec){
        if (x instanceof special.Complex == false){
            x = special.complex(x, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        let S = special.complex(0, 0)
        for (let n=0;n<=100;n++){
            let old = special.complex(S.real, S.imaginary)
            let S2 = 0
            let re = (-1) ** n / (special.factorial(2*n+1)* (2*n+1))
            let im = special.pow(x, 2*n+1)

            S = S.add(special.mult(re, im))
            if (special.abs(special.sub(S, old)) < prec){
                break
            }
        }
        return S
    },
    Ci: function(x, prec){
        if (x instanceof special.Complex == false){
            x = special.complex(x, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        let S = special.complex(0, 0)
        for (let n=1;n<=100;n++){
            let old = special.complex(S.real, S.imaginary)
            let re = (-1) ** n / (special.factorial(2*n)* (2*n))
            let im = special.pow(x, 2*n)

            S = S.add(special.mult(re, im))
            if (special.abs(special.sub(S, old)) < prec){
                break
            }
        }
        return special.add(0.577215664901532, special.log(x)).add(S)
    },
    erf: function(x, prec){
        if (x instanceof special.Complex == false){
            x = special.complex(x, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        return special.mult(Math.PI ** (-0.5), special.lower_incomplete_gamma(0.5, special.pow(x, 2), prec))
    },
    LambertW: function(x, prec){
        if (x instanceof special.Complex == false){
            x = special.complex(x, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        let e_x = special.mult(Math.E, x)
        let sqrt_e_x = special.sqrt(special.add(1, e_x))
        //approximate guess
        let g = special.div(special.mult(e_x, special.log(special.add(1, sqrt_e_x))), special.add(1, e_x).add(sqrt_e_x))
        if (isNaN(g.real) || isNaN(g.imaginary)){
            g = special.log(x)
        }
        //halleys interation
        for (let i=0;i<=500;i++){
            let old = special.complex(g.real, g.imaginary)
            let e_g = special.exp(g)
            let inner_frac = special.div(special.mult(special.add(g, 2), special.sub(special.mult(g, e_g), x)), special.add(special.mult(g, 2), 2)) 
            let full_frac = special.div(special.sub(special.mult(g, e_g), x), special.sub(special.mult(e_g, special.add(g, 1)), inner_frac))
            g = g.sub(full_frac)
            if (special.abs(special.sub(g, old)) < prec){
                console.log(i)
                break
            }
        }
        return g
    },
    AGM: function(a, b, prec){
        if (a instanceof special.Complex == false){
            a = new special.Complex(a, 0)
        }
        if (b instanceof special.Complex == false){
            b = new special.Complex(b, 0)
        }
        prec = prec == undefined ? 1e-3 : prec

        for (let i=0;i<=500;i++){
            let old_a = special.mult(0.5, special.add(a, b))
            b = special.sqrt(special.mult(a, b))
            a = old_a

            if (special.abs(special.sub(a, b)) < prec){
                break
            }
        }
        return a
    }
    
    
    

    
}

special.Complex = function(real, imaginary) {
    this.real = real == undefined ? 0 : real;
    this.imaginary = imaginary == undefined ? 0 : imaginary;
    
    this.get_conjugate = function(){
        return new special.Complex(this.real, -this.imaginary);
    }
    this.conjugate = function(){
        this.imaginary = -this.imaginary;
    }
    this.magnitude = function(){
        return Math.sqrt(this.real**2 + this.imaginary**2);
    }
    this.add = function(other){
        if (other instanceof special.Complex == false){
            other = new special.Complex(other);
        }
        return new special.Complex(this.real + other.real, this.imaginary + other.imaginary);
    }
    this.sub = function(other){
        if (other instanceof special.Complex == false){
            other = new special.Complex(other);
        }
        return new special.Complex(this.real - other.real, this.imaginary - other.imaginary);
    }
    this.mult = function(other){
        if (other instanceof special.Complex == false){
            other = new special.Complex(other);
        }
        return new special.Complex(this.real*other.real - this.imaginary*other.imaginary, this.real*other.imaginary + this.imaginary*other.real);
    }
    this.div = function(other){
        if (other instanceof special.Complex == false){
            other = new special.Complex(other);
        }
        let c = this.mult(other.get_conjugate())
        c.real /= other.magnitude()**2;
        c.imaginary /= other.magnitude()**2;
        return c;
    }
}

special.i = new special.Complex(0, 1);
//Old non functional stuff

// prime_counting: function(x, z, iter){
    //     function SUM(f,l,u){
    //         let s = 0
    //         for (let i=l;i<=u;i++){
    //           s += f(i)
    //         }
    //         return s;
    //     }
    //     function li(x){
    //         function first(n){
    //           return (Math.pow(-1,n-1)*Math.pow(Math.log(x),n))/(special.factorial(n)*Math.pow(2,n-1))
    //         }
    //         function second(n){
    //           return 1/(2*n+1)
    //         }
          
    //         let S = SUM(function(n){
    //           return first(n)*SUM(second,0,Math.floor((n-1)/2))
    //         },0,50)
    //         return 0.57721 + Math.log(Math.log(x))+Math.sqrt(x)*S
    //     }
    //     function f(x,ZS){
    //         let S = 0
    //         for (let i=0;i<ZS;i++){
    //             let Z = special.complex(0.5,special.zeta_zeros[i])
    //             let EZ = special.sub(1,special.complex(0.5,special.zeta_zeros[i]))
    //             let a = special.div(special.exp(special.mult(Z,Math.log(x))),  special.mult(Z,Math.log(x)) )
    //             let b = special.div(special.exp(special.mult(EZ,Math.log(x))),  special.mult(EZ,Math.log(x)) )
    //             S = special.add(S,special.add(a,b))
    //         }
    //         S = S.real
    //         if (ZS == 0){
    //             S = 0
    //         }
            
    //         let I = NumCalc.NIntegrate(function(t){
    //             return 1/(t*(t**2-1)*Math.log(t))
    //         },x,Infinity)
    //         return li(x) - S - Math.log(2) + I
    //     }
    //     let S = 0
    //     for (let i=1;i<iter;i++){
    //         S += special.mobius(i)/i * f(x**(1/i), z)
    //     }
    //     return Math.round(S)

    // },

// harmonic: function(x, n, prec){
    //     if (x instanceof special.Complex == false){
    //         x = special.complex(x,0)
    //     }
    //     prec = prec == undefined ? 1e-3 : prec

    //     if (n != undefined){
    //         if (n instanceof special.Complex == false){
    //             n = special.complex(n, 0)
    //         }
    //         let a = special.harmonic(x.add(n).sub(1))
    //         let b = special.harmonic(n.sub(1))
            
    //         return special.nCr(x.add(n).sub(1), n.sub(1)).mult(special.sub(a, b))
    //     }

    //     let s = special.complex(0, 0)
    //     for (let k=1;k<=50;k++){
    //         let old = special.complex(s.real, s.imaginary)
    //         s = special.add(s, 
    //             special.nCr(x, k).mult(special.div((-1)**k, k))
    //         )
    //         if (special.abs(special.sub(s, old)) <= prec){
    //             break
    //         }
    //     }
    //     return s.mult(-1)

        
    // },
    // polygamma: function(n, x, prec){
    //     if (x instanceof special.Complex == false){
    //         x = special.complex(x, 0)
    //     }
    //     if (n instanceof special.Complex == false){
    //         n = special.complex(n, 0)
    //     }
    //     let numer = special.mult(special.pow(-1, n.add(1)), special.gamma(n.add(1)))
    //     let zeta = special.riemann_zeta(n.add(1))
    //     let harm = special.pow(special.harmonic(x.sub(1)), n.add(1))

    //     return special.mult(numer, special.sub(zeta, harm))
    // }