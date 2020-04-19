/*****************************************
 Author: Adam Lenart
 File: skewnormal.js

 Provides skew normal pdf and cdf and
 truncated skew normal pdf 
******************************************/

function erf(x) {
    /**
      * Returns the value of the error function by a rational approximation according to 
      * Abramowitz and Stegun formula 7.1.26 
      * https://stackoverflow.com/questions/14846767/std-normal-cdf-normal-cdf-or-error-function
      */
    let sign = (x >= 0) ? 1 : -1;
    let z = Math.abs(x);   
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;
    let t = 1.0/(1.0 + p*x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-Math.pow(z,2));
    return sign * y; // erf(-x) = -erf(x);
};


function standardNormalPdf(x) {
    /**
      * Returns the standard normal pdf evaluated at x
     */
    let phi = 1/( Math.sqrt(2*Math.PI) ) * Math.exp(-Math.pow(x,2)/2);
    return phi;
};

function standardNormalCdf(x) {
    /**
      * Returns the standard normal cdf evaluated at x
      */
    let Phi = 1/2 * ( 1 + erf(x/Math.sqrt(2)) );
    return Phi;
};

function skewNormalPdf(x, loc, scale, shape) {
    /**
      * Returns the skew normal pdf evaluated at x, for shape=0 it is the same as the normal
      * distribution, and for loc=1 and scale=1, it reduces to the standard normal
      *
      * @param {number} x Evaluate pdf at x
      * @param {number} loc Location parameter of the skew normal distribution
      * @param {number} scale Scale parameter of the skew normal distribution, non-negative
      * @param {number} shape Shape parameter of the skew normal distribution
      * @returns Pdf of the skew normal distribution
      */
    let f1 = 2/scale * standardNormalPdf((x - loc)/scale);
    let f2 = standardNormalCdf(shape * (x - loc)/scale);
    return f1 * f2;
};

function skewNormalCdf(x, loc, scale, shape, n_grid=8) {
    /**
     * Returns the cdf of a skew normal evaluated approximately at n_grid points.
     * @param {number} x Evaluate cdf at point x
     * @param {number} loc Location parameter of the distribution
     * @param {number} scale Scale parameter of the distribution, non-negative
     * @param {number} shape Shape parameter of the distribution
     * @param {number} n_grid Number of grid points to evaluate Owen's T function
     * by the trapezoidal rule
     * @returns {number} Cdf of the distribution
     */
    let Phi = standardNormalCdf( (x - loc / scale) );
    let T = OwensT( (x - loc) / scale , shape, n_grid);
    return Phi - 2 * T;
};

function truncatedSkewNormalPdf(x, loc, scale, shape, lower=0, upper=1, n_grid=8) {
    /**
     * Returns a truncated skew normal pdf
     *
     * @param {number} x Evaluate pdf at point x
     * @param {number} loc Location parameter of the distribution
     * @param {number} scale Scale parameter of the distribution, non-negative
     * @param {number} shape Shape parameter of the distribution
     * @param {number} lower Lower value of truncation, x in [lower, upper]
     * @param {number} upper Upper value of trancuation, x in [lower, upper]
     * @param {number} n_grid Number of grid points to evaluate Owen's T function
     * in skewNormalCdf by the trapezoidal rule, defaults to 8
     * @returns {number} pdf of the distribution
       */
    let phi = skewNormalPdf(x, loc, scale, shape);
    let F_lower = skewNormalCdf(lower, loc, scale, shape, n_grid);
    let F_upper = skewNormalPdf(upper, loc, scale, shape, n_grid);
    return phi / (F_upper - F_lower);
};

function OwensTh0(a) {
    /**
      * Owen's T function for h=0
      *
      * The truncated skew normal pdf includes Owen's T function
      * T(h,a) = 1/(2*pi) * int_0^a( exp(-1/2*h^2 * ( 1 + x^2 ) ) )/( 1 + x^2 ) dx
      * Here is the evaluation for when T(h,a) := T(0,a), that is at point 0.
      * @param {number} a Parameter a of Owen's T function
      * @returns Value of Owen's T function for T(0,a)  
      */
    return 1/(2*Math.PI) * Math.atan(a);
};

function OwensT(h,a,n_grid) {
    /**
      * Owen's T function, T(h,a) approximated by the trapezoidal rule
      *
      * @param {number} h Parameter h in T(h,a)
      * @param {number} a Parameter a in T(h,a)
      * @param {number} n_grid Number of grid points to evaluate on by the trapezoidal rule
      * @returns The value of Owen's T(h,a)
      */
    // The truncated skew normal pdf includes Owen's T function
    // T(h,a) = 1/(2*pi) * int_0^a( exp(-1/2*h^2 * ( 1 + x^2 ) ) )/( 1 + x^2 ) dx
    if (a === 0) {
        return 0;
    };
    if (h === 0) {
        return OwensTh0(a);
    };
    return 1/( 2 * Math.PI ) * trapezoidal(auxOwensT, h, 0, a, n_grid);
}

function auxOwensT(x, h) {
    // auxiliary function for calculating Owen's T function
    return Math.exp(-1/2 * Math.pow(h,2) * ( 1 + Math.pow(x,2) ) )/( 1 + Math.pow(x,2) )
}

function trapezoidal(f,h,a,b,n) {
    // Performs trapezoidal rule for approximate value of definite integral for function f
    // with extra argument h
    // ported from
    // https://www.geeksforgeeks.org/trapezoidal-rule-for-approximate-value-of-definite-integral/
    // set grid spacing
    let g = (b - a) / n;
    // Computing sum of first and last terms 
    let s = (f(a,h) + f(b,h));     
    let i = 1;
    while (i < n) {
        // adding sum of middle terms
        s = s + 2 * f(a + i * g,h); 
        i ++;
    };
    return ((g / 2) * s);
}; 

export function skewNormalPdf;
export function skewNormalCdf;
export function truncatedSkewNormalPdf;


