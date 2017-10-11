# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 16:07:42 2017

@author: Feng Li
"""
from math import *
TOL = 1.0E-15
NMAX = 100

def sign(x, y):
    r"""
    Return the magnitude of x times the sign of y.
    """
    if y == 0.0:
        return 0.0
    else:
        _y = y / abs(y)
        return abs(x) * _y

def sgn(x):
    r"""
    The signum function of a real number x is defined as follows:
                  ┌  -1   if x < 0
        sgn(x) := ┤   0   if x = 0 
                  └   1   if x > 0
                  
                  ┌   0   if x = 0
                = ┤
                  └ |x|/x  if x != 0 
    
    Reference:
        https://en.wikipedia.org/wiki/Sign_function
    """
    if x == 0.0:
        return 0.0
    else:
        return abs(x)/x
        
def bisection(f, a, b, epsilon=TOL,LoopMax=NMAX):
    r"""
    Bisection method (二分法) solving the function f(x) = 0 at cloased
    interval [a,b] where f(a)*f(b) < 0 and f(x) is monotonous
    at this closed interval.
    
    Reference:
        [1] https://en.wikipedia.org/wiki/Bisection_method
        [2] William et al., Numerical Recipes in C The Art of Scientific 
            Computing, 2nd Edition, 2002: 354-357.
    """
    if a >= b:
        b, a = a, b
    
    fa, fb = f(a), f(b)
    if fa*fb >= 0:
        raise ValueError("End points values have the same sign.")
    else:
#        for N in range(LoopMax):
#            c = (a + b)/2.0
#            fa, fc = f(a), f(c)
#            if fc == 0.0 or abs(b - a)/2.0 < epsilon:
#                print("Bisection method:%d iterations."%N)
#                return c
#            if sign(fc, fa) == fc:
#                # f(a) and f(c) have the same sign.
#                a = c
#            else:
#                b = c
        # Tranlate from C code in reference [2]
        fl, fmid = f(a), f(b)
        if fl*fmid >= 0.0:
            raise ValueError("Root must be bracketed for bisection method.")
        if fl < 0.0:
            dx = b - a
            rt = a
        else:
            dx = a - b
            rt = b
        for N in range(LoopMax):
            dx *= 0.5
            xmid = rt + dx
            fmid = f(xmid)
            if fmid <= 0.0:
                rt = xmid
            if abs(dx) < epsilon or fmid == 0.0:
                print("Bisection method:%d iterations."%N)
                return rt
        raise ValueError("Method failed, iterations exceed the max loop number. "+\
                         "Try other method instead or use a different interval.")

def secant(f,a,b,epsilon=TOL,LoopMax=NMAX):
    r"""
    Secant method (割线法) for root finding at closed interval [a,b] where
    f(a)*f(b) < 0.
    The recurrence relation of secant method is
                                x[n-1] - f(x[n-2])
    x[n] = x[n-1] - f(x[n-1])*----------------------
                               f(x[n-1] - f(x[n-2]))
                               
            x[n-2]*f(x[n-1]) - x[n-1]*f(x[n-2])      b*f(a) - a*f(b)
         = ------------------------------------- = -------------------
                   f(x[n-1] - f(x[n-2]))                f(a) - f(b)
    As for comparsion with other method, wikipedia says:
        'The secant method does not require that the root remain bracketed like 
        the bisection method does, and hence it does not always converge. The 
        false position method (or regula falsi) uses the same formula as the 
        secant method. However, it does not apply the formula on x[n-1] and 
        x[n-2], like the secant method, but on x[n-1] and on the last iterate 
        x[k] such that f(x[k]) and f(x[n-1]) have a different sign. This means 
        that the false position method always converges.
        The recurrence formula of the secant method can be derived from the 
        formula for Newton's method by using the finite difference 
        approximation.'
        
    Reference:
        [1] https://en.wikipedia.org/wiki/Root-finding_algorithm
        [2] https://en.wikipedia.org/wiki/Secant_method
        [3] William et al., Numerical Recipes in C The Art of Scientific 
        Computing, 2nd Edition, 2002: 354-357.
    """
    
#    for N in range(LoopMax):
#        fa, fb = f(a), f(b)
#        if fa != fb:
#            rt = (fa*b - fb*a) / (fa - fb)
#            if abs(b - a) < epsilon:
#                print("Secant method:%d iterations."%N)
#                return rt
#        else:
#            print("Secant method:%d iterations."%N)
#            return rt
#        a = b
#        b = rt

    # Code translated from C in reference [3]
    fl, fx = f(a), f(b)
    # Pick the bound with the smaller function value as the most 
    # recent guess.
    if abs(fl) < fx:
        rt, xl = a, b
        fl, fx = fx, fl
    else:
        xl, rt = a, b
    
    for N in range(LoopMax):
        dx = (xl - rt)*fx/(fx - fl) # Increment with respect to latest value.
        xl, fl = rt, fx
        rt += dx
        fx = f(rt)
        if abs(dx) < epsilon or fx == 0.0:
            print("Secant method:%d iterations."%N)
            return rt
    raise ValueError("Method failed, iterations exceed the max loop number. "+\
                     "Try other method instead or use a different interval.")

def falsePosition(f,a,b,epsilon=TOL,LoopMax=NMAX):
    r"""
    False position method (试位法) of finding the zero piont of a function at
    the closed interval [a,b] where f(a)*f(b) < 0.
    False position is similar to the secant method, except that, instead of 
    retaining the last two points, it makes sure to keep one point on either 
    side of the root.
    The recurrence relation of secant method is
                                x[n-1] - f(x[n-2])
    x[n] = x[n-1] - f(x[n-1])*----------------------
                               f(x[n-1] - f(x[n-2]))
                              
            x[n-2]*f(x[n-1]) - x[n-1]*f(x[n-2])
         = -------------------------------------
                   f(x[n-1] - f(x[n-2]))               
    
    Reference:
        https://en.wikipedia.org/wiki/False_position_method
        https://en.wikipedia.org/wiki/Root-finding_algorithm
        William et al., Numerical Recipes in C The Art of Scientific 
        Computing, 2nd Edition, 2002: 354-357.
    """
    # Translate from C language source code to python code
    fl, fh = f(a), f(b)
    if fl * fh >= 0.0:
        raise ValueError( "Function values at end points have the same sign,"+\
                         "try other method or use other interals instead.")
    # Identify the limits so that xl corresponds to the low side.
    if fl < 0.0:
        xl, xh = a, b
    else:
        xl, xh = b, a
        fl, fh = fh, fl
    dx = xh - xl
    for N in range(LoopMax):
        rt = xl + dx*fl/(fl - fh)
        fx = f(rt)
        if fx < 0.0:
            dif = xl - rt
            xl, fl = rt, fx
        else:
            dif = xh - rt
            xh, fh = rt, fx
        dx = xh - xl
        if abs(dif) < epsilon or fx == 0.0:
            print("False position method:%d iterations."%N)
            return rt
    
    raise ValueError("Method failed, iterations exceed the max loop number. "+\
                     "Try other method instead or use a different interval.")

falsi = falsePosition

def illinois(f,a,b,epsilon=TOL,LoopMax=NMAX):
    r"""
    Improved false position method - Illinois algorithm:
    When the new y-value has the same sign as the previous one, meaning that 
    the data point before the previous one will be retained, the Illinois 
    version halves the y-value of the retained data point.
    
    Reference:
        https://en.wikipedia.org/wiki/False_position_method
    """
    fa, fb = f(a), f(b)
    if fa*fb >= 0.0:
        raise ValueError("Function values at end points have the same sign, "+\
                         "try other method or use other interals instead.")
    side = 0
    for N in range(LoopMax):
        r = (fa*b - fb*a) / (fa - fb)
        if abs(b - a) < epsilon * abs(a + b):
            print("Illinois method:%d iterations."%N)
            break
        fr = f(r)
        
        if fr * fb > 0.0:
            # fr and fb have the same sign, copy r to b
            b, fb = r, fr
            if side == -1:
                fa /= 2.0
            side = -1
        elif fa * fr > 0:
            # fr and fa have the same sign, copy r to a
            a, fa = r, fr
            if side == +1:
                fb /= 2.0
            side = +1
        else:
            # fr * fa or fr * fb very small (looks like zero)
            break
    return r

def ridders(f,a,b,epsilon=TOL,LoopMax=NMAX):
    r"""
    Reference:
        [1] William et al., Numerical Recipes in C The Art of Scientific 
        Computing, 2nd Edition, 2002: 358-359.
    """
    # Translate directly from the C code in reference [1]
    UNUSED = -1.11e30
    
    fl, fh = f(a), f(b)
    if fl * fh < 0.0:
        xl, xh = a, b
        ans = UNUSED # Any highly unlikely value, to simplify logic below
        for N in range(LoopMax):
            xm = 0.5*(xl + xh)
            fm = f(xm)
            s = sqrt(fm*fm - fl*fh)
            if s == 0.0: return ans
            if fl >= fh:
                signs = 1.0
            else:
                signs = - 1.0
            xnew = xm + (xm - xl)*signs*fm/s
            if abs(xnew - ans) <= epsilon:
                print("Ridders' method:%d iterations."%N)
                return ans
            ans = xnew
            fnew = f(ans)
            if fnew == 0.0: return ans
            if sign(fm, fnew) != fm:
                xl, fl = xm, fm
                xh, fh = ans, fnew
            elif sign(fl, fnew) != fl:
                xh, fh = ans, fnew
            elif sign(fh, fnew) != fh:
                xl, fl = ans, fnew
            else:
                pass # Never get here.
            if abs(xh - xl) <= epsilon:
                print("Ridders' method:%d iterations."%N)
                return ans
    else:
        if fl == 0.0: return a
        if fh == 0.0: return b
        raise ValueError("Root must be bracketed in Ridders' method.")

def brent(f, x1, x2, epsilon=TOL,LoopMax=NMAX):
    r"""
    Van Wijngaarden–Dekker–Brent Method:
    Brent’s method is guaranteed (by Brent) to converge, so long as the 
    function can be evaluated within the initial interval known to contain 
    a root. It combines root bracketing, bisection, and inverse quadratic 
    interpolation to converge from the neighborhood of a zero crossing. It
    combines the sureness of bisection with the speed of a higher-order 
    method when appropriate. 
    
    Reference:
        [1] William et al., Numerical Recipes in C The Art of Scientific 
            Computing, 2nd Edition, 2002: 359-362.
    """
    # Tranlated directly from C code in reference [1]
    a, b, c = x1, x2, x2
    fa, fb = f(a), f(b)
    
    if fa * fb > 0.0:
        raise ValueError("Root must be bracketed in Brent's method.")

    fc = fb
    for N in range(LoopMax):
        if fb * fc > 0.0:
            # Rename a,b,c and adjust bounding interval d
            c, fc = a, fa
            e = d = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, a
            fa, fb, fc = fb, fc, fa
        tol1 = 2*epsilon*abs(b) + 0.5*TOL # Convergence check.
        xm = 0.5*(c - b)
        if abs(xm) <= tol1 or fb == 0.0:
            print("Brent's method:%d iterations."%N)
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            # Attempt inverse quadratic interpolation.
            s = fb/fa
            if a == c:
                p, q = 2.0*xm*s, 1.0 - s
            else:
                q, r = fa/fc, fb/fc
                p = s * (2.0 * xm * q * (q - r) - (b - a)*(r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            # Check whether in bounds.
            if p > 0.0: q = -q
            p = abs(p)
            min1 = 3.0 * xm * q - abs(tol1 * q)
            min2 = abs(e * q)
            if 2.0*p < min([min1, min2]): 
                # Accept interpolation.
                e, d = d, p/q
            else: # Interpolation failed, use bisection.
                d = xm
                e = d
        else: # Bounds decreasing too slowly, use bisection.
            d = xm
            e = d
        # Move last best guess to a.
        a, fa = b, fb
        # Evaluate new trial root.
        if abs(d) > tol1:
            b += d
        else:
            b += sign(tol1, xm)
        fb = f(b)
    raise ValueError("Method failed, iterations exceed the max loop number. "+\
                     "Try other method instead or use a different interval.")

def newton_raphson(f,df,a,b,epsilon=TOL,LoopMax=NMAX):
    r"""
    Newton-Raphson method with its recurrence relation 
                         f(x[i])            f(x[i])
        x[i+1] = x[i] - --------- = x[i] - ----------.
                         f'(x[i])           df(x[i])
    where f is the function and df is its first order derivative.
    
    Reference:
        [1] William et al., Numerical Recipes in C The Art of Scientific 
            Computing, 2nd Edition, 2002: 362-366.     
        
    """
    rt = (a + b)/2.0  # Initial guess.
    for N in range(LoopMax):
        fx = f(rt)
        dfx = df(rt)
        dx = fx/dfx
        if (a - rt)*(rt - b) < 0.0:
            raise ValueError("Jumped out of brackets in Newton-Raphson method.")
        
        if abs(dx) <= TOL:
            print("Newton-Raphson method:%d iterations." %N)
            return rt
    raise ValueError("Method failed, iterations exceed the max loop number. "+\
                     "Try other method instead or use a different interval.")   

newton = newton_raphson

def bis_nowton(f,df,a,b,epsilon=TOL,LoopMax=NMAX):
    r"""
    A combination of bisection and NewtonRaphson. The hybrid algorithm takes a
    bisection step whenever Newton-Raphson would take the solution out of 
    bounds, or whenever Newton-Raphson is not reducing the size of the 
    brackets rapidly enough.
    
    Reference:
        [1] William et al., Numerical Recipes in C The Art of Scientific 
            Computing, 2nd Edition, 2002: 366-367.  
    """
    fl, fh = f(a), f(b)
    if fl*fh > 0.0:
        raise ValueError("Root must be bracketed in bisection-Nowton-Raphaon method.")
    if fl == 0.0: return a
    if fh == 0.0: return b
    # Orient the search so that f(xl) < 0.
    if fl < 0.0:
        xl, xh = a, b
    else:
        xh, xl = a, b
    rt = (a + b)/2.0    # Initialize the guess for root
    dxold = abs(b - a)  # the "stepsize before last"
    dx = dxold          # the last step
    frt, dfrt = f(rt), df(rt)
    for N in range(LoopMax):
        ## Bisect if Nowton out of range or not decreasing fast enough.
        if (((rt - xh)*dfrt-  frt) * ((rt - xl) * dfrt - frt) > 0.0) or \
            (abs(2.0 * frt) > abs(dxold * dfrt)):
            dxold = dx
            dx = (xh - xl)/2.0
            rt = xl + dx
            # Change in root is negligible.
            if (xl == rt):
                print("Bisection-Nowton-Raphson method:%d iterations."%N)
                return rt
        else: ## Nowton step acceptable. Take it.
            dxold = dx
            dx = frt / dfrt
            temp = rt
            rt -= dx
            if temp == rt:
                print("Bisection-Nowton-Raphson method:%d iterations."%N)
                return rt
        if abs(dx) < epsilon: # Convergence criterion.
            print("Bisection-Nowton-Raphson method:%d iterations."%N)
            return rt
        frt, dfrt = f(rt), df(rt)
        # Maintain the bracket on the root.
        if frt < 0.0:
            xl = rt
        else:
            xh = rt
    raise ValueError("Method failed, iterations exceed the max loop number. "+\
                     "Try other method instead or use a different interval.")   

if __name__ == "__main__":
    def f1(x):
        return cos(x) - x**3
    def f1_(x):
        return sin(x) - 3*x**2

    print(bisection(f1,0,1,5e-15,100))
    print(falsi(f1,0,1,5e-15,100))
    print(illinois(f1,0,1,5e-15,100))
    print(secant(f1,0,1,5e-15,100))
    print(ridders(f1,0,1,5e-15,100))
    print(brent(f1,0,1,5e-15,100))
    #print(newton(f1,f1_,0,31,5e-15,100))
    print(bis_nowton(f1,f1_,0,31,5e-15,100))