# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 12:31:46 2017

@author: Feng Li
"""

from math import sqrt,exp,e,pi,log,cosh,tan,tanh,atan,erf
import roots
#from time import process_time

DEFAULT_PDA = 1e-10 # default precision degree of accuracy(pda)

# Euler–Mascheroni constant referenced from 
# https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant
EM_GAMMA = 0.57721566490153286060651209008240243104215933593992

## --------------------------------------------------------
## Fast way to calculate factorial is referred to
## http://blog.csdn.net/gxqcn/article/details/18810
def factorial(n, m=1):
    '''
    Factorial of an non-negtive integer n.
    fac(n) = n!
           = fac2(n)*fac2(n-1)
    n≥m≥1
    fac(n,m) = n*(n-1)*...*m =
             2^((n-m)//2) * fac((n-1)//2,(m+1)//2) * fac2(n,m), odd n & odd m
             2^((n-m+1)//2) * fac((n-1)//2,m//2) * fac2(n,m+1), odd n & even m
             2^((n-m+1)//2) * fac(n//2,(m+1)//2) * fac2(n-1,m), even n & odd m
             2^((n-m+2)//2) * fac(n//2,m//2) * fac2(n-1,m+1), even n & even m
    when n is 0 then fac(0,m) is set to 0 whatever m is.
    '''
    try:
        if not isinstance(n, int):
            n = int(n)
        if n < 0:
            raise ValueError("number n must be an non-negtive integer.")
    except:
        raise ValueError("Unknown parameter n = "+str(n)+" given.")
        
    if n == 0: return 1
    
    try:
        if not isinstance(m, int):
            m = int(m)
        if m < 0 or m > n:
            print("m must be an non-negtive integer and m ≤ n.")
            raise ValueError
    except:
        print(n, m)
        raise ValueError("Unknown parameter m = "+str(m)+" of type "+str(type(m))+"  given.")
        
    if m == 0: m = 1    
    if m == n: return n

    if n % 2 == 1:
        if m % 2 == 1: # n奇m奇
            return 2**((n-m)//2)*fac((n-1)//2, (m+1)//2)*fac2(n, m)
        else:         # n奇m偶
            return 2**((n-m+1)//2)*fac((n-1)//2, m//2)*fac2(n, m+1)
    else:
        if m % 2 == 1: # n偶m奇
            return 2**((n-m+1)//2)*fac(n//2, (m+1)//2)*fac2(n-1, m)
        else:          # n偶m偶
            return 2**((n-m+2)//2)*fac(n//2, m//2)*fac2(n-1, m+1)

fac = factorial
## ----------------------------------------
def fac2(n,m=1):
    '''
    Factorial of an non-negtive integer n with an interval of 2.
    fac2(n) = n!! =
              n*(n-2)*...*2 = 2^(n//2) * fac(n//2), even n
              n*(n-2)*...*1,                        odd n
    fac2(n,m) =
            for odd n:
              n*(n-2)*...*(m+1), even m
              n*(n-2)*...*m    , odd m
            for even n: 
              n*(n-2)*...*m = 2^((n-m+2)//2)*fac(n//2,m//2), even m (m>=2)
              n*(n-2)*...*(m+1) = 2^((n-m+1)//2)*fac(n//2,(m+1)//2), odd m            
    when n is 0 then fac(0,m) is set to 0 whatever m is.
    '''
    #print(n,m)
    try:
        if not isinstance(n, int):
            n = int(n)
        if n < 0:
            raise ValueError("number n must be an non-negtive integer.")
    except:
        raise ValueError("Unknown parameter n = "+str(n)+" given.")
    if n == 0: return 1
    
    try:
        if not isinstance(m, int):
            m = int(m)
        if m < 0 or m > n:
            raise ValueError("m must be an non-negtive integer and m≤n.")
    except:
        raise ValueError("Unknown parameter m = "+str(m)+" given.")
    
    if m == n: return n
    if m == 0: m = 1
    
    result = 1
    if n % 2 == 1: # n为奇数
        # 使 end 的奇偶性与 n 保持一致
        if m % 2 == 1:
            end = m
        else:
            end = m + 1
        for i in range(n, end-2, -2):
            result *= i
        return result
    else: # n 为偶数
        if m % 2 == 1:
            return 2**((n-m+1)//2)*fac(n//2, (m+1)//2)
        else:
            return 2**((n-m+2)//2)*fac(n//2, m//2)

def arrangement(n,m):
    '''
    Calculate the number of arrangements to get m articles in a line from n articles.
    A(n,m) = n!/(n-m)! = fac(n)/fac(n-m) = fac(n,n-m+1)
    '''
    if m == 0: 
        return 1
    else:
        return fac(n, n-m+1)

A = arrangement


def combination(n,m):
    '''
    Calculate the number of combinations to get m articles from n articles.
    C(n,m) = A(n,n)/A(m,m) = n!/[m!*(n-m)!] = fac(n)/[fac(m)*fac(n-m)]
           = fac(n,n-m+1)//fac(m)
    '''
    return A(n,m)//fac(m)

C = combination
## ----------------------------------------
def gammaf(alpha):#,step=1e-5, epsilon=1e-5):
    """
    Gamma function with parameter α > 0, where it is
           +∞ 
    Γ(α) = ∫ x^(α-1)*exp(-x) dx
           0 
    Γ(1) = 1, Γ(0.5) = sqrt(pi)
    Γ(α+1) = αΓ(α) for all α > 0
    Γ(α) = (α-1)Γ(α-1) for all α > 1
    Thus if α = α0 + n, where α0 is in range (0,1), then
    Γ(α) = Γ(α0 + n) = (α-1)Γ(α-1) # recursive expression
                     = (α-1)(α-2)...(α-(n-1))α0Γ(α0) # function expression
    So, what we need is just to calculate the cases where 0<α<1.
    Γ(n+1) = nΓ(n) = n!
    """
    try:
        if not isinstance(alpha, int) or not isinstance(alpha, float):
            alpha = float(alpha)
        if alpha <= 0:
            raise ValueError("α must be positive number.")
    except:
        raise ValueError("Unknown parameter α = "+str(alpha)+"is given.")
        
    if alpha < 1:
        """
        The following analysis method is much much more efficient that the 
        integration method and has a better accuracy!
        Referred to http://www.sosmath.com/calculus/improper/gamma/gamma.html
                                +∞
                           1   ─┬─┬─   e^(x/n)
        Γ(x) = e^(-C*x) * --- * | | -------------
                           x    │ │    1 + x/n
                                n=1
        where the constant
                     1     1           1 
        C = lim 1 + --- + --- + ... + --- - ln(n)
            n→+∞     2     3           n
        is the Euler–Mascheroni constant γ whose value is 
        0.57721566490153286060651209008240243104215933593992…
        Referred to 
        https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant
        """
#        n = 2
#        const = []
#        C0 = 1
#        const.append(C0)
#        while 1:
#            C1 = C0 + 1/n
#            C2 = C1 - log(n)
#            const.append(C2)
#            if abs(const[n-1] - const[n-2]) > DEFAULT_PDA: 
#                # Using DEFAULT_PDA garenteed its high accuracy.
#                C0 = C1
#                n += 1
#            else:
#                C = const[-1]
#                print("n = %d, C = "%n,end='')
#                print(C)
#                break
        C = EM_GAMMA
        n = 1
        result = exp(-C*alpha)/alpha
        while 1:
            term = exp(alpha/n)/(1 + alpha/n)
            if abs(term - 1) > DEFAULT_PDA:
                result *= term
                n += 1
            else:
                #print("nMax = %d"%n)
                return result
    elif alpha == 1.0:
        return 1
    else:
        if isinstance(alpha, int):
            return fac(alpha - 1)
        if isinstance(alpha, float):
            return (alpha - 1)*gammaf(alpha - 1) 

## ----------------------------------------
def digamma(x):
    r"""
    Digamma function, it is defined as the logarithmic derivative of the 
    gamma function:
                d               Γ'(x)
        ψ(x) = ---- ln[Γ(x)] = -------
                dx              Γ(x)
            
    It has the following properties:
        ψ(x + 1) = ψ(x) + 1/x      ... Eq(1)
    thus,
        ψ(x) = ψ(x + 1) - 1/x      ... Eq(1')
    For x > 6.0, the following approximation is used
        ψ(x) = ln(x) - (x^-1)/2 - (x^-2)/12 + (x^-4)/120 - (x^-6)/252 
               + (x^-8)/240 - (x^-10)*5/660 + (x^-12)*691/32760 - (x^-14)/12
               + O(x^-16)
    this formula yields "more than enough precision" (at least 12 digits 
    except near the zeroes), while for x ≤ 6.0, use Eq(1') repeatedly to 
    calculate.
    
    Reference:
        [1] https://en.wikipedia.org/wiki/Digamma_function
        [2] https://en.wikipedia.org/wiki/Digamma_function#Computation_and_approximation
    """
    try:
        if not isinstance(x, float):
            x = float(x)
        if x <= 0.0:
            raise ValueError
    except:
        errmsg = \
        '''
        Independent variable x in digamma function must be positive real number.
        '''
        raise ValueError(errmsg)
    
    if x > 6.0:
        return log(x) - 0.5/x - x**(-2)/12.0 + x**(-4)/120.0 - x**(-6)/252 +\
                x**(-8)/240.0 - 5.0*x**(-10)/660 + 691.0/32760*x**(-12) -\
                x**(-14)/12.0
    else:
        return digamma(x + 1.0) - 1.0/x

psif = digamma
## ----------------------------------------
def Betaf(a, b):
    '''
    Beta function with parameter a > 0 and b > 0, where it is
                  1                              Γ(a)Γ(b)
    Betaf(a, b) = ∫ x^(a-1)*(1 - x)^(b-1) dx = ------------
                  0                              Γ(a + b)
    '''
    try:
        if not isinstance(a, float):
            a = float(a)
        if a <= 0:
            raise ValueError("parameter a must be greater than 0.")
    except:
        raise ValueError("Unknown parameter a = "+str(a)+"is given.")
    
    try:
        if not isinstance(a, float):
            a = float(a)
        if a <= 0:
            raise ValueError("parameter a must be greater than 0.")
    except:
        raise ValueError("Unknown parameter a = "+str(a)+"is given.")
        
    return gammaf(a)*gammaf(b)/gammaf(a+b)
    
def sech(x):
    return 1.0/cosh(x)
    
    
## --------------------------------------------------------

def E(X,k=1):
    """
    k th moment of random variable X, which is equal to 
    X.E(), by default it is equal to X.EX.
    """
    return X.E(k)

def Var(X):
    r"""
    Return the variance of X.
    """
    return X.var()
    
def F(X,x):
    """
    Cumulative distribution function F of random variable X at X ≤ x, 
    which is equal to X.dis_function(x) or X.F(x)
    F(x) = P(X ≤ x)
    """
    return X.F(x)

def Q(X, p):
    """
    Quantile function (the inverse cumulative density function (cdf) ).
    """
    return X.Q(p)
    
def Pr(X, x, xlow=None, RangeType=0):
    """
    Probability of random variable X at X=x.
    Support three type of output according to parameter RangeType
    0: Pr(x) = P(X=x)
    1: Pr(x) = P(X≤x)
    2: Pr(x) = P(xlow≤X≤x)
    """
    return X.Pr(x, xlow, RangeType)
    
#def Pr(X,xupper,xlower=None):
#    """
#    Probability of random variable X at X=x.
#    Pr(X,x) = X.P(X=x)
#    If 
#    """
#    if xlower == None:
#        return X.Pr(xupper)
#    else:
#        try:
#            if not isinstance(xlower, float):
#                xlower = float(xlower)
#            if xlower >= xupper:
#                raise ValueError("Lower bound should be less than upper lower.")
#            return X.Pr(xupper)-X.Pr(xlower)
#        except:
#            raise ValueError("Unknower lower bound is given.")

    
## -------------- Error handling functions ----------------
def numCheck(x, Type, flag):
    r"""
    """
    if Type == int: 
        typeName = "integer"
    elif Type == float: 
        typeName = "real number"
        
    if flag == 1:
        try:
            if not isinstance(x, Type):
                x = Type(x)           
            if x < 0:
                raise ValueError("Parameter must be non-negative "+ typeName +".")
        except:
            raise ValueError("Unknown parameter " + str(x) + "is given.")
    elif flag == 2:
        try:
            if not isinstance(x, Type):
                x = Type(x)        
            if x <= 0 or x >= 1:
                raise ValueError("Parameter must be float in range (0,1).")
        except:
            raise ValueError("Unknown parameter " + str(x) + "is given.")
    return x
## --------------------------------------------------------

############################################################################
####             Random Variant (Discrete and continuous)               ####
#### ================================================================== ####
class RandomVariant():
    r"""
    Random variables which is either discrete or continuous. Other types of 
    random variables are currently not supportable.
    """
    def __init__(self, VarType=None):
        self.VarType = VarType
        
        if VarType == "discrete":
#            self.dis_serList = []
#            self.nMax = 0
            #self.defaultRangeType = 0
            self.EX = self.expectation() # 期望值
            self.VarX = self.var() # 方差
            self.sigma = sqrt(self.VarX) # 标准差
            
            self.F = self.cdf # 分布函数 cumulative distribution function,F(x)
            self.E = self.moment # k阶原点矩函数（默认时为期望函数）
            self.Var = self.var # 方差函数[值](对外提供一个统一的操作接口)
        elif VarType == "continuous":
            #self.defaultRangeType = 1
            self.EX = self.expectation() # 期望值
            self.VarX = self.var() # 方差
            self.sigma = sqrt(self.VarX) # 标准差
            self.F_x05 = self._median() # 中位数 median
            self.median = self._median() # 提供对外统一的接口
            self.H = self.entropy() # Shannon entropy
            
            self.F = self.cdf # 分布函数 cumulative distribution function,F(x)
            self.E = self.moment # k阶原点矩函数（默认时为期望函数）
            self.Var = self.var # 方差函数[值] (对外提供一个统一的操作接口)
            self.Q = self.quantile = self._quantile # p分位数函数 p quantile
        else:
            raise ValueError("Unsupported type of random variable.\nSupport"/
                             +" only discrete and continuous random variables.")
            
        if self.EX != 0:
            self.CvX = self.sigma/self.EX # 变异系数 
        if self.VarX != 0:
            self.betaS = self.central_moment(3)/self.VarX**1.5 # 偏度系数 skewness
            self.betaK = self.central_moment(4)/self.VarX**2 - 3 # 峰度系数
        
    def pdf(self, x):
        r"""
        Probability density function, Pr(x) = P(X=x) = p(x).
        Need to be overloaded when enherited.
        """
        pass
    
    def Pr(self, x, xlow=None, RangeType=0):
        """
        Probability of random variable X
        Support three type of output according to parameter RangeType
        0: Pr(x) = P(X=x)
        1: Pr(x) = P(X≤x)
        2: Pr(x) = P(xlow≤X≤x)
        """
        # 概率密度函数 probability density function, Pr(x)
        try:
#            if self.VarType == "discrete":
#                if not isinstance(x, int):
#                    x = int(x)
#            else:
#                if not isinstance(x, float):
#                    x = float(x)
            if not isinstance(x, float):
                    x = float(x)
        except:
            raise ValueError("Unknown parameter "+str(x)+" is given.")
        
        rangeMsg = '''\t0: Pr(x) = P(X=x)
                      \t1: Pr(x) = P(X≤x)
                      \t2: Pr(x) = P(xlow≤X≤x)
                   '''        
        if RangeType == 0:
            # Pr(x) = P(X=x)
            if self.VarType == "discrete":
                return self.pdf(x)
            elif self.VarType == "continuous":
                return 0.0
            else:
                raise ValueError(self.VarType+\
                                 " variable is currently unsupported.")                
        elif RangeType == 1:
            # Pr(x) = P(X≤x)
            return self.cdf(x)
        elif RangeType == 2:
            # Pr(x) = P(xlow≤X≤x)
            try:
                if not isinstance(xlow, float):
                    xlow = float(xlow)
            except:
                raise ValueError("Low bound parameter "+str(xlow)+" should be a real number.")
            
            if x >= xlow:
                return self.cdf(x)- self.cdf(xlow)
            else:
                return self.cdf(xlow) - self.cdf(x)
        else:
            raise ValueError("Unknown range type is given.\n"+rangeMsg)
            
    def cdf(self, x, x0=0.0, n=1000):
        r"""
        Cumulative distribution function, F(x) = P(X≤x).
        For discrete variable, enherited without overloading.
        For continuous one, enherite directly or overload (for a better
        accuracy acquired by analytical results) it as needed.
        When X is a continuous random variant, 
        using Simpson's rule to calculate the integral 
               x
        F(x) = ∫ f(t) dt
               x0
                h
             = ---[f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + ... + 4f(x[n-1]) + f(x[n])]
                3
        where n must be even and h = x/n.
        Its accuracy is |E| ≤ M(b-a)*h^4/180, where M = max{|f''''(x)|}
                                                        [a,b]
        For better accurracy, it can be overloaded if there is an analytical
        form for its cdf.
        """
        try:
            if not isinstance(x, float):
                x = float(x)
        except:
            raise ValueError("Unknown parameter "+str(x)+" is given.")
            
        if self.VarType == "discrete":        
            Sum = 0.0
            for i in range(self.nMax + 1):
                if i <= x:
                    Sum += self.dis_serList[i]
                else:
                    break
            return Sum
            
        elif self.VarType == "continuous":
            # Check: n must be even integer and x must be greater than 0.            
            n = n // 2 * 2
            h = x / n
            
            xi = []
            xi.append(x0)
            for i in range(1,n):
                xi.append(xi[i-1] + h)
            xi.append(x)
            
            Sum = 0.0
            for i in range(n+1):
                if i == 0 or i == n:
                    Sum += self.pdf(xi[i])
                elif i % 2 == 1:
                    Sum += 4*self.pdf(xi[i])
                else:
                    Sum += 2*self.pdf(xi[i])
            
            integral = Sum * h / 3
            return integral

    def _quantile(self, p=0.5):
        r"""
        p quantile of continuous Random variant X, where F(xp) = P(X≤xp) = p.
        By default, p = 0.5 which is the median of X.
        For discrete enherited class, it cannot be called directly.
        """
        print("Method quantile: Not implemented.")
        return None
    
    def _median(self): 
        return self._quantile()
        
    def expectation(self):
        return self.moment()
     
    def moment(self, k=1):
        """
        Calculate the k th moment of the random variable X, that is
                       +∞
    	   μ[k] = E(X^k) = ∫ x^k*P(X=x)], k∈N+.
                       -∞
        When the random variable X is discrete, it's enherited directly,
        while it need to be overloaded when X is continuous.
    	  """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Parameter k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k given in " + str(k))
            
        if k == 0:
            return 1
        else:
            if self.VarType == "discrete":
                # For discrete X, enherited 
                Sum = 0.0
                for i in range(self.Count+1):
                    Sum += i**k * self.dis_serList[i]
                return Sum
            elif self.VarType == "continuous":
                # Need to be overloaded when enherinted.
                print("Not implemented.")
                pass
    
    def var(self):
        return self.central_moment()
    
    def central_moment(self, k=2):        
        r"""
        Calculate the k th central moment of the random variable X, that is
                        +∞             k
        E((X-E(X))^k) = ∫x^k*P(X=x)] = Σ C(k,i) * μ[k] * (- EX)^(k-i)
                        -∞            i=0
        By default, k = 2 which gives the variance of X.
        """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Parameter k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k = " + str(k) + "is given.")
        #k = numCheck(k,int,1)
        Sum = 0.0
        for i in range(k+1):
            Sum += C(k,i) * self.moment(i) * (-self.EX)**(k-i)
        return Sum

    def entropy(self):
        r"""
        Shannon entropy (information entropy) is defined as the average amount 
        of information produced by a probabilistic stochastic source of data.
        """
        print("Entropy: Not implemeted.")
        pass
        
#### ====================== Discrete random variable X ==================####
class Binomial(RandomVariant):
    """
    Binomial distribution with parameters n and p, where its distribution series
    P(X=i) = C(n,i)*p^i*(1-p)^(n-i), i = 0,1,2,3,...,n
    
    Reference:
        [1] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, n, p):
        try:
            if not isinstance(n, int):
                n = int(n)
            if n < 0:
                raise ValueError("n must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter n given in "+str(n))
        try:
            if not isinstance(p, float):
                p = float(p)
            if p < 0 or p > 1:
                raise ValueError("Value error: 0≤p≤1")
        except:
            raise ValueError("Unknown parameter p given in "+str(p))
            
        self.n = n
        self.p = p
        self.dis_serList = []
        for i in range(self.n + 1):
            self.dis_serList.append(self.pdf(i))
        self.nMax = self.n
        self.Count = self.nMax
        
        super().__init__(VarType="discrete")
        
    def pdf(self, i):
        r"""
        Probability density function, Pr(x) = P(X=x) = p(x).
        """
        return C(self.n,i)*(self.p**i)*((1-self.p)**(self.n-i))
    
b = Binomial

## ---------------------------------------
class Poisson(RandomVariant): 
    """
    Poisson distribution with parameters n and p, where its distribution series
    P(X=i) = λ^i * exp(-λ)/i!, i = 0,1,2,...,+∞ [nMax]
    with default precision degree of accuracy(pda) of 1e-10 for distribution 
    series calculations.
    
    Reference:
        [1] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, lamd, pda=DEFAULT_PDA):
        try:
            if not isinstance(lamd, float):
                lamd = float(lamd)
            if lamd < 0:
                raise ValueError("Value error:λ must be greater than 0.")
        except:
            raise ValueError("Unknown parameter λ given in " + str(lamd))

        self.lamd = lamd
        self.pda = pda # 指定的计算精度
        # 计算其分布列至指定精度pda
        self.dis_serList = []
        i = 0
        while 1:
            prob = self.dis_series(i)
            if prob < self.pda:
                self.nMax = i - 1
                break
            else:
                self.dis_serList.append(prob)
                i += 1
        
        self.Count = self.nMax        
        super().__init__(VarType="discrete")

    def pdf(self, i):
        r"""
        Probability density function, Pr(x) = P(X=x) = p(x).
        """
        return (self.lamd**i)*exp(-self.lamd)/fac(i)

P = Poisson

## ---------------------------------------
class Geometric(RandomVariant):
    """
    Geometric distribution with parameters p, where its distribution series
    P(X=i) = (1-p)^(i-1) * p, i = 1,2,3...,+∞ [nMax in fact]
    with default precision degree of accuracy(pda) of 1e-10 for distribution 
    series calculations.
    
    Reference:
        [1] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, p, pda=DEFAULT_PDA):
        try:
            if not isinstance(p, float):
                p = float(p)
            if p < 0 or p > 1:
                raise ValueError("Value error:0≤p≤1.")
        except:
            raise ValueError("Unknown parameter p given in " + str(p))

        self.p = p            
        self.pda = pda
        # 分布列 P(X=i), i=1,2,3,...,nMax
        self.dis_serList = []
        i = 1
        while 1:
            prob = self.dis_series(i)
            if prob < self.pda:
                self.nMax = i - 1
                break
            else:
                self.dis_serList.append(prob)
                i += 1
                
        self.Count = self.nMax
        super().__init__(VarType="discrete")

    def pdf(self, i):
        r"""
        Probability density function, Pr(x) = P(X=x) = p(x).
        """
        return (1 - self.p)**(i-1) * self.p

Ge = Geometric

## ---------------------------------------
class Hypergeometric(RandomVariant):
    """
    Hypergeometric distribution with parameters (n,N,M), 
    where its distribution series is
    P(X=i) = C(M,i) * C(N-M,n-i) / C(N,n), i = 0,1,2,...,r
    r = min{M,n}
    """
    def __init__(self, n, N, M):
        try:
            if not isinstance(n, int):
                n = int(n)
            if not isinstance(N, int):
                N = int(N)
            if not isinstance(M, int):
                M = int(M)
            if N < M or N < n or N < 0 or M <0 or n < 0:
                raise ValueError("n,N,M must be positive integers and " +\
                                 "M<=N, n<=N.")
        except:
            raise ValueError("n,N,M must be positive integers.")

        self.n = n
        self.N = N
        self.M = M
        self.r = min(M,n)
        self.nMax = self.r
        #计算分布列
        self.dis_serList = [0]*(self.r+1)
        for i in range(self.r+1):
            self.dis_serList[i] = self.dis_series(i)
        
        self.Count = self.nMax
        super().__init__(VarType="discrete")

    def pdf(self, i):
        r"""
        Probability density function, Pr(x) = P(X=x) = p(x).
        """
        return C(self.M,i)*C(self.N-self.M,self.n-i)/C(self.N,self.n)

h = Hypergeometric

## ---------------------------------------
class Negtive_binomial(RandomVariant):
    """
    Negtive binomial distribution with parameters r and p, 
    where its distribution density function is 
    P(X=i) = C(i-1,r-1)*p^r*(1-p)^(i-r), i = r,r+1,...,+∞ [nMax in fact]
    r = 1,2,3,...
    with default precision degree of accuracy(pda) of 1e-10 for distribution 
    series calculations.
    """
    def __init__(self, r, p, pda=DEFAULT_PDA):
        try:
            if not isinstance(r, int):
                r = int(r)
            if r < 1:    
                raise ValueError("r must be a positive integer.")
        except:
            raise ValueError("Unknown parameter r given in "+str(r))
        try:
            if not isinstance(p, float):
                p = float(p)
            if p < 0 or p > 1:
                raise ValueError("Value error: 0≤p≤1")
        except:
            raise ValueError("Unknown parameter p given in "+str(p))
            
        self.r = r
        self.p = p
        self.pda = pda
        # 分布列 P(X=i)
        self.dis_serList = []
        i = r
        while 1:
            dis = self.dis_series(i)
            if dis < self.pda:
                self.nMax = i - 1
                break
            else:
                self.dis_serList.append(dis)
                i += 1
                
        self.Count = self.nMax -self.r
        super().__init__(VarType="discrete")

    def pdf(self, i):
        r"""
        Probability density function, Pr(x) = P(X=x) = p(x).
        """
        return C(i-1,self.r-1)*(self.p**self.r)*((1-self.p)**(i-self.r))

Nb = Negtive_binomial    

#### ==================== Continuous random variable X ==================####
class Uniform(RandomVariant):
    r"""
    Uniform distribution of X,where its 
    distribution series is
                 ┌ 1 / (b - a), a < x < b
        P(X=x) = ┤
                 └ 0          , other x
    distribution function is
                        ┌       0       , x ≤ a
        F(x) = P(X≤x) = ┤(x - a)/(b - a), a < x < b
                        └       1       , x ≥ b
    with default precision degree of accuracy(pda) of 1e-10 for distribution 
    series calculations.
    
    Reference:
        [1] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, a, b, pda=DEFAULT_PDA):
        try:
            if not isinstance(a, float):
                err = a
                a = float(a)
            if not isinstance(b, float):
                err = b
                b = float(b)
            if b < a:
                raise ValueError("b must be greater than a.")
        except:
            raise ValueError("Unknown parameter "+str(err)+" given.")
        
        self.a, self.b = a, b
        super().__init__(VarType="continuous")
    
    def pdf(self, x):
        r"""
        Probability density function, Pr(x) = P(X=x) = p(x) = 1 / (b - a).
        """
        if x < self.b and x > self.a:
            return 1/(self.b - self.a)
        else:
            return 0.0
        
    def cdf(self, x):
        r"""
        The analytical form of cdf for Uniform distribution is
        F(x) = P(X≤x) = (x - a)/(b - a)
        """
        return (x - self.a) / (self.b - self.a)

    def _quantile(self, p=0.5):
        r"""
        p quantile of continuous Random variant X, where F(xp) = P(X≤xp) = p.
        That is (xp - a)/(b - a) = p, thus
        xp = p*b + (1-p)*a
        """
        return p * self.b + self.a * (1 - p)
        
    def moment(self, k=1):
        """
        Calculate the k th moment of the random variable X, that is
                        b
    	  μ[k] = E(X^k) = ∫x^k*P(X=x)] = [b^(k+1)-a^(k+1)]/[(k+1)*(b-a)]
                        a
    	  """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Power k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k given in " + str(k))
            
        result = (self.b**(k+1)-self.a**(k+1))/((k+1)*(self.b-self.a))
        return result

U = Uniform

## ---------------------------------------
class stdN(RandomVariant):
    """
    Standard normal distribution of X with μ=0 and σ=1, where its 
    distribution density function is 
        P(X=x) = exp(-x^2/2)/sqrt(2π), -∞< x <+∞
    distribution function is
                        x
        F(x) = P(X≤x) = ∫ exp(-t^2/2)/sqrt(2π) dt = Φ(x)
                        -∞
    Reference:
        [1] https://en.wikipedia.org/wiki/Normal_distribution
        [2] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self):
        self.mu = 0.0
        self.sigma2 = 1.0
        self.sigma = 1.0

        super().__init__(VarType="continuous")

    def pdf(self, x):
        return exp(-x**2/2)/sqrt(2*pi)
    
    def cdf(self, x, n=1000):
        # Enherited from RandomVariant and change the lower bound x0 to 0.0
        # RandomVariant only give the right half of the integral of pdf
        # Thus, a value of 0.5 must be plused.
        # The enherited Simpson's method is okey, the erf method can also be
        # used for calculations for the math module has build-in erf method.
        #return super().cdf(x,0.0,n) + 0.5 # Accuarate and fast enough.
        '''
        Using the relative error function to calculate the cdf of stdN
        Φ(x) = 0.5*[1 + erf(x/sqrt(2))]
        For normal distribution with parameters μ and σ, its cdf is
        F(x) = 0.5*{1 + erf((x - μ)/(sqrt(2)*σ))}
        
        Reference:
            https://en.wikipedia.org/wiki/Normal_distribution#Cumulative distribution function
        '''
        return 0.5 + 0.5 * erf((x-self.mu)/(sqrt(2)*self.sigma))
    
    def _quantile(self, p=0.5):
        r"""
        p quantile of continuous Random variant X, where F(xp) = P(X≤xp) = p.
        It is acturally the inverse function of F(x) for F(x) is monotonic 
        function. By default, p = 0.5 which is the median of X.
        -------------- Algorithm --------------
        Using Peter J. Acklam's lower tail quantile for standard normal 
        distribution function at the website bellow
        https://web.archive.org/web/20070505093933/http://home.online.no/~pjacklam/notes/invnorm/
        and the python source code at
        https://web.archive.org/web/20061005013356/http://home.online.no/~pjacklam/notes/invnorm/impl/field/ltqnorm.txt
        About the algorithm:
            It's said by the author that 
            'he algorithm uses a minimax approximation by rational functions
            and the result has a relative error whose absolute value is less
            than 1.15e-9.'
        Here, the code is copied directly (only changed the function name).
        """
        if p <= 0 or p >= 1:
            # The original perl code exits here, we'll throw an exception instead
            raise ValueError( "Argument to ltqnorm %f must be in open interval (0,1)" % p )

        # Coefficients in rational approximations.
        a = (-3.969683028665376e+01,  2.209460984245205e+02, \
             -2.759285104469687e+02,  1.383577518672690e+02, \
             -3.066479806614716e+01,  2.506628277459239e+00)
        b = (-5.447609879822406e+01,  1.615858368580409e+02, \
             -1.556989798598866e+02,  6.680131188771972e+01, \
             -1.328068155288572e+01 )
        c = (-7.784894002430293e-03, -3.223964580411365e-01, \
             -2.400758277161838e+00, -2.549732539343734e+00, \
              4.374664141464968e+00,  2.938163982698783e+00)
        d = ( 7.784695709041462e-03,  3.224671290700398e-01, \
              2.445134137142996e+00,  3.754408661907416e+00)
    
        # Define break-points.
        plow  = 0.02425
        phigh = 1 - plow
    
        # Rational approximation for lower region:
        if p < plow:
           q  = sqrt(-2*log(p))
           return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                   ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    
        # Rational approximation for upper region:
        if phigh < p:
           q  = sqrt(-2*log(1-p))
           return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    
        # Rational approximation for central region:
        q = p - 0.5
        r = q*q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    
    def moment(self, k=1):
        """
        Calculate the k th moment of the random variable X, that is
                       +∞
    	  μ[k] = E(X^k) = ∫x^k*P(X=x)]
                       -∞
                       ┌
                     = |  0      , k 为奇数
                       |  (k-1)!!, k 为偶数
                       └
    	  """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Value error: k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k given in " + str(k))
          
        if k == 0:
            return 1
        elif k % 2 == 1:
            return 0
        else:
            return fac2(k-1)
    
    def entropy(self):
        r"""
            H(X) = ln(2πe σ^2)/2
        """
        return log(2.0 * pi * e * self.sigma2)/2.0
## ---------------------------------------
class Normal(RandomVariant):
    """
    Normal distribution of X,where its pdf is 
        P(X=x) = exp[-(x - μ)^2 / (2*σ^2)]/[sqrt(2π)*σ], -∞< x <+∞
    and its cdf is
                        x
        F(x) = P(X≤x) = ∫ exp[- (t - μ)^2/(2σ^2)] dt / [sqrt(2π)*σ]
                        -∞
    with default precision degree of accuracy(pda) of 1e-10 for distribution 
    series calculations.
    
    Reference:
        [1] https://en.wikipedia.org/wiki/Normal_distribution
        [2] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, mu, sigma2):   
        self.mu = mu
        self.sigma2 = sigma2
        self.sigma = sqrt(sigma2)
        
        super().__init__(VarType="continuous")
        
    def pdf(self, x):
        return exp(-(x - self.EX)**2/(2*self.sigma2))/(sqrt(2*pi)*self.sigma)
    
    def cdf(self, x, n=1000):
        U = stdN()
        u = (x - self.mu) / self.sigma
        return U.cdf(u, n)
        
    def _quantile(self, p=0.5):
        r"""
        p quantile of continuous Random variant X, where F(xp) = P(X≤xp) = p.
        By default, p = 0.5 which is the median of X.
        """
        return self.mu + self.sigma * stdN().quantile(p)
        
    def moment(self, k=1): 
        r"""
        Calculate the k th moment of the random variable X, that is
                       +∞
    	  μ[k] = E(X^k) = ∫x^k*P(X=x)]
                       -∞
               k
             = Σ C（k,i）*σ^i*μ^(k-i)* E(U^i)
              i=0
        where
            U = (X - μ)/σ,                       
                     ┌
            E(U^i) = |    0    , i 为奇数
                     |  (i-1)!!, i 为偶数
                     └
    	  """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Value error: k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k given in " + str(k))
            
        if k == 0:
            return 1
        else:
            Sum = 0.0
            U = stdN()
            for i in range(k+1):
                Sum += C(k,i) * self.sigma**i * self.mu**(k-i) * E(U,i)
            return Sum
    
    def entropy(self):
        r"""
            H(X) = ln(2πe σ^2)/2
        """
        return log(2.0 * pi * e * self.sigma2)/2.0  
  
N = Normal
Gaussian = Normal
## ---------------------------------------
class LogNormal(Normal):
    r"""
    Logarithmic normal distribution of X, its pdf is
             ┌ exp[-(ln x - μ)^2 / (2*σ^2)]/[x*sqrt(2π)*σ], x > 0
    P(X=x) = ┤
             └                   0                        , x ≤ 0
    
    Reference:
        [1] https://en.wikipedia.org/wiki/Log-normal_distribution
        [2] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, mu, sigma2):
        super().__init__(mu, sigma2)
        
    def pdf(self, x):
        if x >  0.0:
            return exp(-(log(x) - self.EX)**2/(2*self.sigma2))/(x*sqrt(2*pi)*self.sigma)
        else:
            return 0.0
        
    def cdf(self, x, n=1000):
        r"""
        X ~ LN(μ,σ^2), then U = (ln X - μ) / σ ~ N(0,1), X = exp{μ + σU}
        F(x) = P(X≤x) = P(exp(μ+σU) ≤ x) = P(U ≤ (ln x - μ)/σ) = Φ((ln x - μ)/σ)
        """
        if x <= 0.0:
            return 0.0
        else:
            U = stdN()
            u = (log(x) - self.mu)/self.sigma
            return U.cdf(u, n)
        
    def _quantile(self, p=0.5):
        r"""
        p quantile of continuous Random variant X, where F(xp) = P(X≤xp) = p.
        By transformation, U = (ln X - μ)/σ ~ N(0, 1), thus
        p = F(xp) = P(X≤xp) = P(u ≤ exp(μ+σ*xp)) = Φ(exp(μ+σ*xp)) = Φ(up)
        xp = (ln up - μ)/σ
        where up = stdN()._quantile(p).
        By default, p = 0.5 which is the median of X.
        """
        return (log(stdN()._quantile(p)) - self.mu) / self.sigma
        
    def moment(self, k=1): 
        """
        Calculate the k th moment of the random variable X, that is
                       +∞
    	  μ[k] = E(X^k) = ∫x^k*P(X=x)]
                       -∞
             = exp(kμ + k^2*σ^2/2)
    	  """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Value error: k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k given in " + str(k))
            
        return exp(k*self.mu + k**2 * self.sigma2 / 2)
    
    def entropy(self):
        r"""
            H(X) = ln(sqrt(2π)*exp(μ+0.5）*σ)
        """
        return log(sqrt(2.0 * pi) * exp(self.mu + 0.5) * self.sigma)
        
LN = LogNormal

## ---------------------------------------
class Gamma(RandomVariant):
    """
    Gamma distribution of X,where its distribution series is 
     
             ┌   λ^α * x^(α-1) * exp(-λx) 
             | ---------------------------    ,x ≥ 0
    P(X=x) = ┤            Γ(α) 
             |
             └             0                  , x < 0
    and its distribution function is
             ┌   λ^α  x
             | ------ ∫ t^(α-1) * exp(-λt) dt , x ≥ 0
    F(x)   = ┤  Γ(α) -∞
             |
             └              0                 , x < 0
    where shape parameter α > 0 and scale parameter λ > 0.
    
    Reference:
        [1] https://en.wikipedia.org/wiki/Gamma_distribution
        [2] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, alpha, lamd): #, pda=DEFAULT_PDA
        try:
            if not isinstance(alpha, float):
                alpha = float(alpha)
            if alpha <= 0:
                raise ValueError("α must be positive number.")
        except:
            raise ValueError("Unknown parameter α = "+str(alpha)+"is given.")
    
        try:
            if not isinstance(lamd, float):
                lamd = float(lamd)
            if lamd <= 0:
                raise ValueError("λ must be positive number.")
        except:
            raise ValueError("Unknown parameter λ = "+str(lamd)+"is given.")
        
        self.alpha, self.lamd = alpha, lamd
        self._gamma_alpha = gammaf(self.alpha) # To save time.

        super().__init__(VarType="continuous")
        
    def pdf(self, x):
        if x <= 0:
            return 0.0
        return (self.lamd**self.alpha) * (x**(self.alpha-1)) * \
                (exp(-self.lamd*x)) / self._gamma_alpha #gammaf(self.alpha)
    
    def _quantile(self, p=0.5): 
        r"""
        p quantile of continuous Random variant X, where F(xp) = P(X≤xp) = p.
        By default, p = 0.5 which is the median of X.
        """
        # Using root-finding algorithms to get the quantile xp.
        def _cdf_p(x):
            return self.cdf(x) - p
        return roots.secant(_cdf_p, 0.0, 1.0, 1e-15, 100)
        #return roots.illinois(_cdf_p, 0.0, 1.0, 1e-15, 100)
        
    def moment(self, k=1):
        """
        Calculate the k th moment of the random variable X, that is
                       +∞
    	   μ[k] = E(X^k) = ∫ x^k*P(X=x)
                       0
                       Γ(α+k)
                     = ------ * λ^(-k) = α(α+1)...(α+k+1)*λ^(-k)
                        Γ(α)
    	  """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Value error: k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k given in " + str(k))
            
        if k == 0:
            return 1
        else:
            return self.lamd**(-k) * gammaf(self.alpha + k) / self._gamma_alpha #gammaf(self.alpha)
    
    def entropy(self):
        r"""
            H(X) = α - ln(λ) + ln[Γ(α)] - (1 - α)ψ(α)
        where ψ(x) is the digamma function.
        """
        return self.alpha - log(self.lamd) + log(gammaf(self.alpha)) -\
                (1.0 - self.alpha)*psif(self.alpha)

Ga = Gamma
## ---------------------------------------
class Exponential(Gamma):
    """
    Exponential distribution of X,where its distribution density function is 
    P(X=x) = λ * exp(-λx), x ≥ 0
    and it is the Gamma distribution when α = 1, that is
    Exp(λ) = Ga(1, λ).
    
    Reference:
        [1] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, lamd): # , pda=DEFAULT_PDA
        super().__init__(1.0, lamd) #, pda)
    
Exp = Exponential

## ---------------------------------------
class Chi2(Gamma):
    """
    Chi-square distribution of X, where its distribution density function is
                    1
    P(X=x) = --------------- x^(n/2 - 1) * exp(-x/2), x ≥ 0
              2^(n/2)Γ(n/2)
    where n is a positive real number
    and it is the Gamma distribution when α = n/2, λ = 1/2, that is
    Chi2(n) = Ga(n/2, 1/2).
    
    Reference:
        [1] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, n): # , pda=DEFAULT_PDA
        super().__init__(n/2, 0.5) # ,pda)

## ---------------------------------------
class BetaD(RandomVariant):
    """
    Beta distribution of X,where its distribution series is 
             ┌ x^(a-1) * (1-x)^(b-1) * Γ(a+b)/[Γ(a)*Γ(b)], 0<x<1
    P(X=x) = ┥
             └                0          , other
    and its distribution function is
    
    with default precision degree of accuracy(pda) of 1e-10 for distribution 
    series calculations.
    
    Reference:
        [1] 茆诗松,程依明,濮晓龙.概率论与数理统计教程（第二版）,高等教育出版社：2011.2
    """
    def __init__(self, a, b): 
        self.a, self.b = a, b
        
        super().__init__(VarType="continuous")
        
    def pdf(self, x):
        try:
            if not isinstance(x, (float, int)):
                x = float(x)
        except:
            raise ValueError("Unknown parameter x = "+str(x)+" is given.")
            
        if x > 0 and x < 1:
            return x**(self.a-1) * (1-x)**(self.b-1) *\
                   gammaf(self.a+self.b) / (gammaf(self.a)*gammaf(self.b))
        else:
            return 0.0
        
    def _quantile(self, p=0.5):
        # Using root-finding algorithms to get the quantile xp.
        def _cdf_p(x):
            return self.cdf(x) - p
        return roots.secant(_cdf_p, 0.0, 1.0, 1e-15, 100)
            
    def moment(self, k=1):
        """
        Calculate the k th moment of the random variable X, that is
                       +∞
    	  μ[k] = E(X^k) = ∫x^k*P(X=x)]
                       0
                         Γ(a+b)Γ(a+k)        (a+k-1)(a+k-2)...a
                      = -------------- = ----------------------------, k∈N+.
                         Γ(a)Γ(a+b+k)     (a+b+k-1)(a+b+k-2)...(a+b)
    	  """
        try:
            if not isinstance(k, int):
                k = int(k)                
            if k < 0:
                raise ValueError("Value error: k must be non-negtive integer.")
        except:
            raise ValueError("Unknown parameter k given in " + str(k))
            
        if k == 0:
            return 1
        else:
            return gammaf(self.a + self.b) * gammaf(self.a + k) \
                    / gammaf(self.a) / gammaf(self.a + self.b + k)

Be = BetaD

## ---------------------------------------     
class Logistic(RandomVariant):
    r"""
    Logistic distribution of X. Its pdf is
                          exp[-(x - μ)/s]          1          x - μ
    p(x) = P(X=x) = --------------------------- = --- sech^2(-------)
                     s*{1 + exp[-(x - μ)/s]}^2     4s           2s
    x∈R, location parameter μ∈R, scale parameter s∈R and s > 0.
    
    Reference:
        https://en.wikipedia.org/wiki/Logistic_distribution
    """
    def __init__(self, mu, s):
        # Check if s > 0.
        if not isinstance(s, float):
            s = float(s)
        if s <= 0.0:
            raise ValueError("Parameter s must be positive real number.")
        
        self.mu, self.s = mu, s
        
        self.median = mu
        self.VarX = (self.s * pi)*(self.s * pi) / 3
        self.betaS = 0
        self.betaK = 1.2
        
        super().__init__(VarType="continuous")
        
    def pdf(self, x):
        try:
            if not isinstance(x, float):
                x = float(x)
        except:
            raise ValueError("Unknown parameter x = "+str(x)+" is given.")
        
        return (sech((x - self.mu)/(2*self.s)))**2/(4*self.s)
        
    def cdf(self, x):
        """
        Cumulative distribution function
                        1             1     1        x-μ
        F(x) = ------------------- = --- + --- tanh(-----)
                1 + exp[-(x-μ)/s]     2     2         2s
        """
        return 0.5 + 0.5*tanh((x - self.mu)/(2 * self.s))
        
    def _quantile(self, p=0.5):
        """
        The quantile function (the inverse cdf) is
        Q(p) = μ + s*ln[p/(1-p)], 0 < p < 1.
        """
        try:
            if not isinstance(p, float):
                p = float(p)
            if p <= 0 or p >= 1:
                raise ValueError
        except:
            raise ValueError("p must be real number in the open interval (0,1).")
        return self.mu + self.s * log(p / (1 - p))
    
    def moment(self, k=1):
        pass
    
    def entropy(self):
        r"""
            H(X) = ln(s) + 2
        """
        return log(self.s) + 2.0
## ---------------------------------------
class Cauchy(RandomVariant):
    r"""
    Cauchy distribution  of X. Its pdf is
                                1
        p(x;μ,λ) = --------------------------- ,
                       ┌        x - μ   ^2 ┐
                    πλ*| 1 + （---------)   |
                       └          λ        ┘
    and its cdf is
                    1            x - μ      1
        F(x;μ,λ) = --- * arctan(-------) + --- ,
                    π              λ        2
    where x∈R, location parameter μ∈R and scale parameter λ>0.
    
    Reference:
        https://en.wikipedia.org/wiki/Cauchy_distribution
    """
    def __init__(self, mu, lamd):
        try:
            if not isinstance(mu, float):
                mu = float(mu)
        except:
            raise ValueError("Location parameter μ must be real number.")
    
        try:
            if not isinstance(lamd, float):
                lamd = float(lamd)
            if lamd <= 0:
                raise ValueError
        except:
            raise ValueError("Scale parameter λ must be positive real number.")
        
        self.mu, self.lamd = mu, lamd
        
        super().__init__(VarType="continuous")
        self.VarX, self.betaK, self.betaS = [None]*3
    
    def pdf(self, x):
        return 1.0/pi/self.lamd/(1.0 + ((x - self.mu)/self.lamd)**2)
    
    def cdf(self, x):
        return atan((x - self.mu)/self.lamd)/pi + 0.5

    def _quantile(self, p=0.5):
        r"""
        Quantile (inverse cdf) of Cauchy distribution variant X, which is
            Q(x;μ,λ) = μ + λ*tan[π(p - 0.5)]
        in analytical form.
        """
        return self.mu + self.lamd * tan(pi * (p - 0.5))

    def entropy(self):
        r"""
            H(X) = ln(4πλ).
        """
        return log(4.0 * pi * self.lamd)
    
Cau = Cauchy
## ---------------------------------------
class Weibull(RandomVariant):
    r"""
    Weibull distribution of X. Its pdf is
    
                   ┌  α      x                  x
                   | --- * (---)^(α-1) * exp[-(---)^α] , x ≥ 0
        p(x;α,λ) = ┤  λ      λ                  λ
                   |
                   └                   0               , x < 0
    and its cdf is
                   ┌            x
                   | 1 - exp[-(---)^α] , x ≥ 0
        F(x;α,λ) = ┤            λ
                   |
                   └         0         , x < 0
    where x ≥ 0, shape parameter α ≥ 0 and scale parameter λ ≥ 0.
                   
    """
    def __init__(self, alpha, lamd):
        try:
            if not isinstance(alpha, float):
                alpha = float(alpha)
            if alpha <= 0.0:
                raise ValueError
        except:
            raise ValueError("Shape parameter α must be positive real number.")
    
        try:
            if not isinstance(lamd, float):
                lamd = float(lamd)
            if lamd <= 0:
                raise ValueError
        except:
            raise ValueError("Scale parameter λ must be positive real number.")
        
        self.alpha, self.lamd = alpha, lamd
        
        super().__init__(VarType="continuous")
        
    def pdf(self, x):
        if x < 0.0:
            return 0.0
        else:
            return self.alpha/self.lamd * (x/self.lamd)**(self.alpha-1) * \
                    exp(-(x/self.lamd)**self.alpha)

    def cdf(self, x):
        if x < 0.0:
            return 0.0
        else:
            return 1 - exp(-(x / self.lamd)**self.alpha)
    
    def _quantile(self, p=0.5):
        r"""
        Q(p;α,λ) = λ[- ln(1 - p)]^(1/α) for 0≤ p < 1.
        """
        try:
            if not isinstance(p, float):
                p = float(p)
            if p >= 1.0 or p < 0.0:
                raise ValueError
        except:
            raise ValueError("Probability value must lie in [0,1).")
        
        return self.lamd * (- log(1.0 - p)) * (1.0 / self.alpha)

    def moment(self, k=1):
        r"""
                       +∞
        μ[k] = E(X^k) = ∫ x^k*P(X=x)] = λ^k * Γ(1 + k/α), k∈N+.
                        0
        """
        return self.lamd**k * gammaf(1.0 + k/self.alpha)

    def entropy(self):
        r"""
            H(X) = γ(1-1/α) + ln(λ/α) + 1
        where γ is Euler–Mascheroni constant.
        """
        return EM_GAMMA * (1.0 - 1.0/self.alpha) +\
               log(self.lamd/self.alpha) + 1.0
## ---------------------------------------








