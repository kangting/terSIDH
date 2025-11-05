"""
x-only Montgomery isogenies

Even torsion algorithms from: https://ia.cr/2017/1198
Computing Isogenies between Montgomery Curves Using the Action of (0, 0)
Joost Renes

Odd torsion algorithms: https://ia.cr/2017/504.pdf
A simple and compact algorithm for SIDH with arbitrary degree isogenies
Craig Costello and Huseyin Hisil

Codomain computation for velu formula from: https://ia.cr/2018/782
A faster way to the CSIDH
Michael Meyer and Steffen Reith

VeluSqrt from https://velusqrt.isogeny.org/
Faster computation of isogenies of large prime degree
Daniel J. Bernstein, Luca De Feo, Antonin Leroux, Benjamin Smith

TODO: allow composition by defining __mul__ on isogenies to create a composite isogeny
"""

# Sage imports
from sage.all import prod, ZZ, PolynomialRing
# ProductTree lives in two places 
from sage.misc.banner import require_version
if require_version(9,8):
      from sage.rings.generic import ProductTree
else:
      from sage.schemes.elliptic_curves.hom_velusqrt import ProductTree

# Local imports
from montgomery_xz import KummerLine, KummerPoint
import time
import random

def evaluate_factored_kummer_isogeny(phi_list, all_list, P, dummy=True):
    """
    Given a list of isogenies, evaluates the
    point for each isogeny in the list
    """
    if dummy:
        for phi in all_list:
            Q = phi(P)
            if phi in phi_list:          
                P = Q
    else:
        for phi in all_list:
            P = phi(P)
            
    return P
    

def factored_kummer_isogeny(K, oP, oorder, nP, norder, cP, corder, primes, threshold=1000):
    """
    """
    
    def sparse_isogeny_prime_power(P, l, e, split=0.8, threshold=1000):
        """
        Compute chain of isogenies quotienting
        out a point P of order l**e
        https://trac.sagemath.org/ticket/34239
        """

        if l > threshold:
            KummerLineIsogenyAlgorithm = KummerLineIsogeny_VeluSqrt
        else:
            KummerLineIsogenyAlgorithm = KummerLineIsogeny

        def recursive_sparse_isogeny(Q, k):
            assert k
            if k == 1:  # base case
                return [KummerLineIsogenyAlgorithm(Q.parent(), Q, l, check=False)]

            k1 = int(k * split + 0.5)
            k1 = max(1, min(k - 1, k1))  # clamp to [1, k-1]

            Q1 = l**k1 * Q
            L = recursive_sparse_isogeny(Q1, k - k1)

            Q2 = evaluate_factored_kummer_isogeny(L, L, Q, dummy=False)
            R = recursive_sparse_isogeny(Q2, k1)
            return L + R
        
        return recursive_sparse_isogeny(P, e)

    # For computing points 
    cofactor1 = oorder
    cofactor2 = norder    
    cofactor3 = corder
        
    # Ensure P is a point on E

    if oP.parent() != K:
        raise ValueError(f"The supplied kernel must be a point on the line {K}")
    if nP.parent() != K:
        raise ValueError(f"The supplied kernel must be a point on the line {K}")
    if cP.parent() != K:
        raise ValueError(f"The supplied kernel must be a point on the line {K}")
        
    psi_list = []
    phi_list = []
    phi_constant_list = []
    all_list = []
    orders = [oorder, norder, corder]
    random.shuffle(orders)
    
    factors = []
    first1 = True
    first2 = True
    first3 = True
    tmp = 1
    # Deal with isomorphisms
    for order in orders:
        if order == 1:
            raise ValueError("TODO: deal with isomorphisms")  
                
        factors = list(order.factor())
        random.shuffle(factors)
        
        for l, e in factors:
            
            # Map P through chain length e of l-isogenies
            psi_list = []

            D = ZZ(l**e)
            primes //= D
            
            if order == oorder:
                cofactor1 //= D

                # Compute point Q of order l^e
                Q = primes * oP
                
                # Use Q as kernel of degree l^e isogeny
                psi_list = sparse_isogeny_prime_power(Q, l, e, threshold=threshold)        
                phi_list += psi_list
                
                oP = evaluate_factored_kummer_isogeny(psi_list, psi_list, oP, dummy = False)
                nP = evaluate_factored_kummer_isogeny(psi_list, psi_list, nP, dummy = False)
                cP = evaluate_factored_kummer_isogeny(psi_list, psi_list, cP, dummy = False)

            elif order == norder:
                cofactor2 //= D

                # Compute point Q of order l^e
                Q = primes * nP
                
                # Use Q as kernel of degree l^e isogeny
                psi_list = sparse_isogeny_prime_power(Q, l, e, threshold=threshold)    
                phi_list += psi_list
                
                oP = evaluate_factored_kummer_isogeny(psi_list, psi_list, oP, dummy = False)
                nP = evaluate_factored_kummer_isogeny(psi_list, psi_list, nP, dummy = False)
                cP = evaluate_factored_kummer_isogeny(psi_list, psi_list, cP, dummy = False)
                
            else:         
                tmp = 0
                cofactor3 //= D

                # Compute point Q of order l^e 
                Q = primes * cP
                
                # Use Q as kernel of degree l^e isogeny
                psi_list = sparse_isogeny_prime_power(Q, l, e, threshold=threshold)
                phi_constant_list += psi_list    
                
                cP = evaluate_factored_kummer_isogeny(psi_list, psi_list, cP, dummy = False)
                cP = evaluate_factored_kummer_isogeny(psi_list, psi_list, cP, dummy = False)
                cP = evaluate_factored_kummer_isogeny(psi_list, psi_list, cP, dummy = False)

            all_list += psi_list
        
        if cofactor1 == 1 and first1 == True:
            oP *= oorder
            nP *= oorder
            cP *= oorder
            first1 = False
        if cofactor2 == 1 and first2 == True:
            oP *= norder
            nP *= norder
            cP *= norder
            first2 = False
        if cofactor3 == 1 and first3 == True:
            oP *= corder
            nP *= corder
            cP *= corder
            first3 = False

    
    return phi_list, all_list


class KummerLineIsogeny_Generic:
    """
    TODO
    Generic class for x-only isogenies, I'm bad at classes
    so this is probably bad...
    """
    def __init__(self):
        pass

    def __repr__(self):
        return f"Isogeny of degree {(self._degree).factor()} from {self._domain} to {self._codomain}"

    @staticmethod
    def validate_input(domain, kernel, degree, check=True):        

        if not isinstance(domain, KummerLine):
            raise ValueError(f'not a kummer line: {domain}')
    
        if not isinstance(kernel, KummerPoint):
            raise ValueError(f'not a kummer point: {kernel}')
        
        if kernel.parent() != domain:
            raise ValueError(f'Kernel {kernel} is not a point on {domain}')
        
        if check:
            # TODO actually check order with has_order_D function
            assert (degree*kernel).is_zero(), "Input point does not have correct order"


    def domain(self):
        """
        """
        return self._domain
    
    def codomain(self):
        """
        """
        return self._codomain

    def degree(self):
        """
        """
        return self._degree
        
    
class KummerLineIsogenyComposite(KummerLineIsogeny_Generic):
    """
    Computes composite degree isogenies as a chain of prime
    degree isogenies
    """
    
    
    def __init__(self, domain, oB, oP, order, nB, nP, norder, cB, cP, corder, primes, check=True, threshold=1000):
        # Check the input to the isogeny is well-formed
        # self.validate_input(domain, oB*oP, order, check=check)
        # self.validate_input(domain, nB*nP, norder, check=check)
        # self.validate_input(domain, cB*cP, corder, check=check)

        # Compute factored isogeny
        self._phis, self._all_phis = factored_kummer_isogeny(domain, oP, order, nP, norder, cP, corder, primes, threshold=threshold)

        # Make immutable
        self._phis = tuple(self._phis)
        self._all_phis = tuple(self._all_phis)

        # Compute degree, domain and codomain
        self._degree = prod(phi.degree() for phi in self._phis)
        self._domain = self._phis[0].domain()
        self._codomain = self._phis[-1].codomain()

    def __call__(self, P):
        
        """
        """        
        
        return evaluate_factored_kummer_isogeny(self._phis, self._all_phis, P, dummy=True)

    @classmethod
    def from_factors(cls, maps):
        maps = tuple(maps)

        L = maps[0].domain()
        for phi in maps:
            if not isinstance(phi, KummerLineIsogeny_Generic):
                raise TypeError(f'not an kummer-line isogeny: {phi}')
            if phi.domain() != L:
                raise ValueError(f'isogeny has incorrect domain: {phi}')
            L = phi.codomain()

        result = cls.__new__(cls)

        # Make immutable
        result._phis = maps

        # Compute degree, domain and codomain
        result._degree = prod(phi.degree() for phi in result._phis)
        result._domain = result._phis[0].domain()
        result._codomain = result._phis[-1].codomain()
        
        return result

class KummerLineIsogeny(KummerLineIsogeny_Generic):
    """
    Computes prime degree isogenies with Renes/Costello-Hisil formula
    """
    def __init__(self, domain, kernel, degree, check=True):
        # Check the input to the isogeny is well-formed
        self.validate_input(domain, kernel, degree, check=check)

        # Set kernel and degree and domain
        self._degree = degree
        self._kernel = kernel
        self._domain = domain

        # Compute the codomain, we need different formula for even and 
        # odd degree
        if self._degree == 2:
            # We cannot use the point (0 : 0 : 1) as the kernel for
            # these formula
            assert self._kernel.XZ()[0]

        # Compute the codomain
        self._codomain = self._compute_codomain()

    def __call__(self, P):
        if not isinstance(P, KummerPoint):
            raise ValueError
        if self._degree == 2:
            return self._evaluate_isogeny_even(P)
        return self._evaluate_isogeny(P)
    
    def _precompute_edwards_multiples(self, d):
        """
        """
        # Compute the [i]K for i in [1...d]
        K_muls = self._kernel.multiples()
        E_muls = []
        for _ in range(d):
            Ki = next(K_muls)
            KX, KZ = Ki.XZ()
            YE = KX - KZ
            ZE = KX + KZ
            E_muls.append((YE, ZE))
        return E_muls
        
    def _compute_codomain_constants(self):
        """
        """
        # Extract Montgomery constants
        A, C = self._domain.extract_constants()

        # Compute and store pairs of points for later evaluation
        d = (self._degree - 1) // 2
        self._edwards_multiples = self._precompute_edwards_multiples(d)

        # Convert to twisted Edwards curve parameters (Aed,Ded)  
        Ded = C + C
        Aed = A + Ded
        Ded = A - Ded

        prod_Y = 1
        prod_Z = 1
        for EY, EZ in self._edwards_multiples:
            prod_Y *= EY
            prod_Z *= EZ

        # compute prod_Y^8 and prod_Z^8
        prod_Y, prod_Z = prod_Y**2, prod_Z**2
        prod_Y, prod_Z = prod_Y**2, prod_Z**2
        prod_Y, prod_Z = prod_Y**2, prod_Z**2

        # A_new = A_old^ell * prod_Z^8
        # D_new = D_old^ell * prod_Y^8
        Aed = Aed**self._degree * prod_Z
        Ded = Ded**self._degree * prod_Y

        # Change back to Montgomery-parameters
        A = Aed + Ded
        C = Aed - Ded           
        A = A + A
                
        return A, C
    
    def _compute_codomain_constants_even(self):
        """
        """
        # Extract kernel point
        XK, ZK = self._kernel.XZ()
        assert XK, "XK Cannot be zero"

        # C = ZK^2
        C = ZK * ZK
        
        # A = 2*(ZK^2 - 2*XK^2)
        A = XK * XK         # A = XK^2
        A = A + A           # A = 2*XK^2
        A = C - A           # A = ZK^2 - 2XK^2
        A = A + A           # A = 2*(ZK^2 - 2*XK^2)
        return A, C
            
    def _compute_codomain(self):
        # Compute the codomain constants, need different formula for 
        # odd and even ell
        if self._degree == 2:
            A_codomain, C_codomain = self._compute_codomain_constants_even()
        else:
            A_codomain, C_codomain = self._compute_codomain_constants()

        # Constuct a new KummerLine
        F = self._domain.base_ring()
        return KummerLine(F, [A_codomain, C_codomain])

    def _evaluate_isogeny(self, P):
        """
        """
        XP, ZP = P.XZ()
        Psum = XP + ZP
        Pdiff = XP - ZP

        # Loop through the d-multiples
        X_new, Z_new = 1, 1
        for EY, EZ in self._edwards_multiples:
            diff_EZ = Pdiff * EZ
            sum_EY  = EY * Psum
            X_new *= (diff_EZ + sum_EY)
            Z_new *= (diff_EZ - sum_EY)

        # Square and multiple with original
        X_new = X_new**2 * XP
        Z_new = Z_new**2 * ZP
        
        return self._codomain((X_new, Z_new))
        
    def _evaluate_isogeny_even(self, P):
        """
        """
        XK, ZK = self._kernel.XZ()
        assert XK, "XK cannot be zero"

        XP, ZP = P.XZ()

        # TODO: make this friendlier looking?
        T0 = XK + ZK
        T1 = XK - ZK
        T2 = XP + ZP
        T3 = ZP - XP # Typo in formula: paper says XP - ZP
        T4 = T3 * T0 # (ZP - XP)(XK + ZK)
        T5 = T2 * T1 # (XP + ZP)(XK - ZK)
        T6 = T4 - T5 # (ZP - XP)(XK + ZK) - (XP + ZP)(XK - ZK)
        T7 = T4 + T5 # (ZP - XP)(XK + ZK) + (XP + ZP)(XK - ZK)
        T8 = XP * T6 # XP * ((ZP - XP)(XK + ZK) - (XP + ZP)(XK - ZK))
        T9 = ZP * T7 # ZP * ((ZP - XP)(XK + ZK) + (XP + ZP)(XK - ZK))

        return self._codomain((T8, T9))


def product_tree_resultant(hI_tree, poly):
    r"""
    Helper function to evaluate a resultant with `h_I` quickly,
    using the product tree, taken from FastEllipticPolynomial
    sage/src/sage/schemes/elliptic_curves/hom_velusqrt.py

    Author: Lorenz Panny (2022)
    """
    rems = hI_tree.remainders(poly)
    r = prod(rems)
    s = -1 if len(hI_tree)%2 == 1 == poly.degree() else 1
    assert r.is_constant()
    return s * r[0]


class KummerLineIsogeny_VeluSqrt(KummerLineIsogeny_Generic):
    def __init__(self, domain, kernel, degree, check=True):
        # Check the input to the isogeny is well-formed
        self.validate_input(domain, kernel, degree, check=check)

        # Set kernel and degree and domain
        self._degree = degree
        self._kernel = kernel
        self._domain = domain

        # We need the domain coefficient for the elliptic
        # resultants.
        # TODO: can we make this projective? 
        self.a = self._domain.a()

        # We need a polynomial ring, so we create it once
        # and store it to self
        k = self._domain.base_ring()
        self.R = PolynomialRing(k, names="Z")
        self.Z = self.R.gen()

        # baby step and giant step params
        b = (self._degree - 1).isqrt() // 2
        c = (self._degree - 1) // (4*b)
        self.stop = self._degree-4*b*c

        # Pre-compute polynomials which are needed
        # throughout. hI is stored as a product tree
        # for faster resultants
        self.hI_tree  = self._hI_precomputation(kernel, b, c)
        self.EJ_parts = self._EJ_precomputation(kernel, b)
        self.hK = self._hK_precomputation(kernel, degree, b, c)

        # Compute the codomain
        self._codomain = self._compute_codomain()

    def __call__(self, P):
        if not isinstance(P, KummerPoint):
            raise ValueError
        return self._evaluate_isogeny(P)

    def _hI_resultant(self, poly):
        """
        """
        return product_tree_resultant(self.hI_tree, poly)

    def _hI_precomputation(self, ker, b, c):
        """
        """
        Q = (2*b)*ker
        step, diff = Q.double(), Q
        leaves = []
        for i in range(c):
            leaves.append(self.Z - Q.x())
            if i < c - 1:
                Q, diff = Q.add(step, diff), Q
        
        return ProductTree(leaves)

    def _Fs(self, X1, X2):
        """
        """
        X1X2 = X1 * X2
        polys = (
                (X1 - X2)**2,
                -2*( (X1X2 + 1)*(X1 + X2) + 2*self.a*X1X2),
                (X1X2 - 1)**2
            )
        return polys

    # def _Fs_new(self, QX, QZ):
    #     """
    #     """
    #     # TODO: using projective like this seems slower?
    #     ZQX = self.Z*QX
    #     ZQZ = self.Z*QZ
    #     polys = (
    #             (ZQZ - QX)**2,
    #             -2*( (ZQX + QZ)*(QX + ZQZ) + 2*self.a*ZQX*QZ),
    #             (ZQX - QZ)**2
    #         )
    #     return polys
    
    def _EJ_precomputation(self, ker, b):
        """
        """
        Q = ker
        step, diff = Q.double(), Q
        EJ_parts = []
        for i in range(b):
            polys = self._Fs(self.Z, Q.x())
            EJ_parts.append(polys)
            if i < b - 1:
                Q, diff = Q.add(step, diff), Q

        return EJ_parts

    def _hK_precomputation(self, ker, ell, b, c):
        """
        """
        hK = []
        Q = ker.double()
        step, next_point = Q, Q.double()
        stop = ell-4*b*c
        for i in range(2, stop, 2):
            QX, QZ = Q.XZ()
            hK.append(QZ * self.Z - QX)
            if i < stop - 1:
                Q, next_point = next_point, next_point.add(step, Q)

        return self.R(prod(hK))
    
    # def _hK_codomain(self):
    #     """
    #     """
    #     hK_0, hk_1 = 1, 1

    #     # Loop through kernel
    #     Q = self._kernel.double()
    #     step, next_point = Q, Q.double()
        
    #     for i in range(2, self.stop, 2):
    #         QX, QZ = Q.XZ()

    #         hK_0 *= ( QZ - QX)
    #         hk_1 *= (-QZ - QX)

    #         if i < self.stop - 1:
    #             Q, next_point = next_point, next_point.add(step, Q)

    #     return hK_0, hk_1

    # def _hK_evaluate(self, alpha, alpha_inv):
    #     """
    #     """
    #     hK_0, hK_1 = 1, 1

    #     Q = self._kernel.double()
    #     step, next_point = Q, Q.double()
    #     for i in range(2, self.stop, 2):
    #         QX, QZ = Q.XZ()
            
    #         hK_0 *= (QZ * alpha_inv - QX)
    #         hK_1 *= (QZ * alpha - QX)

    #         if i < self.stop - 1:
    #             Q, next_point = next_point, next_point.add(step, Q)

    #     return hK_0, hK_1
        
    def _compute_codomain_constants(self):
        """
        """
        # TODO: should this be projectivised?
        # These are the polynomials for alpha = 1 and alpha = -1
        E0J = prod(F0 + F1 + F2 for F0,F1,F2 in self.EJ_parts)
        E1J = prod(F0 - F1 + F2 for F0,F1,F2 in self.EJ_parts)

        # Compute resultants and evaluate hK at 1 and -1 
        R0  = self._hI_resultant(E0J)
        R1  = self._hI_resultant(E1J)
        M0  = self.hK( 1)
        M1  = self.hK(-1)
        # M0, M1 = self._hK_codomain()

        # We have that
        # d = [(A - 2C)(A + 2C)]^ell * (hS(1) / hS(-1))^8
        # hS = hK R / Delta
        # d = [(A - 2C)(A + 2C)]^ell * (hK R(1) / hK R(-1))^8

        # First compute (hS(1) / hS(-1))^8
        num = R0 * M0
        den = R1 * M1
        # num^8, den^8 with three squares
        num, den = num**2, den**2
        num, den = num**2, den**2
        num, den = num**2, den**2

        # [(A - 2)(A + 2)]^ell
        num = (self.a - 2)**self._degree * num
        den = (self.a + 2)**self._degree * den

        # Compute the new curve y^2 = x^3 + (A:C)x^2 + x
        A_new = num + den
        C_new = den - num      # C =   d - n
        A_new = A_new + A_new  # A = 2(n + d)

        return A_new, C_new
    
    def _compute_codomain(self):
        """
        """
        A_codomain, C_codomain = self._compute_codomain_constants()
        F = self._domain.base_ring()
        return KummerLine(F, [A_codomain, C_codomain])  

    def _evaluate_isogeny(self, P):
        """
        """
        if P.is_zero():
            return self._codomain((1, 0))

        # x-coordinate of point to evaluate
        alpha = P.x()
        alpha_inv = 1/alpha

        # TODO: can I make this projective without needing alpha?
        # Compute two polynomials from giant steps
        EJ0 = prod((F0 * alpha_inv + F1) * alpha_inv + F2 for F0,F1,F2 in self.EJ_parts)
        EJ1 = prod((F0 * alpha + F1) * alpha + F2 for F0,F1,F2 in self.EJ_parts)

        # Resultants and evaluations 
        R0 = self._hI_resultant(EJ0)
        R1 = self._hI_resultant(EJ1)
        M0  = self.hK(alpha_inv)
        M1  = self.hK(alpha)
        # M0, M1 = self._hK_evaluate(alpha, alpha_inv)

        # Make new point
        X_new = (R0 * M0)**2 * alpha**self._degree
        Z_new = (R1 * M1)**2

        return self._codomain((X_new, Z_new))