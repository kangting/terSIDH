"""
Montgomery curve Kummer Line
"""
from sage.all import cached_method, Integer, EllipticCurve

from sage.structure.element import RingElement
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field

class KummerLine:
    def __init__(self, *args):
        r"""
        """
        self._curve = None

        # Allow the creation of the Kummer Line from an EllipticCurve
        if len(args) == 1:
            curve, = args
            if not isinstance(curve, EllipticCurve_generic):
                raise TypeError('not an elliptic curve')
            ainvs = curve.a_invariants()
            A, C = ainvs[1], 1
            if ainvs != (0, A, 0, 1, 0):
                raise TypeError('Must use Montgomery model')
            self._curve = curve
            self._base_ring = curve.base_ring()

        # Allow the creation of the Kummer line from a base field and coeffs.
        elif len(args) == 2:
            base_ring, curve_constants = args
            # Extract curve constants
            if isinstance(curve_constants, Integer) or len(curve_constants) == 1:
                A = curve_constants
                C = 1
            elif len(curve_constants) == 2:
                A, C = curve_constants
            else:
                raise ValueError("TODO")
            self._base_ring = base_ring
        else:
            raise ValueError("TODO")

        # init variables
        self._A = self._base_ring(A)
        self._C = self._base_ring(C)

        # Make sure A,B are not bad
        if (self._A**2 - 4*self._C**2).is_zero():
            raise ValueError(f"Constants {curve_constants} do not define a Montgomery curve")
    
    def __eq__(self, other):
        if self.base_ring() != other.base_ring():
            return False
        return self._A * other._C == other._A * self._C

    def __repr__(self):
        r"""
        """
        if self.a():
            return f'Kummer line of the Montgomery curve y^2 = x^3 + {self.a()._coeff_repr()}*x^2 + x over {self.base_ring()}'
        else:
            return f'Kummer line of the Montgomery curve y^2 = x^3 + x over {self.base_ring()}'
        
    def __call__(self, coords):
        return KummerPoint(self, coords)
        
    def base_ring(self):
        r"""
        """
        return self._base_ring

    def extract_constants(self):
        r"""
        """
        return self._A, self._C
    
    def zero(self):
        return self(None)
    
    def curve(self):
        if not self._curve:
            self._curve = self.montgomery_curve()
        return self._curve
    
    @cached_method
    def montgomery_curve(self):
        F = self.base_ring()
        a = self.a()
        return EllipticCurve(F, [0,a,0,1,0])
    
    @cached_method
    def short_weierstrass_curve(self):
        F = self.base_ring()
        A = self.a()

        A_sqr  = A*A
        A_cube = A*A_sqr
        a = 1 - A_sqr/3
        b = (2*A_cube - 9*A) / 27
        return EllipticCurve(F, [a,b])
    
    @cached_method
    def j_invariant(self):
        """
        """
        j_num = 256 * (self._A**2 - 3*self._C**2)**3
        j_den = self._C**4 * (self._A**2 - 4*self._C**2)
        return j_num / j_den

    @cached_method
    def a(self):
        """
        """
        return self._A / self._C

    def is_isomorphic(self, other):
        if not isinstance(other, KummerLine):
            raise TypeError("TODO")
        return self.j_invariant() == other.j_invariant()


class KummerPoint:
    def __init__(self, parent, coords):
        r"""
        """
        if not isinstance(parent, KummerLine):
            raise TypeError('not a Montgomery Kummer line')

        R = parent.base_ring()

        # Point at infinity
        if coords is None:
            coords = (Integer(1), Integer(0))
        # Construct point from P on an elliptic curve in Montgomery form
        elif isinstance(coords, EllipticCurvePoint_field):
            # Make sure point's parent curve matches with Kummer Line
            a = parent.a()
            assert coords.curve().a_invariants() == (0,a,0,1,0)
            coords = coords[0], coords[2]
        # Construct from X coordinate only
        elif isinstance(coords, RingElement):
            coords = coords,
        # Construct from a tuple (X : Z)
        else:
            coords = tuple(coords)

        # Sanitise the input coordinates
        if len(coords) == 1:
            coords += R.one(),
        if len(coords) != 2:
            raise ValueError('not a point on ℙ¹')
        coords = tuple(map(R, coords))

        # TODO: we should make sure the coordinates
        #       are on the curve!
        self._base_ring = R
        self._parent = parent
        self._X, self._Z = coords

    def __repr__(self):
        r"""
        """
        return f'Kummer Point [{self._X} : {self._Z}] on {self._parent}'

    def __bool__(self):
        r"""
        """
        return bool(self._Z)
    
    def __eq__(self, other):
        if self._parent != other._parent:
            return False
        return self._X * other._Z == other._X * self._Z

    def is_zero(self):
        return self._Z.is_zero()

    def base_ring(self):
        r"""
        """
        return self._base_ring
    
    def parent(self):
        r"""
        """
        return self._parent
    
    def XZ(self):
        return self._X, self._Z

    def x(self):
        r"""
        """
        # TODO: how to handle inft?
        if not self._Z:
            return self.base_ring().zero()
        if self._Z.is_one():
            return self._X
        return self._X / self._Z
    
    @cached_method
    def curve_point(self):
        E = self.parent().curve()
        P = E.lift_x(self.x())
        return min(P, -P)

    # =================================== #
    # Addition and multiplication helpers #
    # =================================== #

    @staticmethod
    def xDBL(X, Z, A, C):
        """
        function for Montgomery doubling with projective curve constant
        input:  projective point P=(X:Z), curve constants (A:C) 
        output: projective point [2]P=(X2:Z2)

        cost: 4M+2S+8a
        """

        t0 = X - Z      
        t1 = X + Z
        t0 = t0**2
        t1 = t1**2
        Z2 = C * t0
        Z2 = Z2 + Z2
        Z2 = Z2 + Z2
        X2 = Z2 * t1
        t1 = t1 - t0
        t0 = C + C
        t0 = t0 + A
        t0 = t0 * t1
        Z2 = Z2 + t0
        Z2 = Z2 * t1
        
        return X2, Z2

    @staticmethod
    def xADD(XP, ZP, XQ, ZQ, xPQ, zPQ):
        """
        function for Montgomery differential addition
        input:  projective coordinates P=(XP:ZP), Q=(XQ:ZQ), and their difference P-Q=(xPQ:zPQ) 
        output: coordinates of sum P+Q=(XQP:ZQP)

        cost: 4M+2S+6a
        """
        t0 = XP + ZP     
        t1 = XP - ZP   
        XP = XQ - ZQ   
        ZP = XQ + ZQ   
        t0 = XP * t0   
        t1 = ZP * t1   
        ZP = t0 - t1   
        XP = t0 + t1   
        ZP = ZP**2     
        XQP = XP**2    
        ZQP = xPQ * ZP 
        XQP = XQP * zPQ
        
        return XQP, ZQP

    @staticmethod
    def xDBLADD(XP,ZP,XQ,ZQ,xPQ,zPQ,A24,C24):
        """
        function for step in Montgomery ladder
        simultaneous doubling and differential addition
        input: projective coordinates P=(XP:ZP) and Q=(XQ:ZQ), 
               projective difference P-Q=(xPQ:zPQ) and 
               curve constant A24/C24=(A+2C)/4C.   
        output: projective coordinates of 2P=(X2P:Z2P)
                and Q+P=(XQP:ZQP)

        cost: 4S+8M+8A
        """
        
        t0 = XP + ZP                  
        t1 = XP - ZP 
        X2P = t0**2
        t2 = XQ - ZQ
        XQP = XQ + ZQ
        t0 = t0 * t2
        Z2P = t1**2
        t1 = t1 * XQP
        t2 = X2P - Z2P
        Z2P = Z2P * C24
        X2P = X2P * Z2P
        XQP = A24 * t2
        ZQP = t0 - t1
        Z2P = XQP + Z2P
        XQP = t0 + t1
        Z2P = Z2P * t2
        ZQP = ZQP**2
        XQP = XQP**2
        ZQP = xPQ * ZQP
        XQP = XQP * zPQ

        return X2P, Z2P, XQP, ZQP

    # =========================== #
    # Addition and multiplication #
    # =========================== #
    def _double(self):
        """
        """
        X, Z = self.XZ()
        A, C = self._parent.extract_constants()
        X2, Z2 = self.xDBL(X, Z, A, C)
        return self._parent((X2, Z2))

    def double(self):
        """
        """
        # Deal with identity
        if not self._Z:
            return self
        return self._double()
    
    def _add(self, Q, PQ):
        """
        """
        
        XP, ZP = self.XZ()
        XQ, ZQ = Q.XZ()
        XPQ, ZPQ = PQ.XZ()

        X_new, Z_new = self.xADD(XP, ZP, XQ, ZQ, XPQ, ZPQ)
        return self._parent((X_new, Z_new))

    def add(self, Q, PQ):
        """
        """
        # Adding O + Q = Q
        if not self._Z:
            return Q

        # Adding P + O = P
        if not Q._Z:
            return self

        # Difference is the identity
        # so P = Q and P+Q = [2]P
        if not PQ._Z:
            return self._double()

        return self._add(Q, PQ)
    
    def __mul__(self, m):
        """
        montgomery-ladder
        input: coordinates of P=(XP:ZP) 
               scalar factor m, curve constants (A:C)
        output: KummerPoint [m]P=(X0:Z0)
        """
        
        if not isinstance(m, (int, Integer)):
            raise ValueError("TODO")
        
        if not m:
            return self.parent().zero()
        # [m]P = [-m]P for x-only
        m = abs(m)
       
        R = self.base_ring()
        XP, ZP = self.XZ()
        
        # Initialise for loop
        X0, Z0 = R.one(), R.zero()
        X1, Z1 = XP, ZP
                
        # converting parameters for projective DBLADD -> (A24:C24)=(A+2C:4C)
        A, C = self.parent().extract_constants()
        A24 = C + C                     
        C24 = A24 + A24
        A24 = A24 + A

        # Montgomery-ladder
        for bit in bin(m)[2:]:
            if bit == "0":
                X0, Z0, X1, Z1 = self.xDBLADD(X0, Z0, X1, Z1, XP, ZP, A24, C24)
            else:
                X1, Z1, X0, Z0 = self.xDBLADD(X1, Z1, X0, Z0, XP, ZP, A24, C24)  
                
        return self._parent((X0, Z0))   
    
    def __rmul__(self, m):
        return self * m
    
    def __imul__(self, m):
        self = self * m
        return self
    
    def ladder_3_pt(self, xP, xPQ, m):
        r"""
        (self = xQ)
        Computes xP + [m]*xQ with x-only arithm.
        """
        if not isinstance(m, (int, Integer)):
            try:
                m = Integer(m)
            except:
                raise ValueError("TODO")
        
        if not m:
            return xP
        # [m]P = [-m]P for x-only
        m = abs(m)

        # converting parameters for projective DBLADD -> (A24:C24)=(A+2C:4C)
        A, C = self.parent().extract_constants()
        A24 = C + C                     
        C24 = A24 + A24
        A24 = A24 + A

        # Extract out coordinates
        XQ, ZQ = self.XZ()
        XP, ZP = xP.XZ()
        XPQ, ZPQ = xPQ.XZ()

        for bit in bin(m)[:1:-1]:
            if bit == "1":
                XQ, ZQ, XP, ZP = self.xDBLADD(XQ, ZQ, XP, ZP, XPQ, ZPQ, A24, C24)
            else:
                XQ, ZQ, XPQ, ZPQ = self.xDBLADD(XQ, ZQ, XPQ, ZPQ, XP, ZP, A24, C24)  
        return  self._parent((XP, ZP))
        
    def multiples(self):
        """
        TODO: this could be made more efficient with double add?
              for some reason when I try and implement this I get spaghetti! 
        """
        yield self
        R = self.double()
        # Order 2 case
        if not R:
            return
        
        # Odd order case
        Q = self
        while R:
            yield R
            S = R.add(self, Q)
            Q, R = R, S

        return