from sage.all import (
    ZZ, 
    Matrix,
    Integer,
    factor, 
    prod, 
    cached_function,
    proof,
)
from random import randint
proof.all(False)


# ================================================== #
#  Code to check whether a group element has order D #
# ================================================== #

def batch_cofactor_mul_generic(G_list, pis, group_action, lower, upper):
    """
    Input:  A list of elements `G_list`, such that
                G is the first entry and the rest is empty
                in the sublist G_list[lower:upper]
            A list `pis` of primes p such that
                their product is D
            The `group_action` of the group
            Indices lower and upper
    Output: None
`
    NOTE: G_list is created in place
    """

    # check that indices are valid
    if lower > upper:
        raise ValueError(f"Wrong input to cofactor_multiples()")

    # last recursion step does not need any further splitting
    if upper - lower == 1:
        return

    # Split list in two parts,
    # multiply to get new start points for the two sublists,
    # and call the function recursively for both sublists.
    mid = lower + (upper - lower + 1) // 2
    cl, cu = 1, 1
    for i in range(lower, mid):
        cu = cu * pis[i]
    for i in range(mid, upper):
        cl = cl * pis[i]
    # cl = prod(pis[lower:mid])
    # cu = prod(pis[mid:upper])

    G_list[mid] = group_action(G_list[lower], cu)
    G_list[lower] = group_action(G_list[lower], cl)

    batch_cofactor_mul_generic(G_list, pis, group_action, lower, mid)
    batch_cofactor_mul_generic(G_list, pis, group_action, mid, upper)


@cached_function
def has_order_constants(D):
    """
    Helper function, finds constants to
    help with has_order_D
    """
    D = ZZ(D)
    pis = [p for p, _ in factor(D)]
    D_radical = prod(pis)
    Dtop = D // D_radical
    return Dtop, pis


def has_order_D(G, D, multiplicative=False):
    """
    Given an element G in a group, checks if the
    element has order exactly D. This is much faster
    than determining its order, and is enough for 
    many checks we need when computing the torsion
    bonus.

    We allow both additive and multiplicative groups
    """
    # For the case when we work with elements of Fp^k
    if multiplicative:
        group_action = lambda a, k: a**k
        is_identity = lambda a: a.is_one()
        identity = G.parent()(1)
    # For the case when we work with elements of E / Fp^k
    else:
        group_action = lambda a, k: k * a
        is_identity = lambda a: a.is_zero()
        identity = G.curve()(0)

    if is_identity(G):
        return False

    D_top, pis = has_order_constants(D)

    # If G is the identity after clearing the top
    # factors, we can abort early
    Gtop = group_action(G, D_top)
    if is_identity(Gtop):
        return False

    G_list = [identity for _ in range(len(pis))]
    G_list[0] = Gtop

    # Lastly we have to determine whether removing any prime 
    # factors of the order gives the identity of the group
    if len(pis) > 1:
        batch_cofactor_mul_generic(G_list, pis, group_action, 0, len(pis))
        if not all([not is_identity(G) for G in G_list]):
            return False

    return True

# =========================================== #
# Compute points of order D and Torsion Bases #
# =========================================== #

def generate_random_point(E, seed=None):
    """
    E0.random_point() returns either P
    or -P with equal probability.

    We always select the element with
    smaller y-coordinate to make this
    deterministic.
    """
    # Allow a setting of the seed to
    # ensure the same point is always returned
    if seed is not None:
        set_random_seed(seed)

    P = E.random_element()
    return min(P, -P)


def generate_point_order_D(E, D):
    """
    Input:  An elliptic curve E / Fp2
            An integer D dividing (p +1)
    Output: A point P of order D.
    """
    p = E.base().characteristic()
    n = (p + 1) // D
    for _ in range(1000):
        P = n * generate_random_point(E)

        # Case when we randomly picked
        # a point in the n-torsion
        if P.is_zero():
            continue

        # Check that P has order exactly D
        if has_order_D(P, D):
            P._order = ZZ(D)
            return P

    raise ValueError(f"Never found a point of order D...")


def generate_linearly_independent_point(E, P, D):
    """
    Input:  An elliptic curve E / Fp2
            A point P ∈ E[D]
            An integer D dividing (p +1)
    Output: A point Q such that E[D] = <P, Q>
    """
    for _ in range(2000):
        # Generate a random point of order D
        Q = generate_point_order_D(E, D)

        # Make sure the point is linearly independent
        pair = P.weil_pairing(Q, D, algorithm="pari")
        if has_order_D(pair, D, multiplicative=True):
            Q._order = ZZ(D)
            return Q

    raise ValueError(f"Never found a linearly independent point...")


def torsion_basis(E, D):
    """
    Generate basis of E(Fp^4)[D] of supersingular curve

    Optional   canonical: bool
               Makes torsion generation deterministic.
    """
    p = E.base().characteristic()

    # Ensure D divides the curve's order
    if (p + 1) % D != 0:
        print(f"{factor(D) = }")
        print(f"{factor(p+1) = }")
        raise ValueError(f"D must divide the point's order")

    P = generate_point_order_D(E, D)
    Q = generate_linearly_independent_point(E, P, D)

    return P, Q