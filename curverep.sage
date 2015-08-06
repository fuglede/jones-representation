"""
The Temperley-Lieb/Jones representation of braid groups

Methods for handling the representations of braid groups obtained by their
actions on Temperley-Lieb (TL) algebras and in turn the irreducible
representations of the TL algebras. The irreducible factors are given in terms
of Kauffman's diagrammatic bases in which one considers non-crossing diagrams
connecting (n+d) points, where n is the number of strands and d is a number
of 'drains' (see [Jon]_).

This in turn may be used to easily calculate the Jones polynomial of the
trace closure of a braid.

EXAMPLES:

    To calculate the matrix representation of, in the notation of [EJ]_,
    $\eta_A^{3,1}(\sigma_1 \sigma_2^{-1})$, one may proceed as follows::

        sage: d = 1
        sage: B = BraidGroup(3)
        sage: b = B([1, -2])
        sage: b.TL_matrix(d)
        [(A^8 - A^4)/(-A^4)         A^2/(-A^4)]
        [     (-A^2)/(-A^4)           1/(-A^4)]

    The trace closure of this particular braid ``b`` is the unknot whose Jones
    polynomial may now be evaluated::

        sage: b.jones_polynomial()
        1

REFERENCES:

- [Jon] Vaughan Jones. The Jones Polynomial.
        https://math.berkeley.edu/~vfr/jones.pdf
- [EJ] Jens Kristian Egsgaard and Søren Fuglede Jørgensen. The homological
       content of the Jones representations at $q = -1$.
       http://front.math.ucdavis.edu/1402.6059
"""

from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.matrix.constructor import identity_matrix, matrix
from sage.groups.braid import Braid, BraidGroup, BraidGroup_class


def dim_of_TL_space(self, drain_size):
    """
    Return the dimension of the TL representation summand when the number of
    drains is fixed to be ``drain_size``

    INPUT:

    - ``drain_size`` -- integer between 0 and the number of strands (both
      included)

    OUTPUT:

    The dimension of the representation corresponding to the number of drains
    given by ``drain_size``.

    EXAMPLES::
        sage: B = BraidGroup(8)
        sage: B.dim_of_TL_space(2)
        28
    """
    n = self.strands()
    d = drain_size
    if d > n:
        raise ValueError("Number of drains may not exceed number of strands")
    if mod(n+d, 2) == 1:
        raise ValueError("Parity of strands and drains must agree")
    if mod(n, 2) == 0:
        k = n/2
        g = k-1
        m = (2*k-d)/2
        return binomial(2*g+1, m) - binomial(2*g+1, m-2)
    else:
        k = (n-1)/2
        g = k
        m = (2*k+1-d)/2
        return binomial(2*g, m) - binomial(2*g, m-2)
BraidGroup_class.dim_of_TL_space = dim_of_TL_space


def TL_basis_with_drain(self, drain_size):
    """
    Return the basis elements given by non-intersecting pairings of $n+d$
    points in a square with $n$ points marked 'on the top' and $d$ points 'on
    the bottom' so that every bottom point is paired with a top point. Here,
    $n$ is the number of strands of the braid group, and $d$ is specified
    by ``drain_size``.

    A basis element is specified as a list of integers obtained by considering
    the pairings as obtained as the 'highest term' of trivalent trees marked by
    Jones--Wenzl projectors (see e.g. [Wan]_). In practice, this is a list of
    non-negative integers whose first element is ``drain_size``, whose last
    element is $0$, and satisfying that consecutive integers have difference
    $1$. Moreover, the length of each basis element is $n+1$.

    Given these rules, the list of lists is constructed recursively in the
    natural way.

    INPUT:

    - ``drain_size`` -- integer between 0 and the number of strands (both
      included)

    OUTPUT:

    A list of basis elements, each of which is a list of integers.

    EXAMPLES::

        sage: B = BraidGroup(5)
        sage: B.TL_basis_with_drain(3)
        [[3, 4, 3, 2, 1, 0],
         [3, 2, 3, 2, 1, 0],
         [3, 2, 1, 2, 1, 0],
         [3, 2, 1, 0, 1, 0]]

    REFERENCES:

    - [Wan] Zhenghan Wang. Tolological quantum computation. Providence,
      RI: American Mathematical Society (AMS), 2010. ISBN 978-0-8218-4930-9
    """
    def fill_out_forest(forest, treesize):
        if len(forest) == 0:
            raise ValueError("Forest has to start with a tree")
        if mod(forest[0][0]+treesize, 2) == 0:
            raise ValueError("Parity mismatch in forest creation")
        # Loop over all trees
        newforest = copy(forest)
        for tree in forest:
            if len(tree) < treesize:
                newtreeup = copy(tree)
                newtreedown = copy(tree)
                newforest.remove(tree)  # Cut down the original tree
                # Add two greater trees, admissibly
                if tree[-1] < treesize - len(tree) + 1:
                    newtreeup.append(tree[-1] + 1)
                    newforest.append(newtreeup)
                if tree[-1] > 0:
                    newtreedown.append(tree[-1] - 1)
                    newforest.append(newtreedown)
        # Are we there yet?
        if len(newforest[0]) == treesize:
            return newforest
        else:
            return fill_out_forest(newforest, treesize)

    n = self.strands()
    d = drain_size
    if d > n:
        raise ValueError("Number of drains must not exceed number of strands")
    if mod(n+d, 2) == 1:
        raise ValueError("Parity of strands and drains must agree")
    basis = [[d]]  # Let's start out with no elements and recursively fill out
    forest = fill_out_forest(basis, n-1)
    for tree in forest:
        tree.extend([1, 0])
    return forest
BraidGroup_class.TL_basis_with_drain = TL_basis_with_drain


def create_TL_rep(self, drain_size, variab='A', ring=IntegerRing()):
    """
    Calculate the matrices of the Temperley--Lieb--Jones representation of
    the standard generators of the braid groups in the basis given by
    non-intersecting pairings of $(n+d)$ points, where $n$ is the number of
    strands, and $d$ is given by ``drain_size'', and the pairings satisfy
    certain rules. This basis has the useful property that all resulting
    entries can be regarded as Laurent polynomials.

    We use the convention that the eigenvalues of the standard generators are
    $1$ and $-A^4$, where $A$ is the generator of the Laurent polynomial ring.

    When $d = n-2$ and the variables are picked appropriately, the resulting
    representation is equivalent to the reduced Burau representation.

    Store the result of the calculation as part of the braid group.

    INPUT:

    - ``drain_size`` -- integer between 0 and the number of strands (both
      included)
    - ``variab`` -- string (default: ``'A'``); the name of the
      variable in the entries of the matrices
    - ``ring`` -- ring (default: ``IntegerRing()``); the ring to which the
      coefficients of the polynomial entries belong

    OUTPUT:

    A list of matrices corresponding to the representations of each of the
    standard generators.

    EXAMPLES::
        sage: B = BraidGroup(4)
        sage: B.create_TL_rep(0)
        [
        [   1    0]  [-A^4  A^2]  [   1    0]
        [ A^2 -A^4], [   0    1], [ A^2 -A^4]
        ]

    REFERENCES:

    - [Jon] Vaughan Jones. The Jones Polynomial.
            https://math.berkeley.edu/~vfr/jones.pdf
    """
    n = self.strands()
    d = drain_size
    basis = self.TL_basis_with_drain(d)
    auxmat = matrix(n-1, len(basis))
    for i in range(1, n):
        for v in range(len(basis)):
            tree = basis[v]
            if tree[i-1] < tree[i] and tree[i+1] < tree[i]:
                auxmat[i-1, v] = v
            if tree[i-1] > tree[i] and tree[i+1] > tree[i]:
                newtree = copy(tree)
                newtree[i] += 2
                auxmat[i-1, v] = basis.index(newtree)
            if tree[i-1] > tree[i] and tree[i+1] < tree[i]:
                newtree = copy(tree)
                newtree[i-1] -= 2
                j = 2
                while newtree[i-j] != newtree[i] and i-j >= 0:
                    newtree[i-j] -= 2
                    j += 1
                if newtree in basis:
                    auxmat[i-1, v] = basis.index(newtree)
                else:
                    auxmat[i-1, v] = -1
            if tree[i-1] < tree[i] and tree[i+1] > tree[i]:
                newtree = copy(tree)
                newtree[i+1] -= 2
                j = 2
                while newtree[i+j] != newtree[i] and i+j <= n:
                    newtree[i+j] -= 2
                    j += 1
                if newtree in basis:
                    auxmat[i-1, v] = basis.index(newtree)
                else:
                    auxmat[i-1, v] = -1
    repmat = []
    R = LaurentPolynomialRing(ring, variab)
    A = R.gens()[0]
    for i in range(1, n):
        repmatnew = identity_matrix(R, len(basis))
        for v in range(len(basis)):
            newmatentry = auxmat[i-1, v]
            if newmatentry == v:
                repmatnew[v, v] = -A**4
            elif newmatentry >= 0:
                repmatnew[newmatentry, v] = A**2
        repmat.append(repmatnew)
    return repmat
BraidGroup_class.create_TL_rep = create_TL_rep


def TL_matrix(self, drain_size, variab='A', ring=IntegerRing()):
    """
    Calculate the matrices of the Temperley--Lieb--Jones representation of
    the braidin the basis given by non-intersecting pairings of $(n+d)$ points,
    where $n$ is the number of strands, and $d$ is given by ``drain_size'',
    and the pairings satisfy certain rules.

    We use the convention that the eigenvalues of the standard generators are
    $1$ and $-A^4$, where $A$ is the generator of the Laurent polynomial ring.

    When $d = n-2$ and the variables are picked appropriately, the resulting
    representation is equivalent to the reduced Burau representation.

    INPUT:

    - ``drain_size`` -- integer between 0 and the number of strands (both
      included)
    - ``variab`` -- string (default: ``'A'``); the name of the
      variable in the entries of the matrices
    - ``ring`` -- ring (default: ``IntegerRing()``); the ring to which the
      coefficients of the polynomial entries belong

    OUTPUT:

    The matrix of the TL representation of the braid.

        sage: B = BraidGroup(4)
        sage: b = B([1, 2, -3])
        sage: b.TL_matrix(0)
        [(A^8 - A^4)/(-A^4)         A^2/(-A^4)]
        [              -A^6                  0]

    REFERENCES:

    - [Jon] Vaughan Jones. The Jones Polynomial.
            https://math.berkeley.edu/~vfr/jones.pdf
    """
    R = LaurentPolynomialRing(ring, variab)
    A = R.gens()[0]
    n = self.strands()
    d = drain_size
    B = BraidGroup(n)
    # It is worth noting that making create_TL_rep dynamic seems like it would
    # provide for faster evaluation. In practice it appears not to, though,
    # as create_TL_rep is significantly faster than matrix multiplication.
    rep = B.create_TL_rep(d, variab, ring)
    M = identity_matrix(R, B.dim_of_TL_space(d))
    for i in self.Tietze():
        if i > 0:
            M = M*rep[i-1]
        if i < 0:
            M = M*rep[-i-1]**(-1)
    return M
Braid.TL_matrix = TL_matrix


def exponent_sum(self):
    """
    Return the exponent sum of the braid.

    OUTPUT:

    Integer.

    EXAMPLES::

        sage: B = BraidGroup(5)
        sage: b = B([1, 4, -3, 2])
        sage: b.exponent_sum()
        2
    """
    tietze = self.Tietze()
    return sum([sign(s) for s in tietze])
Braid.exponent_sum = exponent_sum


def components_in_closure(self):
    """
    Return the number of components of the trace closure of the braid.

    OUTPUT:

    Integer.

    EXAMPLES::

        sage: B = BraidGroup(5)
        sage: b = B([1, -3])
        sage: b.components_in_closure()
        3
    """
    n = self.strands()
    perm = self.permutation()
    cycles = perm.to_cycles(singletons=False)
    return n-sum([len(c)-1 for c in cycles])
Braid.components_in_closure = components_in_closure


def markov_trace(self, variab='A', ring=IntegerRing()):
    """
    Calculate the Markov trace of the braid. The normalisation is so that in
    the underlying braid group representation, the eigenvalues of the standard
    generators of the braid group are $1$ and $-A^4$.

    INPUT:

    - ``variab`` -- string (default: ``'A'``); the name of the variable in the
      resulting Laurent polynomial
    - ``ring`` -- ring (default: ``IntegerRing()``); the ring to which the
      coefficients of the polynomial entries belong

    OUTPUT:

    Quotient of Laurent polynomials over ``ring`` in the variable ``variab``.

    EXAMPLES::

        sage: B = BraidGroup(4)
        sage: b = B([1, 2, -3])
        sage: b.markov_trace().factor()
        (A^4) * (A^4 + 1)^-3

    REFERENCES:

    - [Jon] Vaughan Jones. The Jones Polynomial.
            https://math.berkeley.edu/~vfr/jones.pdf
    """
    def qint(i, variab='A', ring=IntegerRing()):
        R = LaurentPolynomialRing(ring, variab)
        A = R.gens()[0]
        return (A**(2*i) - A**(-2*i))/(A**2 - A**(-2))

    def weighted_trace(b, d, variab='A', ring=IntegerRing()):
        return qint(d+1, variab, ring)*b.TL_matrix(d, variab, ring).trace()

    R = LaurentPolynomialRing(ring, variab)
    A = R.gens()[0]
    delta = -A**2 - A**(-2)
    n = self.strands()
    drains = [d for d in range(n+1) if mod(n+d, 2) == 0]
    traces = [weighted_trace(self, d, variab, ring) for d in drains]
    return sum(traces)/((-delta)**n)
Braid.markov_trace = markov_trace


def jones_polynomial(self, skein_variable=True):
    """
    Return the Jones polynomial of the trace closure of the braid, normalised
    so that the unknot has Jones polynomial $1$. If ``skein_variable'' is True,
    give the result in terms of a variable ``'A'`` so that the result agrees
    with the conventions of [Lic]_ (which in particular differs slightly from
    the conventions used otherwise in this class). If ``skein_variable'' is
    False, return the result in terms of the variable ``'t'``, also used
    in [Lic]_.

    The computation uses the representation of the braid group on the
    Temperley--Lieb algebra.

    INPUT:

    - ``skein_variable`` -- boolean (default: ``True``); determines the
      variable of the resulting polynomial.

    OUTPUT:

    Laurent polynomial in the variable ``'A'``or ``'t'`` depending on the value
    of ``skein_variable``. Might have fractional powers if ``skein_variable``
    is False and the closure of the braid is not a know.

    EXAMPLES:

    The unknot::

        sage: B = BraidGroup(9)
        sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
        sage: b.jones_polynomial()
        1

    Two different representations of the trefoil and one of its mirror::

        sage: B = BraidGroup(2)
        sage: b = B([1, 1, 1])
        sage: b.jones_polynomial()
        1/A^4 + 1/A^12 - 1/A^16
        sage: B = BraidGroup(3)
        sage: b = B([1, 2, 1, 2])
        sage: b.jones_polynomial()
        1/A^4 + 1/A^12 - 1/A^16
        sage: B = BraidGroup(3)
        sage: b = B([-1, -2, -1, -2])
        sage: b.jones_polynomial()
        -A^16 + A^12 + A^4

    K11n42 (the mirror of the "Kinoshita-Terasaka" knot) and K11n34 (the
    mirror of the "Conway" knot)::

        sage: B = BraidGroup(4)
        sage: b11n42 = B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2])
        sage: b11n34 = B([1, 1, 2, -3, 2, -3, 1, -2, -2, -3, -3])
        sage: cmp(b11n42.jones_polynomial(), b11n34.jones_polynomial())
        0

    REFERENCES:

    - [Lic] William B. Raymond Lickorish. An introduction to knot theory,
            volume 175 of Graduate Texts in Mathematics. Springer-Verlag, New
            York, 1997. ISBN 0-387-98254-X
    """
    variab = 'A'
    ring = IntegerRing()
    R = LaurentPolynomialRing(ring, variab)
    A = R.gens()[0]
    delta = -A**2 - A**(-2)
    n = self.strands()
    exp_sum = exponent_sum(self)
    trace = self.markov_trace(variab, ring)
    jones_pol = (-delta)**(n-1) * A**(2*exp_sum) * trace
    jones_pol = jones_pol.factor().expand()
    if skein_variable:
        return jones_pol.subs(A=var('A')**(-1))
    else:
        return jones_pol.subs(A=-var('t')**(1/4))
Braid.jones_polynomial = jones_polynomial
