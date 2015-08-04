###################################################
#
# Notes on usage:
# The function curve_rep(b,d) returns the image
# of b under the rescaled representation
#   \sigma_i \mapsto A \eta_A^{n,d}(sigma_i)
# Here, b is an element in B_n, using the existing
# braid group implementation in sage.
#
# Example of use: To evaluate
#   \eta_A^{3,1}(\sigma_1 \sigma_2^{-1})
# use
#   sage: d = 1
#   sage: B = BraidGroup(3)
#   sage: b = B([1, -2])
#   sage: curve_rep(b, d)
#
###################################################

from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.matrix.constructor import identity_matrix, matrix
from sage.groups.braid import Braid, BraidGroup, BraidGroup_class


def dim_of_rep_space(n, d):
    if mod(n+d, 2) == 1:
        raise ValueError("Arguments must have similar parity")
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
    # if len(newforest) == dim_of_rep_space(treesize + 1, forest[0][0]):
    if len(newforest[0]) == treesize:
        return newforest
    else:
        return fill_out_forest(newforest, treesize)


def basis_with_drain(num_of_strands, drain_size):
    n = num_of_strands
    d = drain_size
    if mod(n+d, 2) == 1:
        raise ValueError("Arguments must have similar parity")
    basis = [[d]]  # Let's start out with no elements and recursively fill out
    forest = fill_out_forest(basis, n-1)
    for tree in forest:
        tree.extend([1, 0])
    return forest


def create_rep(num_of_strands, drain_size, variab='A', ring=IntegerRing()):
    n = num_of_strands
    d = drain_size
    if mod(n+d, 2) == 1:
        raise ValueError("Arguments must have similar parity")
    basis = basis_with_drain(n, d)
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


def curve_rep(braid, drain_size, var='A', ring=IntegerRing()):
    R = LaurentPolynomialRing(ring, var)
    A = R.gens()[0]
    n = braid.strands()
    d = drain_size
    # TODO: Make this part dynamic. For now, the entire
    # representation is recalculated in every run. It's
    # fast enough to not be a problem in many use cases,
    # but it's still a bit silly.
    rep = create_rep(n, d, var, ring)
    M = identity_matrix(R, dim_of_rep_space(n, d))
    for i in braid.Tietze():
        if i > 0:
            M = M*rep[i-1]
        if i < 0:
            M = M*rep[-i-1]**(-1)
    return M
