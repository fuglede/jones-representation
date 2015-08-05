###################################################
#
# Methods for computing the Jones representations
# of certain specific braids. For instance,
# the function matrix_to_appendix_plots() outputs
# the matrices used for some of the plots of the
# appendix of 1402.6059, reformatted for use in
# Mathematica.
#
###################################################

load("curverep.sage")


def matrix_to_mathematica(mat):
    string = "{"
    for i in range(0, mat.nrows()):
        string = string + "{"
        for j in range(0, mat.ncols()):
            string = string + str(mat[i, j])
            if (j < mat.ncols()-1):
                string = string + ","
        string = string + "}"
        if (i < mat.nrows()-1):
            string = string + ","
    string = string + "}"
    return string


def matrices_for_appendix_plots():
    B = BraidGroup(3)
    b = B([1, -2])
    print "mat31 = " + matrix_to_mathematica(simplify(b.TL_matrix(1)))

    B = BraidGroup(4)
    b = B([1, 2, -3])
    print "mat40 = " + matrix_to_mathematica(simplify(b.TL_matrix(0)))
    print "mat42 = " + matrix_to_mathematica(simplify(b.TL_matrix(2)))

    B = BraidGroup(5)
    b = B([1, 2, 3, 1, 2, 3, 4, -3])
    print "mat51 = " + matrix_to_mathematica(simplify(b.TL_matrix(1)))
    print "mat53 = " + matrix_to_mathematica(simplify(b.TL_matrix(3)))

    B = BraidGroup(6)
    b = B([2, 1, 2, 1, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5])
    print "mat60 = " + matrix_to_mathematica(simplify(b.TL_matrix(0)))
    print "mat62 = " + matrix_to_mathematica(simplify(b.TL_matrix(2)))
    print "mat64 = " + matrix_to_mathematica(simplify(b.TL_matrix(4)))

    B = BraidGroup(7)
    b = B([-4, -4, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6])
    print "mat71 = " + matrix_to_mathematica(simplify(b.TL_matrix(1)))
    print "mat73 = " + matrix_to_mathematica(simplify(b.TL_matrix(3)))
    print "mat75 = " + matrix_to_mathematica(simplify(b.TL_matrix(5)))

    B = BraidGroup(8)
    b = B([3, 4, 5, 6, 7] + [1, 2, 3, 4, 5, 6, 7]*4)
    print "mat80 = " + matrix_to_mathematica(simplify(b.TL_matrix(0)))
    print "mat82 = " + matrix_to_mathematica(simplify(b.TL_matrix(2)))
    print "mat84 = " + matrix_to_mathematica(simplify(b.TL_matrix(4)))
    print "mat86 = " + matrix_to_mathematica(simplify(b.TL_matrix(6)))


def browns_pA():
    B = BraidGroup(6)
    eta = B([4, 4, 5, 4, 4, 5, 5, 4, 4, 5, 4, 4])
    delta = B([4, 5, 3, 2, 1, -2, -3, -5, -4])
    gamma = delta*B([3, -2, -1])*delta*B([1, 2, -3])*(delta**(-1))
    phi = eta*gamma*(eta**(-1))*(gamma**(-1))

    for d in [0,2,4]:
        etamatrix = simplify(eta.TL_matrix(0))
        etamatrixinv = simplify((eta**(-1)).TL_matrix(0))
        gammamatrix = simplify(gamma.TL_matrix(0))
        gammamatrixinv = simplify((gamma**(-1)).TL_matrix(0))
        product = simplify(etamatrix*gammamatrix*etamatrixinv*gammamatrixinv)
        print "brown6"+str(d)+" = " + matrix_to_mathematica(product)


def bigelow_element():
    B = BraidGroup(5)
    psi1 = B([-3, 2, 1, 1, 2, 4, 4, 4, 3, 2])
    psi2 = B([-4, 3, 2, -1, -1, 2, 1, 1, 2, 2, 1, 4, 4, 4, 4, 4])
    el1 = (psi1**(-1))*B([4])*psi1
    el2 = (psi2**(-1))*B([4, 3, 2, 1, 1, 2, 3, 4])*psi2

    # el1matrix0 = simplify(TL_matrix(el1,0))
    # el1matrixinv0 = simplify(TL_matrix(el1**(-1), 0))
    # el2matrix0 = simplify(TL_matrix(el2,0))
    # el2matrixinv0 = simplify(TL_matrix(el2**(-1), 0))
    # matrix0 = simplify(el1matrix0*el2matrix0)
    # andenmatrix0 = simplify(el1matrixinv0*el2matrixinv0)
    # print "bigelow60 = " + matrix_to_mathematica(simplify(matrix0*andenmatrix0))

    el1matrix2 = simplify(el1.TL_matrix(1))
    el1matrixinv2 = simplify((el1**(-1)).TL_matrix(1))
    el2matrix2 = simplify(el2.TL_matrix(1))
    el2matrixinv2 = simplify((el2**(-1)).TL_matrix(1))
    matrix2 = simplify(el1matrix2*el2matrix2)
    andenmatrix2 = simplify(el1matrixinv2*el2matrixinv2)
    print "bigelow641 = " + matrix_to_mathematica(simplify(matrix2))
    print "bigelow642 = " + matrix_to_mathematica(simplify(andenmatrix2))
