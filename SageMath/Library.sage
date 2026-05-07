from sage.all import (ZZ, matrix, Permutation, vector, det,
                      identity_matrix, zero_matrix, LatticePolytope, factor,
                      Polyhedron, binomial, factorial, ToricVariety,
                      lcm, QQ, mrange, PolynomialRing, gcd, divisors,
                      diagonal_matrix, PointConfiguration,
                      power, sign, ReflexivePolytopes, Combinations, Tuples)
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.fraction_field import FractionField_generic
from sage.rings.polynomial.flatten import FlatteningMorphism
from sage.rings.polynomial.laurent_polynomial_ring_base import (
    LaurentPolynomialRing_generic)
from sage.rings.polynomial.laurent_polynomial_ring import (
    LaurentPolynomialRing)


# Matrix computations

def matrix_reduce(matrix):
    """
    Returns the row echelon form of a matrix over ZZ with omitted zero rows.

    Arguments:
         matrix : A matrix over ZZ.

    Returns:
         A matrix over ZZ.

    Example:
        sage: M = matrix(ZZ, [[2, 0, 0],
        ....:                 [2, 1, 1],
        ....:                 [0, 1, 1]])
        sage: matrix_reduce(M)
        [2 0 0]
        [0 1 1]
    """
    # Sanity check
    if not isinstance(matrix.parent(), MatrixSpace):
        raise TypeError("The argument should be a matrix.")

    if matrix.base_ring() is not ZZ:
        raise ValueError("Matrix should be defined over ZZ.")

    return matrix.echelon_form()[:matrix.rank()]


def matrix_splitting(matrix, decomposition_list):
    """
    Given a linear system A*X = 0, construct the equivalent linear system
    of the form L*Y = R*Z, where
    Y = [X_i for i not in decomposition_list]
    Z = [X_i for i in decomposition_list].

    Arguments:
        matrix : A matrix over ZZ.
        decomposition_list :
            A list of integers in range(1, matrix.ncols() + 1)
            of length matrix.ncols() - matrix.rank().

    Returns:
        A list [I, L, R], where
        - I is the complement to the decomposition list
        - L is the left-hand side matrix
        - R is the right-hand side matrix.

    Example:
        sage: M = matrix(ZZ, ((2, 0, 2, 0, 0, 0, 0, -1),
        ....:                 (0, 2, 0, 2, 0, 0, 0, -1),
        ....:                 (0, 0, 0, 0, 0, 0, 2, -1),
        ....:                 (1, 1, 0, 0, 1, 1, 0, -1)))
        sage: L = (1, 2, 3, 5)
        sage: matrix_splitting(M, L)[0]
        [4, 6, 7, 8]
        sage: matrix_splitting(M, L)[1]
        [2 0 0 0]
        [0 1 0 0]
        [0 0 2 0]
        [0 0 0 1]
        sage: matrix_splitting(M, L)[2]
        [ 2 -2  2  0]
        [ 1 -1  2 -1]
        [ 2  0  2  0]
        [ 2  0  2  0]
    """
    # Sanity check
    if not isinstance(matrix.parent(), MatrixSpace):
        raise TypeError("The argument should be a matrix.")

    if matrix.base_ring() is not ZZ:
        raise ValueError("Matrix should be defined over ZZ.")

    for i in decomposition_list:
        if (i - 1) not in range(matrix.ncols()):
            raise ValueError("Incorrect list of indexes.")

    if len(decomposition_list) != matrix.ncols() - matrix.rank():
        raise ValueError("Incorrect list of indexes.")

    # Construct the permutation
    total_list = range(1, matrix.ncols() + 1)
    complement_list = [x for x in total_list if x not in decomposition_list]
    permutation = Permutation(complement_list +
                              list(decomposition_list)).to_matrix()

    # Compute the reduced row echelon form of the permuted matrix
    permuted_matrix = matrix*permutation
    permuted_matrix_echelon = matrix_reduce(permuted_matrix)

    # Output
    return [complement_list,
            permuted_matrix_echelon[:, :matrix.rank()],
            (-1)*permuted_matrix_echelon[:, matrix.rank():]]


def relations_matrix_bases(matrix, frac_field):
    """
    Returns the list of possible integral or rational bases
    w.r.t. a given relations matrix.

    Arguments:
        matrix : A matrix over ZZ.
        frac_field : A boolean.

    Returns:
        If frac_field == False, returns the list of integral bases.
        If frac_field == True, returns the list of rational bases.

    Example:
        sage: M = matrix(ZZ, ((2, 0, 2, 0, 0, 0, 0, -1),
        ....:                 (0, 2, 0, 2, 0, 0, 0, -1),
        ....:                 (0, 0, 0, 0, 0, 0, 2, -1),
        ....:                 (1, 1, 0, 0, 1, 1, 0, -1)))
        sage: relations_matrix_bases(M, frac_field = False)
        []
    """
    # Sanity check
    if not isinstance(matrix.parent(), MatrixSpace):
        raise TypeError("The argument should be a matrix.")

    if matrix.base_ring() is not ZZ:
        raise ValueError("Matrix should be defined over ZZ.")

    if not isinstance(frac_field, bool):
        raise TypeError("The second argument should be a Boolean.")

    # Basis matrix of ker(M)
    kernel = matrix.right_kernel()
    kernel_matrix = kernel.matrix()
    kernel_nrows = kernel_matrix.nrows()
    kernel_ncols = kernel_matrix.ncols()
    kernel_rank = kernel_matrix.rank()

    # List of basic minors of the kernel matrix
    basic_minors = [[kernel_matrix[rows, cols], cols]
                    for cols in Combinations(kernel_ncols, kernel_rank)
                    for rows in Combinations(kernel_nrows, kernel_rank)]

    # Run over Z-invertible (or Q-invertible) basic minors of the kernel matrix
    output = []
    for minor_matrix, minor_indexes in basic_minors:
        if frac_field:
            if det(minor_matrix) == 0:
                continue
        else:
            if not (det(minor_matrix) in [1, -1]):
                continue
        basis = vector(minor_indexes) + vector([1]*kernel_rank)

        # Sanity check
        LHS_matrix = matrix_splitting(matrix, basis)[1]
        if frac_field:
            if LHS_matrix.det() == 0:
                continue
        else:
            if (LHS_matrix !=
                    identity_matrix(ZZ, kernel_ncols - kernel_rank)):
                continue

        output.append(basis)

    return output


# Laurent polynomial computations

def immutable(obj):
    """
    Returns the immutable copy of a mutable object.

    Arguments:
        obj : A mutable object.
    """
    obj.set_immutable()
    return obj


def GL_generators(n) -> tuple:
    """
    Returns a tuple of generators of the group GL(n, Z).

    Returns:
        A tuple of invertible matrices over ZZ.

    Example:
        sage: GL_generators(3)
        (
        [1 0 0]  [0 0 1]  [0 1 0]  [1 1 0]  [-1  0  0]
        [0 1 0]  [1 0 0]  [1 0 0]  [0 1 0]  [ 0  1  0]
        [0 0 1], [0 1 0], [0 0 1], [0 0 1], [ 0  0  1]
        )
    """
    if n == 1:
        return tuple([identity_matrix(ZZ, n), -identity_matrix(ZZ, n)])

    s_1 = zero_matrix(ZZ, n)
    for i in range(n):
        s_1[i, i - 1] = 1

    s_2 = identity_matrix(ZZ, n)
    s_2[0, 1] = 1
    s_2[1, 0] = 1
    s_2[0, 0] = 0
    s_2[1, 1] = 0

    s_3 = identity_matrix(ZZ, n)
    s_3[0, 1] = 1

    s_4 = identity_matrix(ZZ, n)
    s_4[0, 0] = -1

    return tuple([identity_matrix(ZZ, n), s_1, s_2, s_3, s_4])


def GL_bounded(matrix_size, matrix_bound):
    """
    Returns the frozenset of immutable invertible matrices over ZZ of
    size matrix_size with the coefficients bounded by matrix_bound.

    Arguments:
        matrix_size : An integer.
        matrix_bound : An integer.

    Returns:
        A frozenset of immutable invertible matrices over ZZ.

    Example:
        sage: GL_bounded(1, 1)
        frozenset({[-1], [1]})
    """
    # Sanity check
    if not ((matrix_size in ZZ) and matrix_size > 0):
        raise ValueError("Incorrect matrix size.")

    if not ((matrix_bound in ZZ) and matrix_bound > 0):
        raise ValueError("Incorrect matrix bound.")

    # Construct the set of invertible bounded matrices
    matrix_set = set()
    for i in Tuples(range(-matrix_bound, matrix_bound + 1), matrix_size**2):
        m = matrix(ZZ, matrix_size, matrix_size, i)
        if m.is_invertible():
            matrix_set.add(immutable(m))

    return frozenset(matrix_set)


def matrix_words(generators_tuple, words_length):
    """
    Returns a frozenset of immutable matrices obtained as a words of length
    words_len in a tuple of matrices generators_tuple and its inverses.

    Arguments:
        generators_tuple : A tuple of invertible matrices of same size.
        words_length : A positive integer.

    Returns:
        A frozenset of immutable invertible matrices.

    Example:
        sage: matrix_words(GL_generators(3), 1)
        frozenset({[-1  0  0]
                   [ 0  1  0]
                   [ 0  0  1],
                   ...
                   [1 1 0]
                   [0 1 0]
                   [0 0 1]})
    """
    # Sanity check
    if not isinstance(generators_tuple, tuple):
        raise TypeError("Incorrect generators tuple.")

    matrix_space = generators_tuple[0].parent()

    for i in generators_tuple:
        if not (isinstance(i.parent(), MatrixSpace) and
                (i.parent() == matrix_space)):
            raise ValueError("Incorrect generators tuple.")

    if not (words_length in ZZ and words_length > 0):
        raise ValueError("Incorrect words length.")

    # Compute all words of length 1
    words_gens = set()
    for gen in generators_tuple:
        words_gens.add(immutable(gen))
        words_gens.add(immutable(gen.inverse()))

    # Compute all words of length words_length
    words_set = words_gens
    for i in range(1, words_length):
        words_set = {immutable(m * g) for m in words_gens for g in words_set}

    return frozenset(words_set)


def reflexive_polytope_type(laurent_polynomial):
    """
    Returns the reflexive polytope that is combinatorially
    isomorphic to the Newton polytope of a Laurent polynomial.

    Arguments:
        laurent_polynomial : A Laurent polynomial.

    Returns:
        An element of ReflexivePolytopes(n) for
        n = (laurent_polynomial.parent()).ngens().

    Example:
        sage: C.<a_1,a_2> = PolynomialRing(ZZ)
        sage: R.<x,y,z> = LaurentPolynomialRing(C)
        sage: reflexive_polytope_type(x + y + a_1*x^-1*y^-1 + z + a_2*z^-1)
        3-d reflexive polytope #4 in 3-d lattice M
    """
    parent_ring = laurent_polynomial.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError("The argument should be a Laurent polynomial.")

    if parent_ring.ngens() not in (2, 3):
        raise NotImplementedError("ReflexivePolytopes(n) is not implemented "
                                  "for n != 2, 3.")

    # Construct the Newton polytope
    monomials_list = []
    for i in laurent_polynomial.monomials():
        monomials_list.append([i.degree(gen) for gen in parent_ring.gens()])
    newton_polytope = LatticePolytope(monomials_list)

    # Build the list of combinatorially isomorphic polytopes
    polytopes_list = []
    for i in ReflexivePolytopes(parent_ring.ngens()):
        if (newton_polytope.normal_form() == i.normal_form()):
            polytopes_list.append(i)

    if len(polytopes_list) == 1:
        return polytopes_list[0]

    raise ValueError("The Newton polytope of the Laurent polynomial is "
                     "not reflexive.")


def period_sequence(laurent_polynomial, sequence_length):
    """
    Returns the list of the first terms of the period sequence
    of a Laurent polynomial.

    Arguments:
        laurent_polynomial : A Laurent polynomial.
        sequence_length : A positive integer.

    Returns:
        A list of non-negative integers of length sequence_len.

    Example:
        sage: C.<a_1,a_2> = PolynomialRing(ZZ)
        sage: R.<x,y,z> = LaurentPolynomialRing(C)
        sage: period_sequence(x + y + a_1*x^-1*y^-1 + z + a_2*z^-1, 6)
        [1, 0, 2*a_2, 6*a_1, 6*a_2^2, 120*a_1*a_2]
    """
    parent_ring = laurent_polynomial.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError("The first argument should be a Laurent polynomial.")

    if not (sequence_length in ZZ and sequence_length > 0):
        raise TypeError("The second argument should be a positive integer.")

    # Build the list of the first terms of the period sequence
    period_list = []

    for i in range(sequence_length):
        polynomial_power = power(laurent_polynomial, i)
        period_list.append(polynomial_power.constant_coefficient())

    return period_list


def GL_action(function, monomial_matrix):
    """
    Returns the image of a Laurent polynomial or an element of a fraction field
    under a monomial change of variables.

    Arguments:
        function : A Laurent polynomial or an element of a fraction field.
        monomial_matrix : An invertible matrix over ZZ of size
            (function.parent()).ngens().

    Returns:
        An element of the parent Laurent polynomial ring or the fraction field.

    Example:
        sage: C.<a_1,a_2> = PolynomialRing(ZZ)
        sage: R.<x,y,z> = LaurentPolynomialRing(C)
        sage: M = diagonal_matrix(ZZ, (1, 1, -1))
        sage: GL_action(x + y + a_1*x^-1*y^-1 + z + a_2*z^-1, M)
        x + y + a_2*z + z^-1 + a_1*x^-1*y^-1
    """
    parent_ring = function.parent()
    parent_ngens = parent_ring.ngens()

    # Sanity check
    if not isinstance(parent_ring, (LaurentPolynomialRing_generic,
                                    FractionField_generic)):
        raise TypeError("The argument should be a Laurent polynomial "
                        "or an element of a fraction field.")

    matrix_space = MatrixSpace(ZZ, parent_ngens, parent_ngens)
    if monomial_matrix not in matrix_space:
        raise ValueError('Incorrect base ring of the monomial matrix.')

    # Formally define a Laurent polynomial ring
    # to construct a monomial change of variables
    laurent_ring = LaurentPolynomialRing(parent_ring.base_ring(),
                                         parent_ring.gens())
    arguments_list = [laurent_ring.monomial(*monomial_matrix[i])
                      for i in range(parent_ngens)]
    monomial_change = tuple(arguments_list)

    return function(monomial_change)


def largest_mutation_factor(laurent_polynomial):
    """
    Returns the largest mutation factor of a Laurent polynomial.

    Arguments:
        laurent_polynomial : A Laurent polynomial over ZZ, QQ, or over
            a polynomial ring over ZZ or QQ.

    Returns:
        An element of the polynomial ring of the parent
        Laurent polynomial ring.

    Example:
        sage: C.<a,b,c> = PolynomialRing(ZZ)
        sage: R.<x,y,z> = LaurentPolynomialRing(C)
        sage: P = x*y*z + x + b*y + c*z + a*x^-1 + c*x^-1*y^-1*z^-1
        sage: M = matrix(ZZ, ((1, -1, 1),
        ....:                 (0, 1, -1),
        ....:                 (0, 0, 1)))
        sage: Q = GL_action(P, M)
        sage: largest_mutation_factor(Q)
        b*x*y + a*y + c
    """
    parent_ring = laurent_polynomial.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError('Incorrect Laurent polynomial.')

    # Define the auxiliary rings
    frac_field = parent_ring.fraction_field()
    poly_ring = parent_ring.polynomial_ring()
    poly_ring_flattening = FlatteningMorphism(poly_ring)
    poly_ring_flattened = poly_ring_flattening.codomain()

    # Sanity check
    if not (poly_ring_flattened.base_ring() in (ZZ, QQ)):
        raise TypeError('Incorrect coefficient ring.')

    # Compute the minimal degree of a Laurent polynomial in the last generator
    parent_ngens = parent_ring.ngens()
    gen = parent_ring.gen(parent_ngens - 1)
    diag_matrix = diagonal_matrix(ZZ, (parent_ngens - 1)*[1] + [-1])
    inverse_poly = GL_action(laurent_polynomial, diag_matrix)
    min_degree = -inverse_poly.degree(gen)

    # If the polynomial has non-negative degree in the last generator,
    # return the trivial mutation factor
    if min_degree > -1:
        return 1

    # Compute the list of numerators
    numerators_list = []
    for i in range(min_degree, 0):
        coefficient = laurent_polynomial.coefficient(gen**i)
        numerator = frac_field(coefficient).numerator()
        numerator_flattened = poly_ring_flattening(numerator)
        numerators_list.append(numerator_flattened)

    # Compute all possible mutations
    divisors_list = []
    for divisor in divisors(gcd(numerators_list)):
        # Check that the divisor defines a mutation
        is_good = 1
        for i in range(min_degree, 0):
            divisor_flattened = poly_ring_flattening(divisor**(-i))
            remainder = numerators_list[i].mod(divisor_flattened)
            if remainder != 0:
                is_good = 0
                break
        # If it holds, append the divisor to the list
        if is_good:
            divisors_list.append(divisor)

    return lcm(divisors_list)


def laurent_polynomial_mutation(laurent_polynomial, mutation_factor):
    """
    Returns the mutation of a Laurent polynomial w.r.t. to a mutation factor.

    Arguments:
        laurent_polynomial : A Laurent polynomial.
        mutation_factor : A polynomial in the parent Laurent polynomial ring.

    Returns:
        An element of the parent Laurent polynomial ring.

    Example:
        sage: C.<a,b,c> = PolynomialRing(ZZ)
        sage: R.<x,y,z> = LaurentPolynomialRing(C)
        sage: P = x*y*z + x + b*y + c*z + a*x^-1 + c*x^-1*y^-1*z^-1
        sage: M = matrix(ZZ, ((1, -1, 1),
        ....:                 (0, 1, -1),
        ....:                 (0, 0, 1)))
        sage: Q = GL_action(P, M)
        sage: F = (b*x*y + a*y + c)(R.polynomial_ring().gens())
        sage: laurent_polynomial_mutation(Q, F)
        b*x^2*y*z + b*x^2*z + (b*c + a)*x*y*z + (a + c)*x*z + a*c*y*z +
        c*x*y^-1*z + c^2*z + x^-1*z^-1
    """
    parent_ring = laurent_polynomial.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError('Incorrect Laurent polynomial.')

    frac_field = parent_ring.fraction_field()

    # Sanity check
    if not (mutation_factor in parent_ring.polynomial_ring()):
        raise TypeError('Incorrect mutation factor.')

    # Compute the maximal and minimal degrees in the last generator
    parent_ngens = parent_ring.ngens()
    gen = parent_ring.gen(parent_ngens - 1)
    max_degree = laurent_polynomial.degree(gen)
    diag_matrix = diagonal_matrix(ZZ, [1]*(parent_ngens - 1) + [-1])
    inverse_poly = GL_action(laurent_polynomial, diag_matrix)
    min_degree = -inverse_poly.degree(gen)

    # Extract the free coefficient w.r.t. the last generator
    temp_ring = LaurentPolynomialRing(parent_ring.base_ring(),
                                      parent_ring.gens()[:-1])
    temp_ring_extended = LaurentPolynomialRing(temp_ring, [gen])
    free_coeff = temp_ring_extended(laurent_polynomial).constant_coefficient()

    # Compute the mutated Laurent polynomial
    output = 0
    for i in range(min_degree, max_degree + 1):
        if (i != 0):
            output += laurent_polynomial.coefficient(gen**i) * \
                frac_field(gen)**i * power(mutation_factor, i)
        else:
            output += free_coeff

    # Return to the parent Laurent polynomial ring
    result = output(parent_ring.gens())

    # Sanity check
    if not (result in parent_ring):
        raise TypeError('Incorrect Laurent polynomial or a mutation factor.')

    return result


def mutation_evaluation(laurent_polynomial, mutation_data):
    """
    Evaluate a mutation given by a mutation data triple on a Laurent polynomial.

    Arguments:
        laurent_polynomial : A Laurent polynomial.
        mutation_data : A tuple (A, f, B), where A and B are square invertible
            matrices over ZZ of size (laurent_polynomial.parent()).ngens(),
            and f is an element of the polynomial ring of
            the parent Laurent polynomial ring.

    Example:
        sage: R.<x,y,z> = LaurentPolynomialRing(ZZ)
        sage: F = y*z + x + y + z + x^-1 + x^-1*z^-1 + x^-1*y^-1
        sage: M = (matrix(ZZ, ((-1, -1, 1),
        ....:                  (0, 1, 0),
        ....:                  (0, 1, -1))),
        ....:      y + 1,
        ....:      matrix(ZZ, ((0, -1, -1),
        ....:                  (1, -1, 0),
        ....:                  (1, -1, -1))))
        sage: mutation_evaluation(F, M)
        x + y + z + x*y^-1 + x^-1 + y^-1*z^-1
    """
    parent_ring = laurent_polynomial.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError('Incorrect Laurent polynomial.')

    parent_ngens = parent_ring.ngens()
    poly_ring = parent_ring.polynomial_ring()

    # Sanity check
    if not (len(mutation_data) == 3):
        raise TypeError('Incorrect mutation data.')

    for i in range(3):
        if i % 2 == 0:
            if not (mutation_data[i] in MatrixSpace(ZZ, parent_ngens)):
                raise TypeError('Incorrect mutation data.')
        else:
            if not (mutation_data[1] in parent_ring):
                raise TypeError('Incorrect mutation data.')

    # Evaluate the mutation
    output = laurent_polynomial
    for i in range(3):
        if i % 2 == 0:
            output = GL_action(output, mutation_data[i])
        else:
            mutation_factor = mutation_data[i](poly_ring.gens())
            output = laurent_polynomial_mutation(output, mutation_factor)

    return output


def find_mutation_factors(laurent_polynomial_input, laurent_polynomial_output,
                          matrices_list, recursion_depth):
    """
    Finds an explicit presentation of a composition of mutations transforming
    one Laurent polynomial into another in the assumption that
    the GL-equivalences are fixed.

    Arguments:
        laurent_polynomial_input : A Laurent polynomial over ZZ or QQ.
        laurent_polynomial_output : A Laurent polynomial over ZZ or QQ.
        matrices_list : A list of invertible matrices over ZZ of size
            n = (laurent_polynomial_input.parent()).ngens().
        recursion_depth : An integer in range(len(matrices_list)).
            It always should be equal to zero for the first call.

    Returns:
        A list of the form [A, f, B, g, ...], where [A,B,...] = matrices_list
        and f,g,... are elements of the flattening of the ring
        (laurent_in.parent()).polynomial_ring().

    Example:
        sage: R.<x,y,z> = LaurentPolynomialRing(ZZ)
        sage: minkowski_0070 = y*z + x + y + z + x^-1 + x^-1*z^-1 + x^-1*y^-1
        sage: minkowski_0021 = x + y + z + x*y^-1 + x^-1 + y^-1*z^-1
        sage: M1 = matrix(ZZ, [[0, 0, 1],
        ....:                  [1, 0, -1],
        ....:                  [0, 1, 0]])
        sage: M2 = matrix(ZZ, [[0, 1, 1],
        ....:                  [1, -1, 0],
        ....:                  [0, 1, 0]])
        sage: find_mutation_factors(minkowski_0070, minkowski_0021, [M1, M2],
        ....:                       recursion_depth = 0)
        [
        [ 0  0  1]         [ 0  1  1]
        [ 1  0 -1]         [ 1 -1  0]
        [ 0  1  0], y + 1, [ 0  1  0]
        ]
    """
    parent_ring = laurent_polynomial_output.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError('Incorrect output Laurent polynomial.')

    if laurent_polynomial_input not in parent_ring:
        raise TypeError('Incorrect input Laurent polynomial.')

    if not (parent_ring.base_ring() in (ZZ, QQ)):
        raise TypeError('Incorrect coefficient ring.')

    # Recursion implementation
    if (recursion_depth == len(matrices_list) - 1):
        GL_image = GL_action(laurent_polynomial_input,
                             matrices_list[recursion_depth])
        if (GL_image == laurent_polynomial_output):
            return [matrices_list[recursion_depth]]
        else:
            return []
    else:
        GL_image = GL_action(laurent_polynomial_input,
                             matrices_list[recursion_depth])
        for divisor in divisors(largest_mutation_factor(GL_image)):
            for D in [-divisor, divisor]:
                input_mutated = laurent_polynomial_mutation(GL_image, D)
                output = find_mutation_factors(input_mutated,
                                               laurent_polynomial_output,
                                               matrices_list,
                                               recursion_depth + 1)
                if (output != []):
                    return ([matrices_list[recursion_depth], D] + output)
    return []


def find_mutation_call(args):
    """
    A wrapper for the function find_mutation_factors.

    Global variables:
        poly_input : A Laurent polynomial over a polynomial ring over ZZ or QQ.
        poly_output : A Laurent polynomial over a polynomial ring over ZZ or QQ.
        words_input : A tuple of invertible matrices over ZZ of size n
            for n = (laurent_in.parent()).ngens().

    Arguments:
        args : A tuple of the elements in range(len(words)).
    """
    return find_mutation_factors(
        *tuple([poly_input, poly_out,
                tuple([words_input[i] for i in args]), 0]))


def find_mutation(laurent_in, laurent_out, words_tuple, steps):
    """
    Finds an explicit presentation of a mutation transforming one Laurent
    polynomial into another.

    WARNING! Could be CPU and memory consuming.

    REMARK. Loads tqdm module for the progress bar.

    Arguments:
        laurent_in : A Laurent polynomial over a polynomial ring over ZZ or QQ.
        laurent_out : A Laurent polynomial over a polynomial ring over ZZ or QQ.
        words_tuple : A tuple of Z-invertible matrices of size n
            for n = (laurent_in).parent().ngens().
        steps : A non-negative integer.

    Returns:
        A list of the form [A, f, B, g, ...], where A,B,... are in words_tuple
        and f,g,... are elements of the flattening of the ring
        (laurent_in.parent()).polynomial_ring().

    Example:
        sage: R.<x,y,z> = LaurentPolynomialRing(ZZ)
        sage: W1 = matrix_words(GL_generators(3), 9)
        sage: W2 = GL_bounded(3, 1)
        sage: W = tuple(W1.intersection(W2))
        sage: minkowski_0070 = y*z + x + y + z + x^-1 + x^-1*z^-1 + x^-1*y^-1
        sage: minkowski_0021 = x + y + z + x*y^-1 + x^-1 + y^-1*z^-1
        sage: find_mutation(minkowski_0070, minkowski_0021, W, 1)
        [
        [ 0  0  1]         [ 1 -1 -1]
        [ 0 -1 -1]         [ 0  0  1]
        [ 1  1  0], x + 1, [-1  0  0]
        ]
    """
    parent_ring = laurent_out.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError('Incorrect output Laurent polynomial.')

    # Sanity check
    if not (laurent_in in parent_ring):
        raise TypeError('The parent Laurent polynomial rings are different.')

    # Sanity check
    if not (steps in ZZ) or not (steps > -1):
        raise TypeError('Incorrect steps number.')

    # Define the global variables for the wrapper find_mutation_call
    global poly_input
    poly_input = laurent_in
    global poly_out
    poly_out = laurent_out
    global words_input
    words_input = words_tuple

    import multiprocessing
    import gc
    from tqdm import tqdm

    # Parallel computing
    if __name__ == '__main__':
        gc.enable()
        p = multiprocessing.Pool(multiprocessing.cpu_count())

        output = []
        indexes_tuple = Tuples(range(len(words_tuple)), steps + 1)
        for result in p.imap_unordered(find_mutation_call,
                                       tqdm(indexes_tuple), chunksize=1):
            if (result != []):
                output.append(result)
                p.terminate()
                break
        p.close()
        p.join()

    if len(output) > 0:
        return output[0]
    else:
        return []


def mutation_presentation_routine(mutation_list, matrices_list):
    """
    Presents a mutation data in the explicit form assuming that
    the GL-equivalences are fixed.

    Args:
        mutation_list : A list of rational functions.
        matrices_list : A list of monomial matrices for a GL-equivalence.

    Returns:
        A list [matrices_list[0], f, matrices_list[1]], where f is
        an element of the ring PolynomialRing(base_ring, par.gens())
        for par = mutation_list[0].parent(), base_ring = par.base_ring().

    Example:
        sage: R.<x,y,z> = LaurentPolynomialRing(ZZ)
        sage: T = (x + y, x/y, (y*z)/(x + y))
        sage: M1 = matrix(ZZ, [[0, -1, 1],
        ....:                  [1, 0, 0],
        ....:                  [-1, 0, -1]])
        sage: M2 = matrix(ZZ, [[1, -1, 0],
        ....:                  [-1, 0, -1],
        ....:                  [1, 1, 1]])
        sage: mutation_presentation_routine(T, [M1, M2])
        [
        [ 0 -1  1]               [ 1 -1  0]
        [ 1  0  0]               [-1  0 -1]
        [-1  0 -1], x*y^2 + y^2, [ 1  1  1]
        ]
    """
    if (len(matrices_list) != 2):
        raise ValueError('Incorrect matrices data.')
    parent_ring = mutation_list[0].parent()
    laurent_ring = \
        LaurentPolynomialRing(parent_ring.base_ring(), parent_ring.gens())
    frac_field = laurent_ring.fraction_field()
    poly_ring_reduced = \
        PolynomialRing(parent_ring.base_ring(), parent_ring.gens()[:-1])

    # Act on the coordinate functions by the product of matrices
    prod = identity_matrix(ZZ, parent_ring.ngens())
    for i in matrices_list:
        prod *= i

    coord_image = []
    for i in range(matrices_list[0].nrows()):
        coord_image.append(GL_action(laurent_ring.gen(i), prod))

    # Twist the original functions by the obtained monomials
    pre_converted = []
    for i in range(len(mutation_list)):
        pre_converted.append(mutation_list[i] / coord_image[i])

    # Act on the obtained functions by the inverse of the second matrix
    converted = []
    inv = matrices_list[1].inverse()
    for i in range(len(mutation_list)):
        converted.append(GL_action(pre_converted[i], inv))

    # Check if we really obtain the required powers of the same polynomial
    bases = set()
    for i in range(len(mutation_list)):
        ind = matrices_list[0].transpose()[-1][i]
        if (ind != 0):
            temp = power(converted[i], sign(ind))
            if (frac_field(temp).denominator() != 1):
                return []
            else:
                temp = frac_field(temp).numerator()

            if temp not in poly_ring_reduced:
                return []

            c = 1
            for j, k in poly_ring_reduced(temp).factor():
                if (k % (ind*sign(ind)) == 0):
                    c *= power(j, (k / (ind*sign(ind))))
                else:
                    return []
            bases.add(c)
        else:
            if converted[i] != 1:
                return []

    # Return the result
    if len(bases) != 1:
        return []

    return [matrices_list[0],
            poly_ring_reduced(tuple(bases)[0]), matrices_list[1]]


def mutation_presentation_call(args):
    """
    A wrapper for the function mutation_presentation_routine.

    Global variables:
        mutation_input : A list of Laurent polynomial over
                         a polynomial ring over ZZ or QQ.
        words_input : A tuple of Z-invertible matrices of size n
            for n = (laurent_in).parent().ngens().

    Arguments:
        args : A tuple of the elements in range(len(words)).
    """
    return mutation_presentation_routine(
        *tuple([mutation_input, tuple([words_input[i] for i in args])]))


def mutation_presentation(mutation_data, words_tuple):
    """
    Presents a mutation data in the explicit form.

    WARNING! Could be CPU and memory consuming.

    REMARK. Loads tqdm module for the progress bar.

    Arguments:
        mutation_data : A list of rational functions.
        words_tuple : A tuple of possible monomial matrices of GL-equivalences.

    Returns:
        A list [A, f, B], where A, B are invertible matrices over ZZ,
        and f is a polynomial.

    Example:
        sage: R.<x,y,z> = LaurentPolynomialRing(ZZ)
        sage: W1 = matrix_words(GL_generators(3), 9)
        sage: W2 = GL_bounded(3, 1)
        sage: W = tuple(W1.intersection(W2))
        sage: T = tuple([x + y, x/y, (y*z)/(x + y)])
        sage: mutation_presentation(T, W)
        [
        [ 0  0  1]           [ 1 -1  0]
        [ 1  0  0]           [ 1  0  1]
        [-1  1 -1], x*y + y, [-1  1 -1]
        ]
    """
    # Define the global variables for the wrapper mutation_presentation_call
    global mutation_input
    mutation_input = mutation_data
    global words_input
    words_input = words_tuple

    import multiprocessing
    import gc
    from tqdm import tqdm

    # Parallel computing
    if __name__ == '__main__':
        gc.enable()
        p = multiprocessing.Pool(multiprocessing.cpu_count())

        output = []
        indexes_tuple = Tuples(range(len(words_tuple)), 2)
        for result in p.imap_unordered(mutation_presentation_call,
                                       tqdm(indexes_tuple), chunksize=1):
            if result:
                output.append(result)
                p.terminate()
                break
        p.close()
        p.join()

    if len(output) > 0:
        return output[0]
    else:
        return []


# Compute the Newton polytope of a Laurent polynomial

def newton_polytope(laurent_polynomial):
    """
    Returns the Newton polytope of a Laurent polynomial.

    Arguments:
        laurent_polynomial : A Laurent polynomial.

    Returns:
        A lattice polytope.

    Example:
        sage: R.<x,y,z> = LaurentPolynomialRing(ZZ)
        sage: newton_polytope(x + y + x^-1*y^-1 + z + z^-1)
        3-d reflexive polytope in 3-d lattice M
    """
    parent_ring = laurent_polynomial.parent()

    # Sanity check
    if not isinstance(parent_ring, LaurentPolynomialRing_generic):
        raise TypeError("The argument should be a Laurent polynomial.")

    monomials_list = []
    for i in laurent_polynomial.monomials():
        monomials_list.append([i.degree(gen) for gen in parent_ring.gens()])
    newton_polytope = LatticePolytope(monomials_list)

    return newton_polytope


# Compute the quantum period of a Fano CI in Gr(2, N)

def gamma(r):
    """
    An auxiliary function from Theorem 4.1 in arXiv:1409.3729.

    Arguments:
        r : A positive integer.

    Returns:
        A rational number.

    Example:
        sage: gamma(5)
        137/60
    """
    # Sanity check
    if r not in ZZ:
        raise TypeError("The argument should be an integer.")

    if r < 1:
        return 0
    return sum(ZZ.one() / i for i in range(1, r + 1))


def period_coeff(k, multidegree, d):
    """
    Computes the coefficient of the period sequence
    according to Theorem 4.1 in arXiv:1409.3729.

    Arguments:
        k : An integer.
        multidegree : A tuple of integers.
        d : An integer.

    Returns:
        An integer.

    Example:
        sage: period_coeff(3, [1, 2], 3)
        58800
    """
    # Codimension
    lmd = len(multidegree)

    # Index
    d_0 = k + 2 - sum(multidegree)

    temp = 0
    for r in range(d + 1):
        temp += power(binomial(d, r), k + 2) * \
            ((k + 2)*(d - 2*r)*(gamma(r) - gamma(d - r)) + 2)

    coeff = power(-1, d) / 2

    multiple = 1
    for i in range(lmd + 1):
        if (i == 0):
            multiple *= factorial(d * d_0)
        else:
            multiple *= factorial(d * multidegree[i - 1])
    multiple = multiple / power(factorial(d), k + 2)

    return temp * coeff * multiple


def period_series_truncated(affine_ambient_dim, multidegree, r):
    """
    Computes the period sequence according to Theorem 4.1 in arXiv:1409.3729.

    Arguments:
        affine_ambient_dim : An integer.
        multidegree : A tuple of integers.
        r : An integer.

    Returns:
        A polynomial in t.

    Example:
        sage: period_series_truncated(5, [1, 2], 3)
        58800*t^6 + 684*t^4 + 12*t^2 + 1
    """
    temp_ring = PolynomialRing(QQ, 't')
    t = temp_ring.gen()

    k = affine_ambient_dim - 2
    d_0 = affine_ambient_dim - sum(multidegree)

    result = 0
    for i in range(r + 1):
        result += period_coeff(k, multidegree, i) * t**(d_0*i)
    return result


# Try to find a smooth fine regular star triangulation of
# the dual of the Newton polytope by bruteforcing
# the nearest mutations of the Laurent polynomial

def poly_divisors(input_polynomial):
    """
    List all divisors (with multiplicities) of a given polynomial.

    Arguments:
        input_polynomial : A polynomial.

    Returns:
        A polynomial.

    Example:
        sage: S.<x,y,z> = PolynomialRing(QQ)
        sage: poly_divisors(x^2 + y^2 + 2*x*y)
        [1, x + y, x^2 + 2*x*y + y^2]
    """
    poly_parent = input_polynomial.parent()

    output = []

    factor_list = factor(poly_parent(input_polynomial))
    L = [B + 1 for _, B in factor_list]

    for I in mrange(L):
        multiple = 1
        for l in range(len(L)):
            multiple *= factor_list[l][0] ** (I[l])
        output.append(multiple)

    return output


def find_smooth_FRST(laurent_poly, n):
    """
    Bruteforces nearest mutations of the Laurent polynomial
    with respect to n, and tries to fined a smooth fine regular
    star triangulation of the dual of the Newton polytope.

    WARNING: could be heavy on CPU and RAM.

    Arguments:
        laurent_poly : A Laurent polynomial.
        n : An integer.

    Returns:
        None.

    Example:
        sage: R.<a_22,a_11,a_21,a_31> = LaurentPolynomialRing(QQ)
        sage: LG = (a_22 * a_11^-1 + a_22*a_21^-1 + a_31^-1 +
        ....:       a_31 * a_21^-1 + a_22^-1 + a_21 + a_11)
        sage: find_smooth_FRST(LG, 1)  # optional - topcom
        Triangulation...
        Vertices:
        [ 1  1 -2 -2  1  1  1  1]
        [-1  2 -1 -1  2 -1  2 -1]
        [ 2  2 -1 -1 -1 -1 -1 -1]
        [ 1  1  1 -2 -2  1  1 -2]
        Triangulation:
        [0 0 0 0 0]
        [1 1 2 2 3]
        [2 2 3 4 4]
        [3 4 4 5 5]
        [4 6 5 6 7]
        Star origin: 0
        (1, -1, 2, 1)
        Monomial change of variables for a mutation:
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        Mutation factor:
        1
        Mutated Laurent polynomial:
        (a_22^-1*a_11^-1*a_21^-1*a_31^-1) * (a_22*a_11^2*a_21*a_31 +
         a_22*a_11*a_21^2*a_31 + a_22^2*a_11*a_31 + a_22^2*a_21*a_31 +
         a_22*a_11*a_31^2 + a_22*a_11*a_21 + a_11*a_21*a_31)
        ...
    """
    poly_parent = laurent_poly.parent()
    poly_ring = poly_parent.polynomial_ring()
    n_gens = poly_parent.ngens()

    # TOPCOM bug workaround
    Polyhedron(LatticePolytope([[0, 1], [1, 0]]).vertices()).\
        triangulate(engine="topcom")

    matrix_list = [identity_matrix(n_gens)]

    if (n > 0):
        matrix_list = matrix_words(GL_generators(n_gens), n)

    for matrix_word in matrix_list:
        laurent_poly_changed = GL_action(laurent_poly, matrix_word)
        largest_factor = largest_mutation_factor(laurent_poly_changed)
        for mutation_factor in poly_divisors(poly_ring(largest_factor)):
            poly = laurent_polynomial_mutation(laurent_poly_changed,
                                               mutation_factor)
            newton = Polyhedron(newton_polytope(poly).vertices())
            if not newton.is_reflexive():
                continue
            newton_polar = newton.polar()
            pointConf = PointConfiguration(newton_polar.vertices())
            newton_polar_length = len(newton_polar.vertices())
            pointConf_restricted = pointConf.\
                restrict_to_regular_triangulations(regular=True).\
                restrict_to_fine_triangulations(fine=True)
            for star_origin in range(newton_polar_length):
                print("Triangulation...")
                triangulation_list = \
                    list(pointConf_restricted.\
                         restrict_to_star_triangulations(star=star_origin).\
                         triangulations())
                for triangulation in triangulation_list:
                    if ToricVariety(triangulation.fan()).is_smooth():
                        print("Vertices:")
                        print(matrix([newton_polar.Vrepresentation(i)
                                      for i in range(newton_polar_length)]).\
                              transpose())
                        print("Triangulation:")
                        print(matrix(triangulation).transpose())
                        print("Star origin: " + str(star_origin))
                        print(vector(newton_polar.\
                                     Vrepresentation(star_origin)))
                        print("Monomial change of variables for a mutation:")
                        print(matrix_word)
                        print("Mutation factor:")
                        print(factor(mutation_factor))
                        print("Mutated Laurent polynomial:")
                        print(factor(poly))
                        print("\n")


# Compute face polynomials of the Laurent polynomial

def face_polynomial(laurent_poly, face_dim):
    """
    Computes face polynomials of the Laurent polynomial
    for faces of dimension face_dim.

    Arguments:
        laurent_poly : A Laurent polynomial.
        face_dim : An integer.

    Returns:
        None.

    Example:
        sage: R.<a_22,a_11,a_21,a_31> = LaurentPolynomialRing(QQ)
        sage: LG = (a_22 * a_11^-1 + a_22*a_21^-1 + a_31^-1 +
        ....:       a_31 * a_21^-1 + a_22^-1 + a_21 + a_11)
        sage: face_polynomial(LG, 1)
        FACE No. 7
        Vertices: 4
        Lattice points: 4
        Interior lattice points: 0
        Original face polynomial:
        a_11^-1*a_21^-1) * (a_22*a_11*a_21^2 + a_22^2*a_21 +
                            a_22*a_11*a_31 + a_11*a_21)
        Face polynomial in reduced coordinates
        X_0 + X_1 + X_2 + 1
    """
    poly_parent = laurent_poly.parent()

    name_list = []
    for i in range(face_dim):
        name_list.append('X_' + str(i))
        temp_ring = LaurentPolynomialRing(names=name_list, base_ring=QQ)

    newton = newton_polytope(laurent_poly)
    # newton_polyhedron = Polyhedron(newton.vertices(), base_ring=ZZ)

    # newton_integral_points = newton_polyhedron.integral_points()
    # newton_interior_points = newton.interior_points()

    faces = newton.faces(face_dim)

    for L in range(len(faces)):
        face = faces[L]
        print("FACE No. " + str(L))
        face_polyhedron = Polyhedron(face.vertices())
        print("Vertices: " + str(len(face_polyhedron.vertices())))
        print("Lattice points: " + str(len(face_polyhedron.integral_points())))
        print("Interior lattice points: " + str(len(face.interior_points())))
        vector_list = []
        monomial_list = []
        poly_restricted = 0

        for point in face_polyhedron.integral_points():
            monomial = poly_parent.monomial(*list(point))
            coeff = laurent_poly.monomial_coefficient(monomial)
            poly_restricted += coeff * monomial
            vector_list.append(vector(point))

        shift = face_polyhedron.vertices()[0]
        for i in range(len(vector_list)):
            vector_list[i] += -vector(shift)
            monomial_list.append(poly_parent.monomial(*vector_list[i]))

        poly_restricted_shifted = \
            poly_restricted / poly_parent.monomial(*list(shift))
        print("Original face polynomial:")
        print(poly_restricted_shifted)
#        print(factor(poly_restricted_shifted))

        relations = matrix(vector_list).kernel().matrix()
        bases = relations_matrix_bases(relations, False)
        for basis in bases:
            splitting = matrix_splitting(relations, basis)
            non_basis = splitting[0]
            RHS = splitting[2]
            result = 0
            for i in range(1, len(monomial_list) + 1):
                if i in basis:
                    monomial = monomial_list[i - 1]
                    coeff = \
                        poly_restricted_shifted.monomial_coefficient(monomial)
                    index = list(basis).index(i)
                    new_monom = temp_ring.gens()[index]
                    result += new_monom * coeff

                if i in non_basis:
                    monomial = monomial_list[i - 1]
                    coeff = \
                        poly_restricted_shifted.monomial_coefficient(monomial)
                    index = non_basis.index(i)
                    new_monom = temp_ring.monomial(*RHS[index])
                    result += new_monom * coeff
            print("Face polynomial in reduced coordinates")
#            print(factor(result))
            print(result)
        print("\n")


# Compute Newton polytopes of irreducible components of face polynomials

def face_minkowski_polytopes(laurent_poly, face_dim, refined):
    """
    Computes Newton polytopes of irreducible components
    of face polynomials of the Laurent polynomial
    (if refined = true, keep face polynomials separately)

    Arguments:
        laurent_poly : A Laurent polynomial.
        face_dim : An integer.
        refined : A boolean

    Returns:
        A list of the form [P,Q,R], where
        - P is a lattice polytope
        - Q is a Laurent polynomial
        - R is a natural number.

    Example:
        sage: R.<a_22,a_11,a_21,a_31> = LaurentPolynomialRing(QQ)
        sage: LG = (a_22 * a_11^-1 + a_22*a_21^-1 + a_31^-1 +
        ....:       a_31 * a_21^-1 + a_22^-1 + a_21 + a_11)
        sage: face_minkowski_polytopes(LG, 1, True)
        [[3-d lattice polytope in 3-d lattice M,
         X_0^2 + X_0*X_1 + X_0*X_2 + X_0 + X_1, 0],
        [3-d lattice polytope in 3-d lattice M,
         X_0^2 + X_0*X_1 + X_0*X_2 + X_1*X_2 + X_0, 1],
        [3-d lattice polytope in 3-d lattice M,
         X_0^2 + X_0*X_1 + X_0*X_2 + X_1*X_2 + X_0, 2],
        [3-d lattice polytope in 3-d lattice M,
         X_0*X_1 + X_1^2 + X_1*X_2 + X_1 + X_2, 3],
        [3-d lattice polytope in 3-d lattice M, X_0 + X_1 + X_2 + 1, 4],
        [3-d lattice polytope in 3-d lattice M, X_0 + X_1 + X_2 + 1, 5],
        [3-d lattice polytope in 3-d lattice M, X_0 + X_1 + X_2 + 1, 6],
        [3-d lattice polytope in 3-d lattice M, X_0 + X_1 + X_2 + 1, 7]]
    """
    poly_parent = laurent_poly.parent()

    name_list = []
    for i in range(face_dim):
        name_list.append('X_' + str(i))
        temp_ring = LaurentPolynomialRing(names=name_list, base_ring=QQ)

    newton = newton_polytope(laurent_poly)
    # newton_polyhedron = Polyhedron(newton.vertices(), base_ring=ZZ)

    # newton_integral_points = newton_polyhedron.integral_points()
    # newton_interior_points = newton.interior_points()

    faces = newton.faces(face_dim)
    if refined:
        Complete_List = []
    else:
        List = []

    for L in range(len(faces)):
        if refined:
            List = []
        face = faces[L]
        face_polyhedron = Polyhedron(face.vertices())
        vector_list = []
        monomial_list = []
        poly_restricted = 0

        for point in face_polyhedron.integral_points():
            monomial = poly_parent.monomial(*list(point))
            coeff = laurent_poly.monomial_coefficient(monomial)
            poly_restricted += coeff * monomial
            vector_list.append(vector(point))

        shift = face_polyhedron.vertices()[0]
        for i in range(len(vector_list)):
            vector_list[i] += -vector(shift)
            monomial_list.append(poly_parent.monomial(*vector_list[i]))

        poly_restricted_shifted = \
            poly_restricted / poly_parent.monomial(*list(shift))

        relations = matrix(vector_list).kernel().matrix()
        bases = relations_matrix_bases(relations, False)
        for basis in [bases[0]]:
            splitting = matrix_splitting(relations, basis)
            non_basis = splitting[0]
            RHS = splitting[2]
            result = 0
            for i in range(1, len(monomial_list) + 1):
                if i in basis:
                    monomial = monomial_list[i - 1]
                    coeff = \
                        poly_restricted_shifted.monomial_coefficient(monomial)
                    index = list(basis).index(i)
                    new_monom = temp_ring.gens()[index]
                    result += new_monom * coeff

                if i in non_basis:
                    monomial = monomial_list[i - 1]
                    coeff = \
                        poly_restricted_shifted.monomial_coefficient(monomial)
                    index = non_basis.index(i)
                    new_monom = temp_ring.monomial(*RHS[index])
                    result += new_monom * coeff
            for F in factor(temp_ring(result)):
                factor_poly = F[0]
                List.append([newton_polytope(factor_poly), factor_poly, L])
        if refined:
            Complete_List.append(List)
    if refined:
        return Complete_List
    else:
        return List
