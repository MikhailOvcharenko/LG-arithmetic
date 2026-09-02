from sage.all import (ZZ, matrix, Permutation, vector, det,
                      identity_matrix, zero_matrix, LatticePolytope, factor,
                      Polyhedron, binomial, factorial, ToricVariety,
                      lcm, QQ, mrange, PolynomialRing, gcd, divisors,
                      diagonal_matrix, PointConfiguration, xgcd,
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
    permuted_matrix = matrix * permutation
    permuted_matrix_echelon = matrix_reduce(permuted_matrix)

    # Output
    return [complement_list,
            permuted_matrix_echelon[:, :matrix.rank()],
            (-1) * permuted_matrix_echelon[:, matrix.rank():]]


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
            if det(minor_matrix) not in [1, -1]:
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
    if not (matrix_size in ZZ and matrix_size > 0):
        raise ValueError("Incorrect matrix size.")

    if not (matrix_bound in ZZ and matrix_bound > 0):
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
    words_length in a tuple of matrices generators_tuple and its inverses.

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
    monomials_list = [[i.degree(gen) for gen in parent_ring.gens()]
                      for i in laurent_polynomial.monomials()]
    newton_polytope = LatticePolytope(monomials_list).normal_form()

    # Build the list of combinatorially isomorphic polytopes
    try:
        return next(i for i in ReflexivePolytopes(parent_ring.ngens())
                    if newton_polytope == i.normal_form())
    except StopIteration:
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
    if poly_ring_flattened.base_ring() not in (ZZ, QQ):
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
    if mutation_factor not in parent_ring.polynomial_ring():
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
        if i != 0:
            output += laurent_polynomial.coefficient(gen**i) * \
                frac_field(gen)**i * power(mutation_factor, i)
        else:
            output += free_coeff

    # Return to the parent Laurent polynomial ring
    result = output(parent_ring.gens())

    # Sanity check
    if result not in parent_ring:
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
    if len(mutation_data) != 3:
        raise TypeError('Incorrect mutation data.')

    MS = MatrixSpace(ZZ, parent_ngens)
    for i in range(3):
        if i % 2 == 0:
            if mutation_data[i] not in MS:
                raise TypeError('Incorrect mutation data.')
        else:
            if mutation_data[1] not in parent_ring:
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

    if parent_ring.base_ring() not in (ZZ, QQ):
        raise TypeError('Incorrect coefficient ring.')

    # Recursion implementation
    if recursion_depth == len(matrices_list) - 1:
        GL_image = GL_action(laurent_polynomial_input,
                             matrices_list[recursion_depth])
        if GL_image == laurent_polynomial_output:
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
                if output:
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
    if laurent_in not in parent_ring:
        raise TypeError('The parent Laurent polynomial rings are different.')

    # Sanity check
    if not (steps in ZZ and steps > -1):
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
    if len(matrices_list) != 2:
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

    coord_image = [GL_action(laurent_ring.gen(i), prod)
                   for i in range(matrices_list[0].nrows())]

    # Twist the original functions by the obtained monomials
    pre_converted = [mutation_list[i] / coord_image[i]
                     for i in range(len(mutation_list))]

    # Act on the obtained functions by the inverse of the second matrix
    inv = matrices_list[1].inverse()
    converted = [GL_action(pre_converted[i], inv)
                 for i in range(len(mutation_list))]

    # Check if we really obtain the required powers of the same polynomial
    bases = set()
    for i in range(len(mutation_list)):
        ind = matrices_list[0].transpose()[-1][i]
        if ind != 0:
            temp = power(converted[i], sign(ind))
            if frac_field(temp).denominator() != 1:
                return []
            else:
                temp = frac_field(temp).numerator()

            if temp not in poly_ring_reduced:
                return []

            c = 1
            for j, k in poly_ring_reduced(temp).factor():
                if k % (ind * sign(ind)) == 0:
                    c *= power(j, (k / (ind * sign(ind))))
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

    monomials_list = [[i.degree(gen) for gen in parent_ring.gens()]
                      for i in laurent_polynomial.monomials()]
    return LatticePolytope(monomials_list)


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
        if i == 0:
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
    with respect to n, and tries to find a smooth fine regular
    star triangulation of the dual of the Newton polytope.

    The triangulation has the star center at the origin,
    as required for maximal partial crepant resolutions.

    We search only for FRS triangulations supported on the set
    [origin + vertices] instead of the whole integral points set
    in order to avoid huge computations.
    In other words, this is only a sufficient condition.

    WARNING: could be heavy on CPU and RAM.

    Arguments:
        laurent_poly : A Laurent polynomial.
        n : An integer.

    Returns:
        None.

    Example:
        sage: R.<a,b,c,d> = LaurentPolynomialRing(QQ)
        sage: LG = (a+b+c+d+1)^3*(1/b/d+1/a/d+1/a/c);
        sage: find_smooth_FRST(LG, 1)
        Vertices:
        [(0, -1, 0, 0, 0, 0, 1, 1, 1),
         (0, -1, 0, 0, 0, 1, 0, 0, 1),
         (0, 0, -1, 0, 1, 0, -1, 0, 0),
         (0, 1, 0, -1, 0, 0, 0, 0, 0)]
        Triangulation:
        [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
         (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3),
         (2, 2, 2, 3, 3, 3, 4, 4, 5, 6, 3, 4, 4, 5, 6),
         (3, 3, 5, 4, 4, 6, 5, 7, 6, 7, 5, 5, 7, 6, 7),
         (5, 6, 6, 5, 7, 7, 8, 8, 8, 8, 6, 8, 8, 8, 8)]
        Star origin: 0
        M(0, 0, 0, 0)
        Monomial change of variables for a mutation:
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        Mutation factor:
        1
        Mutated Laurent polynomial:
        (a^-1*b^-1*c^-1*d^-2) * (c + 1) * (a*c + a + b) * (a*d + b*d + d + 1)^3
    """
    poly_parent = laurent_poly.parent()
    poly_ring = poly_parent.polynomial_ring()
    n_gens = poly_parent.ngens()

    # TOPCOM bug workaround
    Polyhedron(LatticePolytope([[0, 1], [1, 0]]).vertices()).\
        triangulate(engine="topcom")

    matrix_list = [identity_matrix(n_gens)]

    if n > 0:
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
            points = [vector(ZZ, p) for p in newton_polar.integral_points()];
            origin = LatticePolytope(points).interior_points()[0];
            origin_index = points.index(origin);
            points_red = [vector(ZZ, origin)] + [vector(ZZ, p) for p in newton_polar.vertices()];
            origin_red_index = points_red.index(origin);

            # Remove _red below to look for a general FRS triangulation
            # WARNING: TOPCOM memory consumption is usually huge
            
            pointConf = PointConfiguration(points_red)
            pointConf_restricted = pointConf.\
                restrict_to_regular_triangulations(regular=True).\
                restrict_to_fine_triangulations(fine=True).\
                restrict_to_star_triangulations(star=origin_red_index);
            for triangulation in pointConf_restricted.triangulations():
                if ToricVariety(triangulation.fan()).is_smooth():
                    print("Vertices:")
                    print(list(matrix(points_red).transpose()))
                    print("Triangulation:")
                    print(list(matrix(triangulation).transpose()))
                    print("Star origin: " + str(origin_red_index))
                    print(points_red[origin_red_index]);
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
        sage: face_polynomial(LG, 3)
        FACE No. 7
        Vertices: 4
        Lattice points: 4
        Interior lattice points: 0
        Original face polynomial:
        (a_11^-1*a_21^-1) * (a_22*a_11*a_21^2 + a_22^2*a_21 +
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
            print(factor(result))
#            print(result)
        print("\n")


# Compute Newton polytopes of irreducible components of face polynomials

def face_minkowski_polytopes(laurent_poly, face_dim, refined):
    """
    Computes Newton polytopes of irreducible components
    of face polynomials of the Laurent polynomial
    (if refined=true, keep face polynomials separately)

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
        sage: face_minkowski_polytopes(LG, 3, True)
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


# Present a Laurent polynomial whose Newton polytope is of lattice width one
# in the form F(x_0, ..., x_{n - 1}) * x_n + G(x_0, ..., x_{n - 1})
# (and raise an exception if this fails)

def _width_one_exponent_tuple(exponent, n):
    """
    Normalize an exponent key returned by monomial_coefficients().

    SageMath uses an integer key for some univariate Laurent polynomial
    implementations and a tuple-like key for multivariate implementations.
    """
    if n == 1:
        try:
            return (ZZ(exponent),)
        except (TypeError, ValueError):
            return (ZZ(exponent[0]),)
    return tuple(ZZ(a) for a in exponent)


def _width_one_support_exponents(f):
    """
    Return the exponent vectors of the nonzero monomials of f.

    In some SageMath versions f.monomials() method calls
    monomial_coefficients(copy=...), whereas the Laurent polynomial
    implementation does not accept the 'copy' keyword.
    """
    n = f.parent().ngens()
    return [
        vector(ZZ, _width_one_exponent_tuple(exponent, n))
        for exponent in f.monomial_coefficients()
    ]


def _canonical_direction_sign(u):
    """
    Choose one of u and -u: the first nonzero entry is positive.
    """
    for a in u:
        if a > 0:
            return vector(ZZ, u)
        if a < 0:
            return -vector(ZZ, u)
    return vector(ZZ, u)


def lattice_width_one_directions(f):
    """
    Return all lattice-width-one covectors of Newt(f) modulo sign.

    The output consists of primitive vectors u in ZZ^n such that
    max_{a in supp(f)} <a,u> - min_{a in supp(f)} <a,u> = 1.

    The Newton polytope is required to be full-dimensional.

    Algorithm:
        For P = Newt(f), put D = P-P. Then width_P(u) = h_D(u),
        where h_D is the support function. 

        Hence the integral covectors of width at most one are exactly the
        lattice points of the bounded rational polytope
        D^* = {u : <p-q,u> <= 1 for all vertices p,q of P}.

        We enumerate these lattice points and discard zero.

    Example:
        sage: R.<x,y,z> = LaurentPolynomialRing(QQ)
        sage: f = x^2 + x*y^2 + y*z^2 + z
        sage: lattice_width_one_directions(f)
        [(1, 0, 1)]
    """
    R = f.parent()

    if not isinstance(R, LaurentPolynomialRing_generic):
        raise TypeError("f must belong to a Laurent polynomial ring.")
    if f.is_zero():
        raise ValueError("The zero Laurent polynomial has no Newton polytope.")

    n = R.ngens()
    support = _width_one_support_exponents(f)
    P = Polyhedron(vertices=[list(a) for a in support], base_ring=QQ)

    if P.dim() != n:
        raise ValueError(
            "The Newton polytope is not full-dimensional."
        )

    vertices = [vector(ZZ, v) for v in P.vertices_list()]

    # Sage inequalities have the form b + a_1*u_1 + ... + a_n*u_n >= 0.
    # For every unordered pair p,q we impose both
    #     1 - <p-q,u> >= 0  and  1 + <p-q,u> >= 0.
    ieqs = []
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            difference = vertices[i] - vertices[j]
            ieqs.append([ZZ(1)] + list(-difference))
            ieqs.append([ZZ(1)] + list(difference))

    dual_difference_body = Polyhedron(ieqs=ieqs, base_ring=QQ)
    if not dual_difference_body.is_compact():
        raise ArithmeticError(
            "Internal error: the polar of the difference body is not compact."
        )

    directions = set()
    for point in dual_difference_body.integral_points():
        u = vector(ZZ, point)
        if all(a == 0 for a in u):
            continue

        levels = [a.dot_product(u) for a in support]
        if max(levels) - min(levels) == 1:
            u = _canonical_direction_sign(u)
            directions.add(tuple(ZZ(a) for a in u))

    if not directions:
        raise ValueError("The Newton polytope does not have lattice width one.")

    ordered = sorted(
        directions,
        key=lambda u: (sum(a*a for a in u), u)
    )
    return [vector(ZZ, u) for u in ordered]


def _unimodular_matrix_with_last_column(u):
    """
    Complete a primitive vector u to a matrix in GL(n,ZZ),
    with u as its last column.

    This uses elementary extended-gcd row operations, so it does not
    depend on optional normal-form packages.
    """
    u = vector(ZZ, u)
    n = len(u)

    if n == 0:
        raise ValueError("The direction vector must be nonempty.")
    if gcd(list(u)) not in (ZZ(1), ZZ(-1)):
        raise ValueError("The direction vector must be primitive.")

    # Construct U in GL(n,ZZ) with U*u = e_1.
    U = identity_matrix(ZZ, n)
    v = vector(ZZ, u)

    for i in range(1, n):
        a = ZZ(v[0])
        b = ZZ(v[i])
        if b == 0:
            continue

        g, s, t = xgcd(a, b)  # s*a + t*b = g, with g >= 0
        E = identity_matrix(ZZ, n)
        E[0, 0] = s
        E[0, i] = t
        E[i, 0] = -b // g
        E[i, i] = a // g

        U = E * U
        v = E * v

    if v[0] == -1:
        E = identity_matrix(ZZ, n)
        E[0, 0] = -1
        U = E * U
        v = E * v

    if v != vector(ZZ, [1] + [0]*(n - 1)):
        raise ArithmeticError("Failed to complete the primitive vector.")

    # C = U^{-1} has u as its first column. We cyclically move that column
    # to the end. The resulting M is still unimodular.
    C = matrix(ZZ, U.inverse())
    M = matrix(
        ZZ, n, n,
        lambda i, j: C[i, (j + 1) % n]
    )

    if abs(M.det()) != 1 or M.column(n - 1) != u:
        raise ArithmeticError("Failed to construct the monomial matrix.")

    return M


def _validate_width_direction(support, u):
    """
    Validate a custom primitive width-one covector.
    """
    u = vector(ZZ, u)
    if len(u) != len(support[0]):
        raise ValueError("The direction has the wrong dimension.")
    if gcd(list(u)) not in (ZZ(1), ZZ(-1)):
        raise ValueError("The direction must be primitive.")

    levels = [a.dot_product(u) for a in support]
    if max(levels) - min(levels) != 1:
        raise ValueError("The custom direction does not compute width one.")
    return _canonical_direction_sign(u)


def _split_linear_last_variable(h):
    """
    Given h whose last exponents are 0 and 1, return F,G with h = F * x_n + G.
    """
    R = h.parent()
    n = R.ngens()
    x_n = R.gen(n - 1)

    F = R.zero()
    G = R.zero()

    for exponent, coefficient in h.monomial_coefficients().items():
        exponent = _width_one_exponent_tuple(exponent, n)

        last_degree = exponent[-1]
        base_exponent = exponent[:-1] + (ZZ(0),)
        base_monomial = R.monomial(*base_exponent)

        if last_degree == 0:
            G += coefficient * base_monomial
        elif last_degree == 1:
            F += coefficient * base_monomial
        else:
            raise ArithmeticError(
                "The normalized polynomial is not linear in the last variable."
            )

    if h != F*x_n + G:
        raise ArithmeticError("Failed to split the normalized polynomial.")

    # Explicitly verify that F and G do not involve x_n.
    for coefficient_part in (F, G):
        for exponent in _width_one_support_exponents(coefficient_part):
            if exponent[-1] != 0:
                raise ArithmeticError("F or G still depends on the last variable.")

    return F, G


def linearize_width_one(f, direction=None):
    """
    Transform a Laurent polynomial of lattice width one into
    F(x_1,...,x_{n-1})*x_n + G(x_1,...,x_{n-1}).

    Arguments:
        f : a nonzero Laurent polynomial with full-dimensional
            lattice-width-one Newton polytope.
        direction : optional primitive width-one covector.
                    If omitted, one is found automatically.

    Returns:
        A dictionary with keys:
        - polynomial : the normalized polynomial F*x_n + G;
        - F, G : Laurent polynomials independent of x_n;
        - raw_GL_image : GL_action(f, M);
        - matrix : M in GL(n,ZZ), whose last column is the width covector;
        - width_direction : the primitive covector u;
        - minimum_level : m = min <a,u>;
        - laurent_unit : x_n^(-m), so that
          polynomial = laurent_unit * raw_GL_image.

    Conventions:
        GL_action(f,M) sends an exponent row vector a to a*M. Therefore its
        last coordinate is <a,u> when the last column of M is u.

        A GL(n,ZZ) substitution alone gives last exponents {m,m+1}. Multiplying
        by the Laurent unit x_n^(-m) translates them to {0,1}. Multiplication
        by a nonzero Laurent monomial does not change the toric hypersurface.

    Example:
        sage: R.<X_1, X_2, X_3> = LaurentPolynomialRing(QQ)
        sage: f = X_1^2 + X_1*X_2^2 + X_2*X_3^2 + X_3
        sage: output = linearize_width_one(f)
        sage: output["width_direction"]
        (1, 0, 1)
        sage: output["polynomial"] == output["F"]*X_3 + output["G"]
        True
    """
    R = f.parent()

    if not isinstance(R, LaurentPolynomialRing_generic):
        raise TypeError("f must belong to a Laurent polynomial ring.")
    if f.is_zero():
        raise ValueError("f must be nonzero.")

    support = _width_one_support_exponents(f)
    n = R.ngens()

    if direction is None:
        u = lattice_width_one_directions(f)[0]
    else:
        # Also performs the ambient-dimension check through the polytope.
        P = Polyhedron(vertices=[list(a) for a in support], base_ring=QQ)
        if P.dim() != n:
            raise ValueError(
                "The Newton polytope is not full-dimensional."
            )
        u = _validate_width_direction(support, direction)

    M = _unimodular_matrix_with_last_column(u)
    raw = GL_action(f, M)

    raw_support = _width_one_support_exponents(raw)
    last_degrees = [a[-1] for a in raw_support]
    minimum_level = min(last_degrees)
    maximum_level = max(last_degrees)

    if maximum_level - minimum_level != 1:
        raise ArithmeticError(
            "GL_action did not put the width direction last."
        )

    x_n = R.gen(n - 1)
    laurent_unit = x_n**(-minimum_level)

    # GL_action constructs a formal Laurent ring from the generators.  For
    # one generator Sage may switch between its univariate and one-variable
    # multivariate Laurent implementations, which need not coerce negative
    # valuations directly into one another.  Perform the shift in raw's own
    # parent, then convert the normalized (nonnegative in x_n) result to R.
    raw_x_n = raw.parent().gen(n - 1)
    normalized = R(raw_x_n**(-minimum_level) * raw)
    F, G = _split_linear_last_variable(normalized)

    return {
        "polynomial": normalized,
        "F": F,
        "G": G,
        "raw_GL_image": raw,
        "matrix": M,
        "width_direction": u,
        "minimum_level": minimum_level,
        "laurent_unit": laurent_unit,
    }


# Check weak and nested non-degeneracy of a Laurent polynomial
# Weak non-degeneracy condition is relaxed in comparison with arXiv:2307.15607

def _primitive_integer_vector(v):
    """
    Returns the primitive integral generator of a rational ray.

    The orientation of the input vector is preserved.

    Arguments:
        v : A nonzero vector over ZZ or QQ.

    Returns:
        A primitive vector over ZZ spanning the same oriented ray as v.

    Example:
        sage: _primitive_integer_vector(vector(QQ, [2/3, 4/3]))
        (1, 2)
    """
    v = vector(QQ, v)
    if all(a == 0 for a in v):
        raise ValueError("The zero vector does not span a ray.")

    denominator = lcm([a.denominator() for a in v])
    w = vector(ZZ, [ZZ(denominator*a) for a in v])
    common_divisor = gcd(list(w))
    if common_divisor < 0:
        common_divisor = -common_divisor
    return vector(ZZ, [a // common_divisor for a in w])


def _face_reduced_data(laurent_poly, face):
    """
    Computes a face polynomial in coordinates of the saturated face lattice.

    Let delta be the given face.  One vertex of delta is chosen as the
    origin, and the lattice parallel to delta is replaced by its saturation
    in the ambient exponent lattice.  The face polynomial is then expressed,
    up to multiplication by a Laurent monomial, in a Laurent polynomial ring
    whose exponent lattice is this saturated lattice.

    This differs from face_polynomial(): the saturated face lattice
    is needed in order to test relative unimodularity correctly.

    Arguments:
        laurent_poly : A Laurent polynomial over QQ.
        face : A positive-dimensional face of its Newton polytope.

    Returns:
        A dictionary with keys:
        - polynomial : the face polynomial in saturated reduced coordinates,
          up to multiplication by a Laurent monomial;
        - polyhedron : the face in the reduced lattice coordinates;
        - lattice_basis : a matrix whose rows form a basis of the saturated
          lattice parallel to face;
        - shift : the ambient lattice point chosen as the origin.
    """
    poly_parent = laurent_poly.parent()
    if not isinstance(poly_parent, LaurentPolynomialRing_generic):
        raise TypeError("The first argument should be a Laurent polynomial.")

    vertices = [vector(ZZ, v) for v in face.vertices()]
    if len(vertices) == 0:
        raise ValueError("The face should be nonempty.")

    shift = vertices[0]
    differences = matrix(ZZ, [list(v - shift) for v in vertices])
    lattice_basis = differences.saturation()
    face_dim = lattice_basis.rank()

    if face_dim < 1:
        raise ValueError("The face should be positive-dimensional.")

    # saturation() returns a full-rank matrix, but we remove zero rows
    # defensively in case this changes in Sage.
    lattice_basis = matrix(
        ZZ, [list(row) for row in lattice_basis.rows()
             if any(a != 0 for a in row)])

    if lattice_basis.nrows() != face_dim:
        raise ArithmeticError("Failed to compute the saturated face lattice.")

    name_list = ['X_' + str(i) for i in range(face_dim)]
    temp_ring = LaurentPolynomialRing(names=name_list, base_ring=QQ)

    def reduced_coordinates(point):
        """
        Returns coordinates of a face lattice point in the saturated basis.
        """
        difference = vector(QQ, vector(ZZ, point) - shift)
        coordinates = matrix(QQ, lattice_basis.transpose()).solve_right(
            difference)
        if any(a.denominator() != 1 for a in coordinates):
            raise ArithmeticError(
                "A face lattice point has non-integral reduced coordinates.")
        return vector(ZZ, [ZZ(a) for a in coordinates])

    reduced_vertices = [reduced_coordinates(v) for v in vertices]
    face_polyhedron = Polyhedron(vertices=reduced_vertices, base_ring=QQ)
    ambient_face_polyhedron = Polyhedron(vertices=vertices, base_ring=QQ)

    result = temp_ring.zero()
    n = poly_parent.ngens()
    for exponent, coefficient in laurent_poly.monomial_coefficients().items():
        exponent = _width_one_exponent_tuple(exponent, n)
        point = vector(ZZ, exponent)
        if not ambient_face_polyhedron.contains(point):
            continue
        coordinates = reduced_coordinates(point)
        result += QQ(coefficient) * temp_ring.monomial(*list(coordinates))

    if result.is_zero():
        raise ArithmeticError("Failed to construct the face polynomial.")

    return {
        "polynomial": result,
        "polyhedron": face_polyhedron,
        "lattice_basis": lattice_basis,
        "shift": shift,
    }


def _relative_normal_rays(face_polyhedron, subface):
    """
    Computes primitive rays of a relative normal cone.

    The ambient polyhedron is assumed to be full-dimensional in its lattice.
    Sage inequalities have the form b + <A,x> >= 0.  Therefore the vector A
    for a facet containing subface points towards increasing order along the
    corresponding toric normal coordinate and has the required orientation.

    Arguments:
        face_polyhedron : A full-dimensional lattice polyhedron.
        subface : A proper face of face_polyhedron.

    Returns:
        A list of primitive vectors over ZZ generating the extremal rays of
        the relative normal cone N_{subface/face_polyhedron}.
    """
    rays = []
    subface_vertices = list(subface.vertices())

    for inequality in face_polyhedron.inequality_generator():
        if all(inequality.eval(v) == 0 for v in subface_vertices):
            ray = _primitive_integer_vector(inequality.A())
            if ray not in rays:
                rays.append(ray)

    return rays


def _relative_normal_cone_data(face_polyhedron, subface):
    """
    Computes lattice data of a relative normal cone.

    Unimodularity is tested in the saturated lattice generated by the rays,
    not by taking a determinant in the original ambient lattice.

    Arguments:
        face_polyhedron : A full-dimensional lattice polyhedron.
        subface : A proper face of face_polyhedron.

    Returns:
        A dictionary with keys:
        - dimension : the dimension of the relative normal cone;
        - rays : its primitive ray generators;
        - simplicial : whether the rays form a simplicial cone;
        - lattice_index : the index of the ray lattice in its saturation if
          the cone is simplicial, and None otherwise;
        - unimodular : True precisely when the cone is simplicial and has
          lattice index one.
    """
    relative_dim = face_polyhedron.dim() - subface.dim()
    rays = _relative_normal_rays(face_polyhedron, subface)

    if len(rays) == 0:
        raise ArithmeticError("Failed to compute the relative normal cone.")

    ray_matrix = matrix(ZZ, [list(v) for v in rays])
    simplicial = (len(rays) == relative_dim and
                  ray_matrix.rank() == relative_dim)

    if simplicial:
        lattice_index = ZZ(ray_matrix.index_in_saturation())
    else:
        lattice_index = None

    return {
        "dimension": relative_dim,
        "rays": rays,
        "simplicial": simplicial,
        "lattice_index": lattice_index,
        "unimodular": simplicial and lattice_index == 1,
    }


def _unimodular_normal_matrix(rays, ambient_dim):
    """
    Completes primitive generators of a unimodular cone to a lattice basis.

    The rows of the returned matrix give monomial exponents in the new toric
    coordinates.  Its first rows are exactly the supplied normal rays; the
    remaining rows are tangent directions along the corresponding stratum.

    Arguments:
        rays : Linearly independent primitive vectors generating a
            unimodular cone.
        ambient_dim : The rank of the ambient lattice.

    Returns:
        A matrix U in GL(ambient_dim, ZZ) whose first len(rays) rows are rays.

    Example:
        sage: rays = [vector(ZZ, [1, 0]), vector(ZZ, [0, 1])]
        sage: _unimodular_normal_matrix(rays, 2)
        [1 0]
        [0 1]
    """
    ray_matrix = matrix(ZZ, [list(v) for v in rays])
    c = ray_matrix.nrows()

    if ray_matrix.ncols() != ambient_dim:
        raise ValueError("The normal rays have the wrong ambient dimension.")
    if c > ambient_dim or ray_matrix.rank() != c:
        raise ValueError("The normal rays should be linearly independent.")
    if ray_matrix.index_in_saturation() != 1:
        raise ValueError("The relative normal cone is not unimodular.")

    # Since the rows of ray_matrix span a primitive sublattice, the columns
    # of ray_matrix^t extend to a basis.  Hermite reduction gives
    # U*ray_matrix^t = [I_c; 0].
    hermite, transformation = ray_matrix.transpose().hermite_form(
        transformation=True, include_zero_rows=True)
    expected = zero_matrix(ZZ, ambient_dim, c)
    for i in range(c):
        expected[i, i] = 1

    if hermite != expected:
        raise ArithmeticError("Failed to reduce a unimodular normal cone.")

    completion = matrix(ZZ, transformation.inverse()).transpose()

    if abs(completion.det()) != 1:
        raise ArithmeticError("Failed to construct unimodular coordinates.")
    if completion[:c, :] != ray_matrix:
        raise ArithmeticError("Normal rays were not preserved by completion.")

    return completion


def _polynomial_in_normal_chart(laurent_poly, normal_matrix, normal_dim,
                                poly_ring):
    """
    Extends a Laurent polynomial to a toric normal chart.

    The chart is A^c x (G_m)^(d-c), where c = normal_dim.  The first c rows
    of normal_matrix are primitive normal rays and the remaining rows are
    tangent lattice directions.  Exponents are transformed by normal_matrix.
    For each normal coordinate the minimum valuation is subtracted, giving
    the regular extension across its boundary divisor.  Minimum exponents in
    tangent coordinates are also subtracted; this only multiplies the
    polynomial by a Laurent unit on the torus factor and clears denominators.

    Arguments:
        laurent_poly : A Laurent polynomial in d variables over QQ.
        normal_matrix : A matrix in GL(d, ZZ) defining the chart coordinates.
        normal_dim : The number c of normal coordinates.
        poly_ring : A polynomial ring over QQ in d variables.  Its first c
            generators are interpreted as normal variables.

    Returns:
        A polynomial in poly_ring defining the same reduced hypersurface in
        A^c x (G_m)^(d-c), up to multiplication by a torus unit.
    """
    d = laurent_poly.parent().ngens()
    if normal_matrix.nrows() != d or normal_matrix.ncols() != d:
        raise ValueError("Incorrect normal-coordinate matrix.")
    if poly_ring.ngens() != d:
        raise ValueError("Incorrect polynomial ring.")
    if normal_dim not in range(d + 1):
        raise ValueError("Incorrect number of normal variables.")

    transformed_terms = []
    for exponent, coefficient in laurent_poly.monomial_coefficients().items():
        exponent = _width_one_exponent_tuple(exponent, d)
        transformed_exponent = normal_matrix * vector(ZZ, exponent)
        transformed_terms.append((transformed_exponent, QQ(coefficient)))

    coordinate_minima = []
    for j in range(d):
        coordinate_minima.append(
            min(exponent[j] for exponent, _ in transformed_terms))

    result = poly_ring.zero()
    for exponent, coefficient in transformed_terms:
        shifted = [ZZ(exponent[j] - coordinate_minima[j])
                   for j in range(d)]
        result += coefficient * poly_ring.monomial(*shifted)

    return result


def _saturate_by_variables(input_ideal, variables):
    """
    Saturates an ideal by a collection of coordinate variables.

    This removes components supported on coordinate hyperplanes that are
    absent from the torus factor of a normal chart.

    Arguments:
        input_ideal : An ideal in a polynomial ring.
        variables : A list of generators of the same polynomial ring.

    Returns:
        The saturation of input_ideal by the product of variables.  If the
        list is empty, input_ideal itself is returned.
    """
    if len(variables) == 0:
        return input_ideal

    multiple = input_ideal.ring().one()
    for variable in variables:
        multiple *= variable
    return input_ideal.saturation(input_ideal.ring().ideal(multiple))[0]


def _ideal_is_unit(input_ideal):
    """
    Tests whether an ideal is the unit ideal.

    Arguments:
        input_ideal : An ideal in a polynomial ring.

    Returns:
        A boolean.
    """
    return input_ideal.ring().one() in input_ideal


def _irreducible_component_smooth_along_stratum(
        component, normal_variables, torus_variables,
        component_is_torus_saturated=False):
    """
    Tests one irreducible component for smoothness along a toric stratum.

    Arguments:
        component : A prime ideal in a polynomial ring over QQ.
        normal_variables : The coordinate variables vanishing on the stratum.
        torus_variables : The invertible coordinate variables of the chart.
        component_is_torus_saturated : A boolean.  If True, component is
            already known to be saturated by the torus variables, so the
            redundant saturation step is skipped.

    Returns:
        A boolean.  Components removed by torus saturation count as empty and
        hence pass the test.
    """
    if not component_is_torus_saturated:
        component = _saturate_by_variables(component, torus_variables)
        if _ideal_is_unit(component):
            return True

    poly_ring = component.ring()
    ambient_dim = poly_ring.ngens()
    codimension = ambient_dim - component.dimension()

    jacobian_list = []
    for poly in component.gens():
        gradient = [poly.derivative(g) for g in poly_ring.gens()]
        jacobian_list.append(gradient)
    jacobian_matrix = matrix(poly_ring, jacobian_list)

    singular_list = list(component.gens()) + list(normal_variables)
    if codimension > 0:
        singular_list += jacobian_matrix.minors(codimension)

    singular_ideal = poly_ring.ideal(singular_list)
    singular_ideal = _saturate_by_variables(
        singular_ideal, torus_variables)
    return _ideal_is_unit(singular_ideal)


def _component_collection_smooth_along_stratum(
        polynomial_list, normal_dim, prime_chart_seeds=False):
    """
    Tests the recursive componentwise smoothness condition at one stratum.

    The ambient chart is A^c x (G_m)^(d-c), where c = normal_dim and the first
    c variables are normal.  polynomial_list is the list of regular extensions
    of the distinct irreducible factors of one reduced face polynomial.  The
    case c = 0 is the relaxed weak non-degeneracy test on the open torus of
    the face under consideration.

    The collection of strata is generated globally from all these factor
    closures.  Each factor closure is first reduced and decomposed into minimal
    irreducible components.  The collection is then closed under reduced
    pairwise intersection: every new intersection is decomposed into its
    minimal primes, and every new irreducible component is added to the same
    collection.  In the ordinary path an explicit radical is taken first; in
    the optimized chart-seed path it is omitted because an ideal and its
    radical have exactly the same minimal primes.  Iteration therefore
    produces the irreducible components
    of all successive finite reduced intersections.  Every component obtained
    at any stage is required to be smooth along z_1 = ... = z_c = 0.

    This is equivalent to first considering the reduced common zero locus for
    every nonempty subset of the face factors and then closing the union of all
    of their irreducible components under reduced intersections.  Starting
    from the singleton factor closures is enough, since the components of a
    common zero locus are produced when those singleton components are
    intersected and decomposed.

    No transversality or expected-codimension condition is imposed.  Distinct
    smooth components are allowed to meet, provided that every irreducible
    component arising from their successive reduced intersections is smooth
    along the stratum.

    Arguments:
        polynomial_list : A nonempty list of polynomials in one polynomial
            ring over QQ, namely the chart extensions of all distinct
            irreducible factors of a reduced face polynomial.
        normal_dim : The number of initial variables defining the normal
            coordinate hyperplanes.
        prime_chart_seeds : A boolean optimization flag.  If True,
            polynomial_list is known to consist of chart extensions of
            irreducible Laurent factors after a unimodular monomial change of
            coordinates and subtraction of coordinate minima.  Each principal
            ideal is then already prime and saturated with respect to the
            torus variables, so the initial saturation, radical, and
            minimal-prime decomposition are skipped.  Prime components of
            later torus-saturated intersections are likewise not re-saturated.

    Returns:
        A triple (smooth, reason, witness), where smooth is a boolean.  On
        success reason and witness are None.  On failure reason is
        "singular irreducible component of an iterated reduced intersection"
        and witness is a dictionary containing generators of the offending
        prime ideal and the factor indexes whose closures generated it.
    """
    if len(polynomial_list) == 0:
        raise ValueError("At least one polynomial is required.")

    poly_ring = polynomial_list[0].parent()
    ambient_dim = poly_ring.ngens()
    if normal_dim not in range(ambient_dim + 1):
        raise ValueError("Incorrect number of normal variables.")

    normal_variables = list(poly_ring.gens()[:normal_dim])
    torus_variables = list(poly_ring.gens()[normal_dim:])

    components = []
    factor_supports = []

    def add_component(component, support):
        """Adds and checks one irreducible component if it is new."""
        if not prime_chart_seeds:
            component = _saturate_by_variables(component, torus_variables)
            if _ideal_is_unit(component):
                return (True, None)

        support = frozenset(support)
        for known_number in range(len(components)):
            if component == components[known_number]:
                # Keep the smallest known factor support only for diagnostics.
                if len(support) < len(factor_supports[known_number]):
                    factor_supports[known_number] = support
                return (True, known_number)

        if not _irreducible_component_smooth_along_stratum(
                component, normal_variables, torus_variables,
                component_is_torus_saturated=prime_chart_seeds):
            return (
                False,
                {
                    "component_generators": list(component.gens()),
                    "factor_subset": sorted(support),
                },
            )

        components.append(component)
        factor_supports.append(support)
        return (True, len(components) - 1)

    # Seed one common collection with the irreducible components of every
    # singleton factor closure. The collection is shared by all factors
    # rather than being restarted for each factor subset.
    for factor_number in range(len(polynomial_list)):
        input_ideal = poly_ring.ideal([polynomial_list[factor_number]])

        if prime_chart_seeds:
            # In the optimized chart-seed path polynomial_list comes
            # from distinct irreducible Laurent factors.  A unimodular
            # monomial coordinate change and subtraction of coordinate minima
            # preserve irreducibility and remove every common torus-variable
            # factor, so this principal ideal is prime and torus-saturated.
            seed_components = [input_ideal]
        else:
            input_ideal = _saturate_by_variables(
                input_ideal, torus_variables)
            input_ideal = input_ideal.radical()

            if _ideal_is_unit(input_ideal):
                continue

            seed_components = input_ideal.minimal_associated_primes()

        for component in seed_components:
            added, witness = add_component(component, [factor_number])
            if not added:
                return (
                    False,
                    "singular irreducible component of an iterated reduced "
                    "intersection",
                    witness,
                )

    # Close the single global collection under reduced intersections.  When a
    # new component is appended, the while loop automatically intersects it
    # with every component discovered earlier.  Hence all finite successive
    # intersections are eventually generated.
    component_number = 0
    while component_number < len(components):
        for other_number in range(component_number):
            intersection_ideal = poly_ring.ideal(
                list(components[component_number].gens()) +
                list(components[other_number].gens()))
            intersection_ideal = _saturate_by_variables(
                intersection_ideal, torus_variables)
            if not prime_chart_seeds:
                intersection_ideal = intersection_ideal.radical()

            if _ideal_is_unit(intersection_ideal):
                continue

            support = (factor_supports[component_number] |
                       factor_supports[other_number])
            for component in intersection_ideal.minimal_associated_primes():
                added, witness = add_component(component, support)
                if not added:
                    return (
                        False,
                        "singular irreducible component of an iterated "
                        "reduced intersection",
                        witness,
                    )

        component_number += 1

    return (True, None, None)


def _face_factor_list(face_polynomial):
    """
    Returns the distinct nonconstant irreducible factors of a face polynomial.

    Multiplicities are deliberately discarded, in accordance with passing to
    the radical of the face polynomial.

    Arguments:
        face_polynomial : A polynomial or Laurent polynomial over QQ.

    Returns:
        A list of distinct nonconstant irreducible factors.
    """
    output = []
    for factor_poly, _ in factor(face_polynomial):
        if factor_poly.is_constant():
            continue
        output.append(factor_poly)
    return output


def _factor_intersections_smooth_in_torus(factor_list):
    """
    Tests relaxed weak non-degeneracy for Laurent factors in a face torus.

    The distinct factor hypersurfaces are put in one common collection.  The
    collection is closed recursively under reduced intersection and
    irreducible decomposition, and every irreducible component obtained at
    any stage is required to be smooth in the algebraic torus.  Reducible
    intersections are therefore allowed, and no transversality or
    expected-codimension condition is imposed.

    This is the normal_dim = 0 instance of the same componentwise recursive
    condition used by the nested checker.  Since factor_list already consists
    of distinct irreducible Laurent factors, this path uses the optimized
    prime-chart seed path: after clearing Laurent monomial units, the
    singleton factor
    ideals are inserted directly as prime torus-saturated components.  For
    later intersections, torus saturation is retained, but an explicit radical
    is omitted before minimal-prime decomposition.  These changes do not alter
    the resulting irreducible component collection.

    Arguments:
        factor_list : A list of Laurent polynomials in the same Laurent
            polynomial ring over QQ.

    Returns:
        A triple (smooth, subset, reason).  On success this is
        (True, None, None).  On failure, subset is the list of zero-based
        factor indexes generating the offending component and reason
        describes the failure.
    """
    if len(factor_list) == 0:
        return (True, None, None)

    d = factor_list[0].parent().ngens()
    name_list = ['Y_' + str(i) for i in range(d)]
    poly_ring = PolynomialRing(QQ, len(name_list), names=name_list)
    identity = identity_matrix(ZZ, d)
    chart_factors = [_polynomial_in_normal_chart(
        g, identity, 0, poly_ring) for g in factor_list]

    smooth, reason, witness = \
        _component_collection_smooth_along_stratum(
            chart_factors, 0, prime_chart_seeds=True)
    if smooth:
        return (True, None, None)

    return (False, witness["factor_subset"], reason)


def _validate_nondegeneracy_input(laurent_poly, verbose):
    """
    Validates common input for non-degeneracy checks.

    Arguments:
        laurent_poly : A nonzero Laurent polynomial over QQ whose Newton
            polytope is full-dimensional.
        verbose : A boolean.

    Returns:
        A pair (newton, ambient_dim), where newton is the Newton polytope of
        laurent_poly and ambient_dim is its dimension.
    """
    poly_parent = laurent_poly.parent()
    if not isinstance(poly_parent, LaurentPolynomialRing_generic):
        raise TypeError("The first argument should be a Laurent polynomial.")
    if not isinstance(verbose, bool):
        raise TypeError("The second argument should be a Boolean.")
    if laurent_poly.is_zero():
        raise ValueError("The Laurent polynomial should be nonzero.")
    if poly_parent.base_ring() is not QQ:
        raise ValueError("The Laurent polynomial should be defined over QQ.")

    newton = newton_polytope(laurent_poly)
    ambient_dim = newton.dim()
    if ambient_dim != poly_parent.ngens():
        raise ValueError("The Newton polytope should be full-dimensional.")

    return (newton, ambient_dim)


def _weak_nondegeneracy_data(laurent_poly, newton, ambient_dim,
                              verbose=False):
    """
    Computes weak non-degeneracy data and caches reduced face polynomials.

    This is the common internal routine used by weak_nondegeneracy() and
    nested_nondegeneracy().  Weak non-degeneracy is checked in the relaxed
    componentwise form used here for every positive-dimensional proper face
    polynomial.  The zero cone of the normal fan, corresponding to the
    full-dimensional torus and the original Laurent polynomial itself, is
    deliberately not tested.  For every proper face, the irreducible
    components of the reduced factor hypersurfaces and of all their
    successive reduced intersections are required to be smooth in the
    corresponding algebraic torus.

    Arguments:
        laurent_poly : A Laurent polynomial over QQ.
        newton : The full-dimensional Newton polytope of laurent_poly.
        ambient_dim : The dimension of newton.
        verbose : A boolean.  If True, every detected failure is printed.

    Returns:
        A triple (output, faces_by_dimension, face_data_cache).  The first
        entry is a dictionary with keys result and failures.  The second
        stores the proper faces of newton by dimension.  The third stores the
        reduced coordinate data for every positive-dimensional proper face.
    """
    failures = []
    weak_ok = True
    faces_by_dimension = {
        face_dim: newton.faces(face_dim)
        for face_dim in range(1, ambient_dim)
    }

    # The zero cone of the normal fan corresponds to the full-dimensional
    # torus.  By definition it is deliberately excluded from the smoothness
    # check, so only proper faces are considered below.  A reduced divisor in a
    # one-dimensional torus is a finite reduced scheme, hence is smooth in
    # characteristic zero; thus one-dimensional faces need no algebraic test.
    face_data_cache = {}
    for face_dim in range(1, ambient_dim):
        faces = faces_by_dimension[face_dim]
        for face_number in range(len(faces)):
            face = faces[face_number]
            data = _face_reduced_data(laurent_poly, face)
            factors = _face_factor_list(data["polynomial"])
            data["factors"] = factors
            face_data_cache[(face_dim, face_number)] = data

            if face_dim >= 2:
                smooth, subset, reason = _factor_intersections_smooth_in_torus(
                    factors)
                if not smooth:
                    weak_ok = False
                    failure = {
                        "type": "weak non-degeneracy",
                        "face_dimension": face_dim,
                        "face_number": face_number,
                        "factor_subset": subset,
                        "reason": reason,
                    }
                    failures.append(failure)
                    if verbose:
                        print(failure)

    output = {
        "result": weak_ok,
        "failures": failures,
    }
    return (output, faces_by_dimension, face_data_cache)


def weak_nondegeneracy(laurent_poly, verbose=False):
    """
    Checks weak non-degeneracy of a Laurent polynomial.

    For every positive-dimensional proper face delta of the Newton polytope,
    write

        F_delta = g_1^e_1 * ... * g_m^e_m

    with distinct irreducible g_i and discard the multiplicities.  Starting
    with the irreducible components of the reduced factor hypersurfaces,
    close the collection under reduced intersection and irreducible
    decomposition.  Every irreducible component obtained at any stage is
    required to be smooth in the corresponding algebraic torus.  Reducible
    intersections are allowed; no transversality condition is imposed.  The
    zero cone, corresponding to the full-dimensional torus and the original
    Laurent polynomial itself, is deliberately excluded.

    Arguments:
        laurent_poly : A nonzero Laurent polynomial over QQ whose Newton
            polytope is full-dimensional.
        verbose : A boolean.  If True, every detected failure is printed.

    Returns:
        A dictionary with keys:
        - result : True if the Laurent polynomial is weakly non-degenerate;
        - failures : a list of diagnostic dictionaries.

        Each diagnostic dictionary has type "weak non-degeneracy" and records
        the face dimension, face number, factor subset and geometric reason.

    Example:
        sage: R.<x,y> = LaurentPolynomialRing(QQ)
        sage: F = x + y + x^-1*y^-1
        sage: weak_nondegeneracy(F)["result"]
        True
    """
    newton, ambient_dim = _validate_nondegeneracy_input(
        laurent_poly, verbose)
    output, _, _ = _weak_nondegeneracy_data(
        laurent_poly, newton, ambient_dim, verbose)
    return output


def is_weakly_nondegenerate(laurent_poly, verbose=False):
    """
    Returns whether a Laurent polynomial is weakly non-degenerate.

    This is a Boolean wrapper around weak_nondegeneracy().  Use the latter
    function when diagnostic information is required.

    Arguments:
        laurent_poly : A nonzero Laurent polynomial over QQ whose Newton
            polytope is full-dimensional.
        verbose : A boolean.  If True, every detected failure is printed.

    Returns:
        A boolean.

    Example:
        sage: R.<x,y> = LaurentPolynomialRing(QQ)
        sage: F = x + y + x^-1*y^-1
        sage: is_weakly_nondegenerate(F)
        True
    """
    return weak_nondegeneracy(laurent_poly, verbose)["result"]


def nested_nondegeneracy(laurent_poly, verbose=False):
    """
    Checks whether a Laurent polynomial is nested non-degenerate.

    Let P be the Newton polytope of laurent_poly.  First, relaxed weak
    non-degeneracy on the positive-dimensional proper faces of P is checked
    by calling the same internal routine used by weak_nondegeneracy().  The
    zero cone, corresponding to the full-dimensional torus, is not tested.

    The additional nested-face condition is imposed for every proper face
    delta and every positive-dimensional proper subface tau < delta.  Vertex
    substrata are deliberately excluded.  For each such pair the relative
    normal cone N_{tau/delta} is required to be unimodular in the saturated
    relative lattice.  No subdivision is chosen when this condition fails.

    If c = codim_delta(tau), the primitive normal rays are completed to a
    unimodular coordinate system z_1, ..., z_c, y_1, ..., y_(dim(delta)-c).
    The factors of the reduced face polynomial F_delta are extended to the
    resulting chart

        A^c x (G_m)^(dim(delta)-c),

    and their closures are placed in one common collection.  This collection
    is closed recursively under reduced intersection and irreducible
    decomposition.  Every irreducible component obtained at any stage is
    required to be smooth along z_1 = ... = z_c = 0.  Thus reducible reduced
    intersections are allowed when their components and all components of
    their successive reduced intersections are smooth along the stratum.

    Multiplicities of irreducible factors are ignored throughout.

    Arguments:
        laurent_poly : A nonzero Laurent polynomial over QQ whose Newton
            polytope is full-dimensional.
        verbose : A boolean.  If True, every detected failure is printed.

    Returns:
        A dictionary with keys:
        - result : True if all conditions hold;
        - weak_nondegenerate : whether weak non-degeneracy passes;
        - relative_unimodular : whether every relevant relative normal cone
          is unimodular;
        - nested_substrata : whether the global collection generated by the
          face-factor closures and closed under successive reduced
          intersections has only components smooth along every nested
          substratum with a unimodular normal chart;
        - failures : a list of diagnostic dictionaries.

        A diagnostic dictionary has type "weak non-degeneracy",
        "relative unimodularity", or "nested substratum".  Face and subface
        numbers are the zero-based indexes returned by Sage for the relevant
        face lists.  For nested failures the report also records the factor
        subset, normal rays, chart equations and the geometric reason.

    Example:
        sage: R.<x,y> = LaurentPolynomialRing(QQ)
        sage: F = x + y + x^-1*y^-1
        sage: output = nested_nondegeneracy(F)
        sage: output["result"]
        True
        sage: output["weak_nondegenerate"]
        True
    """
    newton, ambient_dim = _validate_nondegeneracy_input(
        laurent_poly, verbose)
    weak_output, faces_by_dimension, face_data_cache = \
        _weak_nondegeneracy_data(
            laurent_poly, newton, ambient_dim, verbose)

    failures = list(weak_output["failures"])
    weak_ok = weak_output["result"]
    unimodular_ok = True
    nested_ok = True

    # The nested-face condition is nonempty only for dim(delta) >= 2.
    for face_dim in range(2, ambient_dim):
        faces = faces_by_dimension[face_dim]
        for face_number in range(len(faces)):
            data = face_data_cache[(face_dim, face_number)]
            face_polyhedron = data["polyhedron"]
            factors = data["factors"]

            for subface_dim in range(1, face_dim):
                subfaces = face_polyhedron.faces(subface_dim)
                for subface_number in range(len(subfaces)):
                    subface = subfaces[subface_number]
                    cone_data = _relative_normal_cone_data(
                        face_polyhedron, subface)

                    if not cone_data["unimodular"]:
                        unimodular_ok = False
                        failure = {
                            "type": "relative unimodularity",
                            "face_dimension": face_dim,
                            "face_number": face_number,
                            "face_vertices": list(face_polyhedron.vertices()),
                            "subface_dimension": subface_dim,
                            "subface_number": subface_number,
                            "subface_vertices": list(subface.vertices()),
                            "rays": cone_data["rays"],
                            "simplicial": cone_data["simplicial"],
                            "lattice_index": cone_data["lattice_index"],
                        }
                        failures.append(failure)
                        if verbose:
                            print(failure)
                        # Unimodular normal coordinates do not exist in
                        # this chart.
                        continue

                    normal_dim = cone_data["dimension"]
                    normal_matrix = _unimodular_normal_matrix(
                        cone_data["rays"], face_dim)
                    name_list = (['z_' + str(i) for i in range(normal_dim)] +
                                 ['y_' + str(i)
                                  for i in range(face_dim - normal_dim)])
                    chart_ring = PolynomialRing(
                        QQ, len(name_list), names=name_list)
                    chart_factors = [_polynomial_in_normal_chart(
                        g, normal_matrix, normal_dim, chart_ring)
                        for g in factors]

                    if len(chart_factors) > 0:
                        smooth, reason, witness = \
                            _component_collection_smooth_along_stratum(
                                chart_factors, normal_dim,
                                prime_chart_seeds=True)
                        if not smooth:
                            nested_ok = False
                            failure = {
                                "type": "nested substratum",
                                "face_dimension": face_dim,
                                "face_number": face_number,
                                "face_vertices": list(
                                    face_polyhedron.vertices()),
                                "subface_dimension": subface_dim,
                                "subface_number": subface_number,
                                "subface_vertices": list(
                                    subface.vertices()),
                                "factor_subset": witness["factor_subset"],
                                "normal_rays": cone_data["rays"],
                                "chart_equations": chart_factors,
                                "component_generators": witness[
                                    "component_generators"],
                                "reason": reason,
                            }
                            failures.append(failure)
                            if verbose:
                                print(failure)

    return {
        "result": weak_ok and unimodular_ok and nested_ok,
        "weak_nondegenerate": weak_ok,
        "relative_unimodular": unimodular_ok,
        "nested_substrata": nested_ok,
        "failures": failures,
    }


def is_nested_nondegenerate(laurent_poly, verbose=False):
    """
    Returns whether a Laurent polynomial is nested non-degenerate.

    This is a Boolean wrapper around nested_nondegeneracy().  Use the
    latter function when diagnostic information is required.

    Arguments:
        laurent_poly : A nonzero Laurent polynomial over QQ whose Newton
            polytope is full-dimensional.
        verbose : A boolean.  If True, every detected failure is printed.

    Returns:
        A boolean.

    Example:
        sage: R.<x,y> = LaurentPolynomialRing(QQ)
        sage: F = x + y + x^-1*y^-1
        sage: is_nested_nondegenerate(F)
        True
    """
    return nested_nondegeneracy(laurent_poly, verbose)["result"]
