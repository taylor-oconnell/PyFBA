import sys
import glpk

"""

Run linear programming using GLPK. This uses the updated pyGLPK library.

This is a generic linear programming wrapper that I wrote that we can
build upon for fba work, but it is not limited to fba.

Do not use the standard pyGLK. The only version that I could get to
compile is https://github.com/bradfordboyle/pyglpk and that will only
compile on python2.x

That means this code is limited to python2.7.


"""

solver = glpk.LPX()


def load(matrix, rowheaders=None, colheaders=None, verbose=0, likelihood_gapfill=False):
    """
    Load the data matrix into the linear programming solver

    :param matrix: the 2D array of data. It should not have row or column
    headers, they can be specified separately
    :type matrix: list of list
    :param rowheaders: (optional) are the row identifiers
    :type rowheaders: list
    :param colheaders: (optional) are the column identifiers
    :type colheaders: list
    :param verbose: verbose turns on some debugging output. The higher the number the more output is generated
    :type verbose: int
    :param likelihood_gapfill: Run in likelihood-based gapfill mode
    :type likelihood_gapfill: bool
    :return: void
    :rtype: void

    """
    global solver

    solver.erase()
    
    if likelihood_gapfill:
        solver.obj.maximize = False
    else:
        solver.obj.maximize = True

    nrows = len(matrix)
    ncols = len(matrix[0])

    if verbose > 0:
        sys.stderr.write("We are loading " + str(nrows) + " rows and " + str(ncols) + " columns\n")

    solver.rows.add(nrows)
    solver.cols.add(ncols)

    # we need to flatten the 2D array before we add it to the lp object
    temp = []
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            temp.append(matrix[i][j])

    if verbose > 4:
        sys.stderr.write("Matrix: " + str(temp) + "\n")
    solver.matrix = temp

    # name the rows and columns
    if rowheaders and len(rowheaders) == nrows:
        for i in range(len(rowheaders)):
            if len(rowheaders[i]) > 255:
                if verbose > 0:
                    sys.stderr.write("WARNING ROW HEADER: " + str(rowheaders[i]) + " truncated to 255 characters\n")
                solver.rows[i].name = rowheaders[i][0:255]
            else:
                solver.rows[i].name = rowheaders[i]
    elif rowheaders:
        raise ValueError("The size of row headers (" + str(len(rowheaders)) +
                         ") does not match the expected number of rows (" + str(nrows) + "\n")

    if colheaders and len(colheaders) == ncols:
        for i in range(len(colheaders)):
            if len(colheaders[i]) > 255:
                if verbose > 0:
                    sys.stderr.write("WARNING ROW HEADER: " + str(colheaders[i]) + " truncated to 255 characters\n")
                solver.cols[i].name = colheaders[i][0:255]
            else:
                solver.cols[i].name = colheaders[i]
    elif colheaders:
        raise ValueError("Warning: the size of col headers (" + str(len(colheaders)) +
                         ") does not match the expected number of cols (" + str(ncols) + "\n")


def row_bounds(bounds):
    """
    Set the bounds for the rows in the linear programming. 
    This should be an array of the same length as the number of rows, 
    and each element should be a tuple of (lower bound, upper bound)

    :param bounds: The bounds as a single tuple for each of the rows
    :type bounds: list of tuples
    :return: void
    :rtype: void

    """

    global solver
    if len(bounds) != len(solver.rows):
        raise ValueError("There must be the same number of bounds as rows bounds:" + str(bounds) + " rows: " + str(len(
            solver.rows)) + "\n")

    for i in range(len(bounds)):
        solver.rows[i].bounds = bounds[i]


def col_bounds(bounds):
    """
    Set the bounds for the columns in the linear programming.
    This should be an array of the same length as the number of columns,
    and each element should be a tuple of (lower bound, upper bound)

    :param bounds: The bounds as a single tuple for each of the columns
    :type bounds: list of tuples
    :return: void
    :rtype: void
    """

    global solver
    if len(bounds) != len(solver.cols):
        raise ValueError("There must be the same number of bounds as cols")

    for i in range(len(bounds)):
        solver.cols[i].bounds = bounds[i]


def objective_coefficients(coeff):
    """
    Set the objective coefficients. coeff should be an array of
    coefficients

    :param coeff: The objective cooefficient for the linear solver
    :type coeff: list of float
    :return: void
    :rtype: void
    """
    global solver
    solver.obj[:] = coeff


def solve():
    """
    Solve the lp and return the status and the objective function
    value

    :return: The status and value of the solution
    :rtype: str, float

    """
    solver.simplex()
    return solver.status, solver.obj.value


def col_primal_hash():
    """
    Return a hash of the column names and the primals (activities)
    associated with those columns. This presumes that you have named
    the columns

    :return: A hash of the column names and their primals
    :rtype: dict
    """

    d = {}
    for c in solver.cols:
        d[c.name] = c.primal
    return d


def col_primals():
    """
    Return an array of the primals (activities), one for each column

    :return: A list of the column primals
    :rtype: list
    """

    d = []
    for c in solver.cols:
        d.append(c.primal)
    return d


def row_primal_hash():
    """ Retrieve a hash of the primals (activity) of the rows. This
    presume that you have named the columns

    :return: A hash of the row names and their primals
    :rtype: dict
    """

    d = {}
    for r in solver.rows:
        d[r.name] = r.primal
    return d


def row_primals():
    """
    Return an array of the primals (activities), one for each column

    :return: A list of the row primals
    :rtype: list
    """

    d = []
    for r in solver.rows:
        d.append(r.primal)
    return d


