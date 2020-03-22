import sys
from random import seed
from random import randint
from random import random
import numpy as np
from numpy import copy
from numpy.linalg import inv
from numpy.linalg import solve
from numpy.linalg import det
import copy
seed(7)
matrix=[]
n=7

def make_square_matrix(n): # n rows m colunms matrix
    """Make an NxN square matrix with non-null determinant
    """
    A=[]
    for i in range(n):
        row=[]
        for j in range(n):
            f=randint(0,n)
            row.append(f)
        A.append(row)
    return A
def test_matrix(A,name): # Compare A with A transpose
    """ Count how many rows are equal in A and A transpose
    If all rows are equal, an exception will be raised
    """
    rows_equal=0
    A_t=[[A[j][i]
        for j in range(len(A))] 
        for i in range(len(A[0]))]
    for row in range(n):
        if A_t[row] == A[row]:
            rows_equal+=1
        else:
            pass
    if rows_equal == n:
        raise Exception('Symmetric Matrix!')
def print_matrix(A,name):
    """ Print each row of a matrix 
    A   : matrix name, list
    name: string of the matrix name to print
    """
    print('          '+ name)
    for i in range(n):
        # Format float numbers when printing matrices
        i_fmt=["%.2f" % col_element for col_element in A[i]]
        print(i_fmt)
    print()
def matrix_mult(A,B):
    """ Multiply two matrices A and B into C
    rA  : A rows
    cA  : A colunms
    rB  : B rows
    cB  : B colunms
    C   :  Initially an NxN zero matrix
    """
    # Number of rows and colunms and a list of 0 in 'C'
    rA=len(A)
    cA=len(A[0])
    rB=len(B)
    cB=len(B[0])
    # List comprehension of a zero matrix NxN
    C=[[0 for row in range(cB)] for col in range(rA)]
    # Nested Loop to substitue values of the zero matrix
    for i in range(rA): # i ranges in A rows
        for j in range(cB): # j range B colunms
            for k in range(cA): # k range in A colunms
                C[i][j]+=A[i][k]*B[k][j]
    return C

def copy_matrix(A):
    """ A   : matrix to be copied
    Returns :  A copy of matrix A
    """
    row=len(A)
    col=len(A[0])
    A_copy=[[0 for i in range(row)] for j in range(col)]

    for i in range(row):
        for j in range(col):
            A_copy[i][j]=A[i][j]
    return A_copy

def determinant(A,t=0):
    """ Calculate a determinant recursively
    A         : Matrix to calculate the determinant
    t         : The total value of the determinant 
    row       : Row loop element
    col       : Col loop element
    sub_matrix: Cut off the first row
    N_cols    : Length of sub_matrix
    sign      : Depends on the index, (even > 0, odd < 0)
    sub_det   : Recursive determinant of sub_matrix
    """
    index=list(range(len(A)))
    # 2X2 Matrix Case
    if len(A) == 2 and len(A[0]) == 2:
        det_2d=A[0][0]*A[1][1]-A[0][1]*A[1][0]
        return det_2d
    for row in index:
        sub_matrix=A[1:] # Cut the first row
        N_cols=len(sub_matrix)
       # print(N_cols)
       # print(sub_matrix,' Submatrix of N_cols = ',N_cols)
        for col in range(N_cols): # loop on colunms
            sub_matrix[col]=sub_matrix[col][0:row]+sub_matrix[col][row+1:] # Remove colunm
        sign=(-1)**(row%2) # sub_det is None?
        sub_det=determinant(sub_matrix)
        t+=sign*A[0][row]*sub_det 
    return t

def transpose_matrix(A):
    """ Calculate the transpose matrix of a 2D array
    A         : Input matrix, 2D array
    A_t       : Tranpose matrix, 2D array

    """
    A_t=[]
    for i in range(n):
        row=[]
        for j in range(n):
            row.append(A[j][i])
        A_t.append(row)
    return A_t

def get_minor(A,i,j):
    """ Get minor of an input matrix A
    A         : Input Matrix
    i,j       : Rows and colunms
    """
    return [row[:j] + row[j+1:] for row in (A[:i]+A[i+1:])]
    

def get_inverse(A):
    """ Calculate the inverse matrix using 
    A         : Input matrix
    cofactors : Elements used to calculate the inverse matrix
    """
    d=determinant(A)
    cofactors=[]
    for r in range(len(A)):
        cofactor_row=[]
        for c in range(len(A)):
            minor=get_minor(A,r,c)
            cofactor_row.append(((-1)**(r+c))*determinant(minor))
        cofactors.append(cofactor_row)
    cofactors=transpose_matrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c]=cofactors[r][c]/d
    return cofactors 

def gauss_elim(A,b):
    """ Do gaussian elimination with matrix A and vector b
    A       : Input matrix
    b       : input vector
    Returns : b (output vector), det (determinant-> float)
    """
    a=copy_matrix(A) # Make a copy of matrix A and vector b
    b=copy.deepcopy(b)
    n=len(a)
    p=len(b[0])
    det=1 # Gaussian elimination
    for i in range(n-1):
        k=i
        for j in range(i+1,n):
            if abs(a[j][i]) > abs(a[k][i]):
                k=j
        if k != i:
            a[i],a[k]=a[k],a[i]
            b[i],b[k]=b[k],b[i]
            det=-det

        for j in range(i+1,n):
            t = a[j][i]/a[i][i]
            for k in range(i+1,n):
                a[j][k]-=t*a[i][k]
            for k in range(p):
                b[j][k]-=t*b[i][k]

    for i in range(n-1,-1,-1):
        for j in range(i+1,n):
            t=a[i][j]
            for k in range(p):
                b[i][k]-=t*b[j][k]
        t=1/a[i][i]
        det*=a[i][i]
        for j in range(p):
            b[i][j]*=t
    return det, b
if __name__ == "__main__":
    A=make_square_matrix(n)
    test_matrix(A,'A') # Matrix A is being tested as singular
    print("{0:.2f}".format(round(determinant(A))),' Calculated determinant')
    print("{0:.2f}".format(round(det(A))),' Numpy determinant')
    A_t=transpose_matrix(A)
# Answers below
    print_matrix(A,'A')
    print_matrix(A_t,'A transpose')
    b=[[1],[2],[3],[4],[5],[6],[7]]
    print('b vector = ',b)
    RESP_A=matrix_mult(A,A_t)
    print_matrix(RESP_A,'Answer A: A times A transpose = ')
    RESP_B=get_inverse(A)
    AiA=(matrix_mult(RESP_B,A))
    print_matrix(AiA,' Inverse of A times A')
    print_matrix(np.dot(A,A_t),' A times A transpose with Numpy = ')
    print_matrix(RESP_B,'Answer B: Inverse of A is = ') 
    det_C,x=gauss_elim(A,b)
    print_matrix(x,' x, Answer C,  with Gaussian Elimination ')
    print_matrix(solve(A,b),' x, Answer C, with Numpy       ')
    AAi=(matrix_mult(A,RESP_B))
    print_matrix(AAi,' A times inverse of A')



