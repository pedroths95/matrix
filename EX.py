import sys
from random import seed
from random import randint
from random import random
import numpy
from numpy import copy
from numpy.linalg import inv
seed(7)
matrix=[]
n=7
def make_square_matrix(n): # n rows m colunms matrix
    """Make an NxN square matrix with non-null determinant
    Loop in i: For each row
    Loop in j: For each colunm
    If i == j: Elements of the diagonal
    else     : Off-diagonal elements
    """
    A=[]
    row=[]
    for i in range(n):
        for j in range(n):
            if i == j: # Non null determinant
                f=random()
                #row.append("{0:.2f}".format(f))
            else:
                seed(9)
                g=random()
                #row.append('{0:+7.3f}'.format(g))
        A.append(row)
    print_matrix(A,'A')
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
        if A_t[row] == A_t[row]:
            rows_equal+=1
        else:
            pass
    if rows_equal == n:
        raise Exception('Symmetric Matrix!')
def print_matrix(A,name):
    """ Print each row of a matrix 
    A: matrix name, list
    name: string of the matrix name to print
    """
    print('          '+ name)
    for i in A:
        print(i)
def matrix_mult(A,B):
    """ Multiply two matrices A and B into C
    rA: A rows
    cA: A colunms
    rB: B rows
    cB: B colunms
    C:  Initially an NxN zero matrix
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

def determinant(A,t=0):
    """ Calculate a determinant recursively
    A         : Matrix to calculate the determinant
    t         : The total value of the determinant
    row       : Index for the rows
    sub_matrix: Cut off the first row
    N_cols    : Length of sub_matrix
    col       : Index for the colunms
    sign      : Depends on the index, (even > 0, odd < 0)
    sub_det   : Recursive determinant of sub_matrix
    """

    index=list(range(len(A)))
    # 2X2 Matrix
    if len(A) == 2 and len(A[0]) == 2:
        det_2d=A[0][0]*A[1][1]-A[0][1]*A[1][0]
        return det_2d
    for row in index:
        sub_matrix=A[1:] # Cut the first row
        N_cols=len(sub_matrix)
        for col in range(N_cols): # loop on colunms
            sub_matrix[col]=sub_matrix[col][0:row]+sub_matrix[col][row+1:] # Remove colunm
        sign=(-1)**(row%2) # sub_det is None?
        sub_det=determinant(sub_matrix)
        t+=sign*A[0][row]*sub_det 
    return t




# Make the matrix and print it
A=make_square_matrix(n)
#test_matrix(A,'A') # Matrix A is being tested
# List comprehension to make transpose matrix
A_t=[[A[j][i] 
    for j in range(len(A))]
    for i in range(len(A[0]))]
print(determinant(A))


# Answers below
#print_matrix(A,'A')
#print_matrix(A_t,'A transpose')
#RESP_A=matrix_mult(A,A_t)
#print_matrix(RESP_A,'Answer A')

