from random import seed
from random import randint
from random import random
import numpy
from numpy import copy
from numpy.linalg import inv
seed(23)
matrix=[]
n=7
def make_square_matrix(n): # n rows m colunms matrix
    A=[]
    for i in range(n):
        row=[]
        for j in range(n):
            f=randint(0,n)
            row.append(f)
        A.append(row)
    return A
def test_matrix(A):
    for i,l1 in enumerate(A):
        for j,l2 in enumerate(l1):
            if A[i][j] == A[j][i]:
                print ("Matrix element ",A[i][j],' is equal to',A[j][i],'for row = ',i,' and colunm = ',j)
            else:
                print("Matrix is Non Symmetric!")
def print_matrix(A,name):
    print('          '+ name)
    for i in A:
        print(i)
def custom_matrix(n):
    A=[]
    for i in range(n):
        k=i+2
        print('Starting loop in i = ', i)
        row=[]
        print(k,' Value of k')
        for j in range(n):
            print(' Starting loop in j = ',j)
            elem=(k)**(k-1)
            row.append(elem)
            k+=1
        A.append(row)
    return A
def matrix_mult(A,B):
    # Number of rows and colunms and an empty list 'C'
    C=[]
    rA=len(A)
    cA=len(A[0])
    rB=len(B)
    cB=len(B[0])
    # Nested Loop
    for i in range(rA): # i ranges in A rows
        for j in range(cB): # j range B colunms
            for k in range(cA): # k range in A colunms
                C[i][j]+=A[i][k] * B[k][j]
    return C



A=make_square_matrix(n)
# List comprehension to make transpose matrix
A_t = [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]
# A=custom_matrix(n)
# print(A)
print_matrix(A,'A')
print_matrix(A_t,'A transpose')
RESP_A=matrix_mult(A,A_t)
print_matrix(RESP_A,'Resposta A')

