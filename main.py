import math
import copy
import time
import matplotlib.pyplot as plt

def print_matrix(MatrixA):
    if type(MatrixA[0]) == int or type(MatrixA[0]) == float:
        for i in range(0, len(MatrixA)):
            if type(MatrixA[i]) == float:
                print("%.7f" % MatrixA[i], end="\t")
            else:
                print(MatrixA[i], end="\t")
            print()
        print()
    else:
        for i in range(0, len(MatrixA)):
            for j in range(0, len(MatrixA[0])):
                if type(MatrixA[i][j]) == float:
                    print("%.7f" % MatrixA[i][j], end="\t")
                else:
                    print(MatrixA[i][j], end="\t")
            print()
        print()


def euclidean_norm(MatrixA):
    if type(MatrixA[0]) == int or type(MatrixA[0]) == float:
        suma = 0
        for i in range(0, len(MatrixA)):
            suma = suma + abs(MatrixA[i] * MatrixA[i])
        wynik = math.sqrt(suma)
        return wynik
    else:
        suma = 0
        for i in range(0, len(MatrixA)):
            for j in range(0, len(MatrixA[0])):
                suma = suma + abs(MatrixA[i][j] * MatrixA[i][j])
        wynik = math.sqrt(suma)
        return wynik


def negative_matrix(MatrixA):
    if type(MatrixA[0]) == int or type(MatrixA[0]) == float:
        wyniki = [0 for i in range(len(MatrixA))]
        for i in range(0, len(MatrixA)):
            wyniki[i] = (-1) * MatrixA[i]
        return wyniki
    else:
        wyniki = [[0 for j in range(len(MatrixA[0]))] for i in range(len(MatrixA))]
        for i in range(0, len(MatrixA)):
            for j in range(0, len(MatrixA[0])):
                wyniki[i][j] = (-1) * MatrixA[i][j]
        return wyniki


def invert_diagonal_matrix(MatrixA):
    if type(MatrixA[0]) == int or type(MatrixA[0]) == float:
        wyniki = [0 for i in range(len(MatrixA))]
        for i in range(0, len(MatrixA)):
            wyniki[i] = 1 / (MatrixA[i])
        return wyniki
    else:
        wyniki = [[0 for j in range(len(MatrixA[0]))] for i in range(len(MatrixA))]
        for i in range(0, len(MatrixA)):
            for j in range(0, len(MatrixA[0])):
                if i == j:
                    wyniki[i][j] = 1 / (MatrixA[i][j])
                else:
                    wyniki[i][j] = 0
        return wyniki


def matrix_addition(MatrixA, MatrixB):
    if type(MatrixA[0]) == int or type(MatrixA[0]) == float or type(MatrixB[0]) == int or type(MatrixB[0]) == float:
        if len(MatrixA) == len(MatrixB):
            wyniki = [0 for i in range(len(MatrixA))]
            for i in range(len(MatrixA)):
                wyniki[i] = MatrixA[i] + MatrixB[i]
            return wyniki
    else:
        wyniki = [[0 for j in range(len(MatrixA[0]))] for i in range(len(MatrixA))]
        for i in range(0, len(MatrixA)):
            for j in range(0, len(MatrixA[0])):
                wyniki[i][j] = MatrixA[i][j] + MatrixB[i][j]
        return wyniki


def matrix_subtraction(MatrixA, MatrixB):
    if (type(MatrixA[0]) == int or type(MatrixA[0]) == float) and (
            type(MatrixB[0]) == int or type(MatrixB[0]) == float):
        if len(MatrixA) == len(MatrixB):
            wyniki = [0 for i in range(len(MatrixA))]
            for i in range(len(MatrixA)):
                wyniki[i] = MatrixA[i] - MatrixB[i]
            return wyniki
    else:
        wyniki = [[0 for j in range(len(MatrixA[0]))] for i in range(len(MatrixA))]
        for i in range(0, len(MatrixA)):
            for j in range(0, len(MatrixA[0])):
                wyniki[i][j] = MatrixA[i][j] - MatrixB[i][j]
        return wyniki


def matrix_multiplication(MatrixA, MatrixB):
    if (type(MatrixA[0]) == int or type(MatrixA[0]) == float) and (
            type(MatrixB[0]) == int or type(MatrixB[0]) == float):
        if len(MatrixA[0]) == len(MatrixB):
            wyniki = [0 for i in range(len(MatrixA))]
            for i in range(len(MatrixA)):
                for j in range(len(MatrixB)):
                    wyniki[i] += MatrixA[i][j] * MatrixB[j]
            return wyniki
    elif type(MatrixA[0]) == int or type(MatrixA[0]) == float:
        wyniki = [0 for i in range(len(MatrixA))]
        for i in range(len(MatrixA)):
            for j in range(len(MatrixB)):
                wyniki[i] += MatrixA[j] * MatrixB[i][j]
        return wyniki
    elif type(MatrixB[0]) == int or type(MatrixB[0]) == float:
        wyniki = [0 for i in range(len(MatrixA))]
        for i in range(len(MatrixA)):
            for j in range(len(MatrixB)):
                wyniki[i] += MatrixA[i][j] * MatrixB[j]
        return wyniki
    else:
        wyniki = [[0 for j in range(len(MatrixB[0]))] for i in range(len(MatrixA))]
        for i in range(len(MatrixA)):
            for j in range(len(MatrixB[0])):
                for k in range(len(MatrixB)):
                    wyniki[i][j] += MatrixA[i][k] * MatrixB[k][j]
        return wyniki


def norm(MatrixA, MatrixB, MatrixR):
    FIRSTM = matrix_multiplication(MatrixA, MatrixR)
    SUB = matrix_subtraction(FIRSTM, MatrixB)
    return euclidean_norm(SUB)


def jacobi_method(MatrixA, MatrixB, MatrixR, MatrixD, MatrixL, MatrixU):
    ID = invert_diagonal_matrix(MatrixD)
    NID = negative_matrix(ID)
    SUM = matrix_addition(MatrixL, MatrixU)
    FIRSTM = matrix_multiplication(NID, SUM)
    THIRDM = matrix_multiplication(ID, MatrixB)
    start = time.time()
    while norm(MatrixA, MatrixB, MatrixR) > 10 ** (-9):
        SECONDM = matrix_multiplication(FIRSTM, MatrixR)
        MatrixR = matrix_addition(SECONDM, THIRDM)
    end = time.time()
    czas = end - start
    return MatrixR, czas


def jacobi_method_iterations(MatrixA, MatrixB, MatrixR, MatrixD, MatrixL, MatrixU):
    iterations = 0
    # normArray = []
    ID = invert_diagonal_matrix(MatrixD)
    NID = negative_matrix(ID)
    SUM = matrix_addition(MatrixL, MatrixU)
    FIRSTM = matrix_multiplication(NID, SUM)
    THIRDM = matrix_multiplication(ID, MatrixB)
    start = time.time()
    while norm(MatrixA, MatrixB, MatrixR) > 10 ** (-9) and iterations < 100:
        SECONDM = matrix_multiplication(FIRSTM, MatrixR)
        MatrixR = matrix_addition(SECONDM, THIRDM)
        # normArray.append(norm(MatrixA, MatrixB, MatrixR))
        iterations = iterations+1
    end = time.time()
    czas = end - start
    return iterations, czas, norm(MatrixA, MatrixB, MatrixR)


def forward_substitution(MatrixA, MatrixB):
    wyniki = [0 for i in range(len(MatrixB))]
    wyniki[0] = MatrixB[0]/MatrixA[0][0]
    for i in range(1, len(MatrixB)):
        temp = MatrixB[i]
        for j in range(0, i):
            temp = temp - MatrixA[i][j] * wyniki[j]
        wyniki[i] = temp/MatrixA[i][i]
    return wyniki


def backward_substitution(MatrixA, MatrixB):
    n = len(MatrixB) - 1
    wyniki = [0 for i in range(len(MatrixB))]
    wyniki[n] = MatrixB[n]/MatrixA[n][n]
    for i in range(1, n+1):
        j = n-i
        temp = 0
        for k in range(0, i):
            temp = temp + MatrixA[j][n-k] * wyniki[n-k]
        wyniki[j] = (MatrixB[j] - temp)/MatrixA[j][j]
    return wyniki


def gauss_seidel_method(MatrixA, MatrixB, MatrixR, MatrixD, MatrixL, MatrixU):
    DL = matrix_addition(MatrixL, MatrixD)
    NDL = negative_matrix(DL)
    SECONDFS = forward_substitution(DL, MatrixB)
    start = time.time()
    while norm(MatrixA, MatrixB, MatrixR) > 10 ** (-9):
        UR = matrix_multiplication(MatrixU, MatrixR)
        FIRSTFS = forward_substitution(NDL, UR)
        MatrixR = matrix_addition(FIRSTFS, SECONDFS)
    end = time.time()
    czas = end - start
    return MatrixR, czas


def gauss_seidel_method_iterations(MatrixA, MatrixB, MatrixR, MatrixD, MatrixL, MatrixU):
    iterations = 0
    # normArray = []
    DL = matrix_addition(MatrixL, MatrixD)
    NDL = negative_matrix(DL)
    SECONDFS = forward_substitution(DL, MatrixB)
    start = time.time()
    while norm(MatrixA, MatrixB, MatrixR) > 10 ** (-9) and iterations < 100:
        UR = matrix_multiplication(MatrixU, MatrixR)
        FIRSTFS = forward_substitution(NDL, UR)
        MatrixR = matrix_addition(FIRSTFS, SECONDFS)
        # normArray.append(norm(MatrixA, MatrixB, MatrixR))
        iterations = iterations + 1
    end = time.time()
    czas = end - start
    return iterations, czas, norm(MatrixA, MatrixB, MatrixR)


def lu_decomposition(MatrixTest):
    lower = create_diagonal(len(MatrixTest), 1)
    '''
    upper = [[0 for j in range(MatrixTest[0])] for i in range(MatrixTest)]
    for i in range(0, len(MatrixTest)):
        for j in range(0, len(MatrixTest[0]))
    '''
    upper = copy.deepcopy(MatrixTest)
    for k in range(0, len(MatrixTest) - 1):
        for j in range(k+1, len(MatrixTest)):
            lower[j][k] = upper[j][k] / upper[k][k]
            for d in range(k, len(MatrixTest)):
                upper[j][d] = upper[j][d] - lower[j][k]*upper[k][d]
    return lower, upper


def elimination_lower(MatrixL, MatrixB):
    wyniki = [0 for i in range(len(MatrixB))]
    wyniki[0] = MatrixB[0]
    for i in range(0, len(MatrixB)):
        suma = 0
        for j in range(0, i):
            suma = suma + MatrixL[i][j] * wyniki[j]
        wyniki[i] = MatrixB[i] - suma
    return wyniki


def elimination_upper(MatrixU, MatrixZ):
    wyniki = [0 for i in range(len(MatrixZ))]
    n = len(MatrixZ)-1
    wyniki[n] = MatrixZ[n]/MatrixU[n][n]
    for i in reversed(range(0, len(MatrixZ))):
        suma = 0
        for j in reversed(range(i+1, len(MatrixZ))):
            suma = suma + MatrixU[i][j] * wyniki[j]
        wyniki[i] = (MatrixZ[i] - suma)/MatrixU[i][i]
    return wyniki


def lu_method(MatrixA, MatrixB):
    LOWER, UPPER = lu_decomposition(MatrixA)
    Y = forward_substitution(LOWER, MatrixB)
    X = backward_substitution(UPPER, Y)
    return X


def lu_method_residual(MatrixA, MatrixB):
    start = time.time()
    LOWER, UPPER = lu_decomposition(MatrixA)
    # TEST = matrix_multiplication(LOWER, UPPER)
    Y = forward_substitution(LOWER, MatrixB)
    X = backward_substitution(UPPER, Y)
    end = time.time()
    czas = end - start
    return norm(MatrixA, MatrixB, X), czas


def create_matrixA(N, a1, a2, a3):
    wyniki = [[0 for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                wyniki[i][j] = a1
            elif i == j + 1:
                wyniki[i][j] = a2
            elif i + 1 == j:
                wyniki[i][j] = a2
            elif i == j + 2:
                wyniki[i][j] = a3
            elif i + 2 == j:
                wyniki[i][j] = a3
    return wyniki


def create_matrixB(N):
    wyniki = []
    for i in range(N):
        wyniki.append(math.sin(9 * i))
    return wyniki


def create_diagonal_matrix(N, MatrixA):
    wyniki = [[0 for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                wyniki[i][j] = MatrixA[i][j]
    return wyniki


def create_diagonal(N, value):
    wyniki = [[0 for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                wyniki[i][j] = value
    return wyniki


def create_lower_triangular_matrix(N, MatrixA):
    wyniki = [[0 for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            if i > j:
                wyniki[i][j] = MatrixA[i][j]
    return wyniki


def create_upper_triangular_matrix(N, MatrixA):
    wyniki = [[0 for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            if i < j:
                wyniki[i][j] = MatrixA[i][j]
    return wyniki


iteracjeJ = [23, 24, 24, 24, 25, 25]
iteracjeGS = [16, 17, 17, 17, 17, 17]
iloscN = [100, 500, 901, 1000, 2000, 3000]

plt.plot(iloscN, iteracjeJ)
plt.plot(iloscN, iteracjeGS)
plt.title("Ilosc iteracji w zależności od metody i N")
plt.xlabel("Rozmiar boku macierzy")
plt.ylabel("Ilość iteracji")
plt.legend(["Ilość iteracji dla metody Jacobiego                                             ", "Ilość iteracji dla metody Gaussa-Seidla                                      "])
plt.savefig("IteracjeJacobiGaussSeidl")
plt.clf()


timeJ = [0.04024362564086914, 0.8497908115386963, 2.7541346549987793, 3.389817476272583, 14.13434648513794, 40.0833740234375]
timeGS = [0.0299835205078125, 0.7420492172241211, 2.4115970134735107, 2.9634106159210205, 12.066675424575806, 33.114846395492554]
plt.plot(iloscN, timeJ)
plt.plot(iloscN, timeGS)
plt.title("Czas w zależności od metody i rozmiaru boku macierzy")
plt.xlabel("Rozmiar boku macierzy")
plt.ylabel("Czas w sekundach")
plt.legend(["Metoda Jacobiego", "Metoda Gaussa-Seidla"])
plt.savefig("CzasJacobiSeidl")
plt.clf()


timeLU = [0.029994964599609375, 3.715615749359131, 22.965853691101074, 29.2647602558136, 232.4151086807251, 890.6725356578827]
plt.semilogy(iloscN, timeLU)
plt.title("Czas faktoryzacji LU w zależności od rozmiaru boku macierzy")
plt.xlabel("Rozmiar boku macierzy")
plt.ylabel("Czas w sekundach")
plt.legend(["Metoda LU"])
plt.savefig("CzasLU")
plt.clf()

tabRozmiarow = [901, 500, 1000, 2000, 3000]


for index in range(0, 0):
    N = tabRozmiarow[index]
    a1 = 3
    a2 = -1
    a3 = -1

    tabA = create_matrixA(N, a1, a2, a3)
    tabB = create_matrixB(N)
    tabDA = create_diagonal_matrix(N, tabA)
    tabLA = create_lower_triangular_matrix(N, tabA)
    tabUA = create_upper_triangular_matrix(N, tabA)
    tabR = [1 for i in range(N)]

    '''
    print("Dla N wynoszącego", N)
    iteracjeJ, czasJ, normaJ = jacobi_method_iterations(tabA, tabB, tabR, tabDA, tabLA, tabUA)
    print("Dla metody Jacobiego błąd rezydualny to", normaJ, "liczba iteracji to", iteracjeJ, "a bylo to w czasie", czasJ)
    iteracjeGS, czasGS, normaGS = gauss_seidel_method_iterations(tabA, tabB, tabR, tabDA, tabLA, tabUA)
    print("Dla metody Gaussa-Seidla błąd rezydualny to", normaGS, "liczba iteracji to", iteracjeGS, "a bylo to w czasie", czasGS)
    '''

    '''
    plt.semilogy(normArrayJ)
    plt.semilogy(normArrayGS)
    plt.title("Błąd rezydualny w kolejnych iteracjach dla podpunktu C")
    plt.xlabel("Numer iteracji")
    plt.ylabel("Wartość błędu rezydualnego")
    plt.legend(["Błąd rezydualny dla metody Jacobiego", "Błąd rezydualny dla metody Gaussa-Seidla"])
    plt.savefig("BladRezydualny")
    plt.clf()
    '''

    c1 = 3
    c2 = -1
    c3 = -1
    tabC = create_matrixA(N, c1, c2, c3)
    tabDC = create_diagonal_matrix(N, tabC)
    tabLC = create_lower_triangular_matrix(N, tabC)
    tabUC = create_upper_triangular_matrix(N, tabC)


    residualLU, czasLU = lu_method_residual(tabC, tabB)
    print("Dla faktoryzacji LU dla przypadku C błąd rezydualny", residualLU, "został on obliczony w czasie", czasLU)
    residualLU, czasLU = lu_method_residual(tabA, tabB)
    print("Dla faktoryzacji LU dla przypadku A błąd rezydualny", residualLU, "został on obliczony w czasie", czasLU)

    print()