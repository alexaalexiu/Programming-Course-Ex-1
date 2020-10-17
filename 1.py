import numpy as np
import math

#calculates eigenvalues of matrix h
def get_evals (h):
    evals, evecs = np.linalg.eig(h)
    return evals

#constructs nxn Huckel matrix for a linear polyene of length n
#a=alpha, b=beta, arbitrarily defined
a=0
b=-1
def huckel_linear (n):
    h=np.ndarray((n, n))
    for i in range(n):
        for j in range(n):
            if (i==j): 
                h[i][j]=a
            elif (j==i+1 or j==i-1):
                h[i][j]=b
            else:
                h[i][j]=0
    return h

#constructs nxn Huckel matrix for a cyclic polyene of length n
def huckel_cyclic (n):
    h=np.ndarray((n, n))
    for i in range(n):
        for j in range(n):
            if (i==j): 
                h[i][j]=a
            elif (j==i+1 or j==i-1 or (j==n-1 and i==0) or (i==n-1 and j==0)):
                h[i][j]=b
            else:
                h[i][j]=0
    return h

#prints the groups of degenerate orbitals
def degeneracies(n,sort_evals):
    i=0
    print("The following orbital pairs are degenerate:")
    while (i<=n-2):
        while ((i<=n-2) and (round(sort_evals[i],3)==round(sort_evals[i+1],3))):
            print("%s, " % str(i+1),end='')
            i+=1
            if (i==n-1):
                print("%s" % str(i+1))
                break
            if(round(sort_evals[i],3)!=round(sort_evals[i+1],3)):
                print("%s" % str(i+1))
                break
        i+=1
    return None

#huckel matrix for tetrahedron            
def huckel_tetrahedron ():
    h=np.ndarray((4, 4))
    for i in range(4):
        for j in range(4):
            if (i==j): 
                h[i][j]=a
            else:
                h[i][j]=b
    return h

#huckel matrix for cube
def huckel_cube ():
    h=np.zeros((8, 8))
    h[0][7]=b
    h[0][3]=b
    h[1][2]=b
    h[1][4]=b
    h[1][7]=b
    h[2][5]=b
    h[2][6]=b
    h[3][6]=b
    h[3][5]=b
    h[4][6]=b
    h[5][7]=b
    h[0][4]=b
    #matrix has to be symmetric
    for i in range(8):
        for j in range(8):
            if (h[i][j]==b):
                h[j][i]=b
    return h

#huckel matrix for dodecahedron
def huckel_dodecahedron ():
    h=np.zeros((20, 20))
    h[0][1]=b
    h[0][2]=b
    h[0][3]=b
    h[1][4]=b
    h[1][5]=b
    h[2][6]=b
    h[2][7]=b
    h[3][8]=b
    h[4][7]=b
    h[5][8]=b
    h[6][9]=b
    h[10][11]=b
    h[12][13]=b
    h[14][15]=b
    h[13][16]=b
    h[14][17]=b
    h[15][18]=b
    h[16][19]=b
    h[17][19]=b
    h[18][19]=b
    i=3
    j=9
    while (i<=12):
        while (j<=18):
            h[i][j]=b
            i+=1
            j+=1
    #matrix has to be symmetric
    for i in range(20):
        for j in range(20):
            if (h[i][j]==b):
                h[j][i]=b
    return h

#huckel matrix for buckminsterfullerene
def huckel_fullerene ():
    h=np.zeros((60, 60))
    for i in range(59):
        h[i][i+1]=b
        h[0][4]=b
        h[0][8]=b
        h[1][11]=b
        h[2][14]=b
        h[3][17]=b
        h[5][19]=b
        h[6][21]=b
        h[7][24]=b
        h[9][25]=b
        h[10][28]=b
        h[12][29]=b
        h[13][32]=b
        h[15][33]=b
        h[16][36]=b
        h[18][37]=b
        h[20][39]=b
        h[22][41]=b
        h[23][43]=b
        h[26][44]=b
        h[27][46]=b
        h[30][47]=b
        h[31][49]=b
        h[34][50]=b
        h[35][52]=b
        h[38][53]=b
        h[40][54]=b
        h[42][56]=b
        h[45][57]=b
        h[48][58]=b
        h[51][59]=b
        h[55][59]=b
        
        #matrix has to be symmetric
        for i in range(60):
            for j in range(60):
                if (h[i][j]==b):
                    h[j][i]=b
    return h

#prints orbital energies in order
def solve_huc(huc,n):
    evals_h=get_evals(huc)
    evals_huc=np.array(evals_h)
    sort_evals=np.sort(evals_huc)
    
    #E=a-bx, x=E for a=0, b=-1
    #requires change of sign for eigenvalues
    for i in range(n):
        if (round(sort_evals[i],3)<0):
            print("%d. \u03B1 +%.3f \u03B2" % ((i+1), abs(sort_evals[i])))
        elif (round(sort_evals[i],3)==0):
            print("%d. \u03B1" % (i+1))
        else:
            print("%d. \u03B1 -%.3f \u03B2" % ((i+1), abs(sort_evals[i])))
    return sort_evals

print("This program determines the Huckel pi-energies and degeneracies for a variety of systems.\n")

option_continue="y"

while (option_continue=="y"):
    print("Please choose one of the available options:")
    print("1 - A linear polyene with n carbons.")
    print("2 - A cyclic polyene with n carbons.")
    print("3 - One of the sp2-hybridized Platonic solids.")
    print("4 - Buckminsterfullerene")
    option=input()

    if (option=="1"):
        print("Please enter the number of atoms in the polyene: ")
        n=input()
        n=int(n)
        huc=huckel_linear(n)
        print("The energies of the molecular orbitals predicted by the Huckel theory are (in order of increasing energy):")
        sort_evals=solve_huc(huc,n)
        print("There are no degenerate orbitals.") 
    
            
    elif (option=="2"):
        print("Please enter the number of atoms in the polyene: ")
        n=input()
        n=int(n)
        huc=huckel_cyclic(n)
        print("The energies of the molecular orbitals predicted by the Huckel theory are (in order of increasing energy):")
        sort_evals=solve_huc(huc,n)
        degeneracies(n,sort_evals)
    
    elif (option=="3"):
        print("Please choose one of the following Platonic solids:")
        print("a - tetrahedron")
        print("b - cube")
        print("c - dodecahedron")
        option2=input()
    
        if (option2=="a"):
            huc=huckel_tetrahedron()
            n=4
        
        elif (option2=="b"):
            huc=huckel_cube()
            n=8
        
        elif (option2=="c"):
            huc=huckel_dodecahedron()
            n=20
        
        else:
            print("The option you selected is invalid.")
        
        print("The energies of the molecular orbitals predicted by the Huckel theory are (in order of increasing energy):")
        sort_evals=solve_huc(huc,n)
        degeneracies(n,sort_evals)

    elif(option=="4"):
        print("The energies of the molecular orbitals predicted by the Huckel theory are (in order of increasing energy):")
        huc=huckel_fullerene()
        n=60
        sort_evals=solve_huc(huc,n)
        degeneracies(n,sort_evals)
    
    else: 
        print("The option you selected is invalid.")
    
    print("\nDo you want to continue? y/n")
    option_continue=input()


