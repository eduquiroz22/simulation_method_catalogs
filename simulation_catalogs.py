import numpy as np
import pandas as pd
from random import uniform
from datetime import datetime
from functools import partial
from multiprocessing import Pool

n_proc = 4
dist = 2
max_dif = 5/3600


# __________________ for data simulation __________________

def cosmic_coordenate(n):		 	                     #(sin(0) ,  sin(30.85)) ~10,575 square degrees)
    return [[uniform(0,360), np.degrees(np.arcsin(uniform(0,0.5127922537212052)))] for i in range(n)]

#def cosmic_coordenate(n):		 	                     #(sin(0) ,  sin(67.34)) ~10,575 square degrees)
#    return [[uniform(0,200), np.degrees(np.arcsin(uniform(0,0.9228072778591087)))] for i in range(n)]


# __________________ for binary search __________________

def calculateDistance(RA_1,RA_2,Dec_1,Dec_2):
    dRA = abs(RA_2-RA_1)
    y = np.sqrt((np.cos(np.radians(Dec_2))*np.sin(np.radians(dRA)))**2+(np.cos(np.radians(Dec_1))*np.sin(np.radians(Dec_2))-np.sin(np.radians(Dec_1))*np.cos(np.radians(Dec_2))*np.cos(np.radians(dRA)))**2)
    x = np.sin(np.radians(Dec_1))*np.sin(np.radians(Dec_2))+np.cos(np.radians(Dec_1))*np.cos(np.radians(Dec_2))*np.cos(np.radians(dRA))
    dist = np.arctan2(y,x)
    dist = np.degrees(dist)*3600.
    return dist

# this function returns a list of 2 values or a empty list (2 values of the start and the end of the range that there is a coincidence)
def binary_search(arr, x, dif=max_dif, mid=0, low=0, high=0):
    if high==0:
        high=len(arr)-1
    result = []
 
    while low <= high:
        mid = round((high + low) / 2)
 
        # If x is at the left of the range, ignore left half
        if arr[mid]-x < -dif:
            low = mid + 1
 
        # If x is at the right of the range, ignore right half
        elif arr[mid]-x > dif:
            high = mid - 1
 
        # means x is between the range
        else:
            if mid>0:
                c=1
                while abs(arr[mid-c]-x)<=dif:
                    c+=1
                    if (mid-c)<0:
                        break
                result.append(mid-c+1)
            else:
                result.append(mid)
            if mid<len(arr)-1:
                c = 1
                while abs(arr[mid+c]-x)<=dif:
                    c+=1
                    if (mid+c)>=len(arr):
                        break
                result.append(mid+c-1)
            else:
                result.append(mid)
            return result
 
    # If we reach here, then (x-dif, x+dif) was not present
    return result


def diference_ra(ra1, ra2):
    dif = abs(ra2-ra1)
    if dif>180:
        dif = 360-dif
    return dif


def binary_search_ra(arr, x, dif=max_dif, mid=0, low=0, high=0):
    if high==0:
        high=len(arr)-1
    result = []
 
    while low <= high:
        mid = round((high + low) / 2)
 
        # If x is greater, ignore left half
        if arr[mid]-x < -dif:
            low = mid + 1
 
        # If x is smaller, ignore right half
        elif arr[mid]-x > dif:
            high = mid - 1
 
        # means x is present at mid
        else:
            c=1
            while diference_ra(arr[mid-c],x)<=dif:
                c+=1
            result.append(mid-c+1)

            l = len(arr)
            c = 1
            while diference_ra(arr[(mid+c)%l],x)<=dif:
                c+=1
            result.append((mid+c-1)%l)
            return result
            
    return result


def get_max_dif_ra(x, max_dif):
    x = np.radians(x)
    return np.degrees(np.arccos((np.cos(x))**(-2)*np.cos(np.radians(max_dif))-(np.tan(x))**2))


def ra_iteration(result, l):
    ini = result[0]
    fin = result[1]

    if ini<=fin:
        iteration = [i for i in range(ini,fin+1)]
    else:
        arr_it = [i for i in range(ini,l)] + [i for i in range(0,ini+1)]
    return iteration


def comparar(vlass, gaia): 
    matches = [] # where we store the coincidences
    for n in range(len(gaia[2])):
        arr = vlass[2] # array of vlass declinations
        dec = gaia[2][n] # value to compare
        
        # getting the coincidences in declination
        resultDec = binary_search(arr, dec)

        if len(resultDec)!=0:
            # creating a copy of vlass value in the range found previously for the declination
            vlass_copy = [e[resultDec[0]:resultDec[-1]+1].copy() for e in vlass]
            order = np.argsort(vlass_copy[1])
            # sorting this data by right ascension ("ra") values
            vlass_copy = [e[order] for e in vlass_copy]
            arr = vlass_copy[1]  # array of vlass "ra" in the range found
            ra = gaia[1][n] # value to compare

            max_dif_ra = get_max_dif_ra(abs(dec)+max_dif,max_dif)
            # if ra exists between the minimum and maximum value of arr
            if (ra<max_dif_ra or ra>(360-max_dif_ra)) or ((arr[0]-max_dif_ra)<ra and ra<(arr[-1]+max_dif_ra)):
                # getting the coincidences in right ascension
                result = binary_search_ra(arr, ra, max_dif_ra)
                if len(result)!=0:
                    for i in ra_iteration(result, len(arr)):
                        distance = calculateDistance(ra, arr[i], dec, vlass_copy[2][i])
                        if distance<dist:
                            #print(n, "gaia: ", gaia[0][n], " vlass: ", vlass[0][i])     # to show in the console the coincidences
                            matches.append([gaia[0][n],      ra,      dec, vlass_copy[0][i], vlass_copy[1][i], vlass_copy[2][i], distance])
                            #matches.append(gaia index, gaia ra, gaia dec,      vlass index,         vlass ra,        vlass dec, distance)
    return matches

def simular(n, i):
    if n%100==0:
        print("Simulacion ",n+1)
    # Generando posiciones de fuentes de radio
    vlass1_ = cosmic_coordenate(n_vlass1)
    vlass1 = [vlass1_index,  np.array(vlass1_)[:,0],  np.array(vlass1_)[:,1]]
    order_Dec = np.argsort(vlass1[2])
    vlass1 = [e[order_Dec] for e in vlass1]

    # Generando de posiciones de estrellas
    gaia_index = np.arange(info_dist['Numero de estrellas'][i])

    gaia_ = cosmic_coordenate(info_dist['Numero de estrellas'][i])
    gaia = [gaia_index,  np.array(gaia_)[:,0],  np.array(gaia_)[:,1]]
    order_Dec = np.argsort(gaia[2])
    gaia = [e[order_Dec] for e in gaia]

    matches1 = comparar(vlass1, gaia)

    n_matches1 = len(matches1)

    return n_matches1

n_vlass1 = 946432
vlass1_index = np.arange(n_vlass1)


info_dist = pd.read_csv('simulaciones_first_250.csv')
fr_1 = [] # final result epoch 1

if __name__=='__main__':
    start = datetime.now()
    times = np.arange(1000)
    for i in range(len(info_dist)):
        start = datetime.now()
        print ("Para distancia de ", info_dist['Distancias'][i], " pc, empezando simulacion")
        with Pool(processes=n_proc) as pool:
            results = np.array(pool.map(partial(simular, i=i), times))
        fr_1.append(results)
        end = datetime.now()
        print(f'tiempo para distancia de {info_dist["Distancias"][i]} pc: {end-start}')

    df1 = pd.DataFrame(np.matrix(fr_1).T, columns = [f'{info_dist["Distancias"][i]}pc' for i in range(len(info_dist))])
    df1.to_csv(f'ss_first_parallel_1000_{dist}arcsec_new_area_250pc.csv',index=False)
