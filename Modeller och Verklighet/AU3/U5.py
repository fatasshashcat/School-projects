import tkinter as tk
import numpy as np
from scipy.optimize import minimize_scalar

tot_eng = 0
tot_h = 0
P = 0
optimalvinkel = 0

def beräknaeffekt(phi, n, t, I_0, alt_p, azimut_p, epsilon, A):
    deklination = -23.44 * np.cos(2 * np.pi / 365 * n)
    timvinkel = 15 * t - 180
    timvinkel = timvinkel if timvinkel != 0 else 0.01
    solalt = np.degrees(np.arcsin(np.sin(np.radians(phi)) * np.sin(np.radians(deklination)) + np.cos(np.radians(phi)) * np.cos(np.radians(deklination)) * np.cos(np.radians(timvinkel))))
    
    '''
    az_west = np.degrees(np.arccos(((np.sin(np.radians(phi)) * np.sin(np.radians(solalt))) - (np.sin(np.radians(deklination)))) / ((np.cos(np.radians(phi)) * np.cos(np.radians(solalt))))))
    azimut = 180 + az_west
    print("West:", az_west)
    '''
    if timvinkel > 0:
        az_east = np.degrees(np.arccos(((np.sin(np.radians(phi)) * np.sin(np.radians(solalt))) - (np.sin(np.radians(deklination)))) / (np.cos(np.radians(phi)) * np.cos(np.radians(solalt)))))
        azimut = 180 - az_east
        #print("East:", az_east)
    else: 
        az_west = np.degrees(np.arccos(((np.sin(np.radians(phi)) * np.sin(np.radians(solalt))) - (np.sin(np.radians(deklination)))) / ((np.cos(np.radians(phi)) * np.cos(np.radians(solalt))))))
        azimut = 180 + az_west
        #print("West:", az_west)

    I = 1.1 * I_0 * (0.7 ** (1 / (np.sin(np.radians(solalt)) ** 0.678)))
    I_p = I * (np.cos(np.radians(0)) * np.cos(np.radians(0)) + (1 - np.cos(np.radians(0))) * np.sin(np.radians(alt_p)) * np.sin(np.radians(solalt)))
    P = epsilon * I_p * A 
    P = P if P >= 0 else 0
    return P

def beräknaeffektmånad(n):    
    global tot_eng, tot_h  
    tot_eng = 0
    tot_mån = 0
    n = n
    for dag in range(1, 31):
        tot_dag, tot_h_dag = beräknaeffektdag(dag, n)
        tot_mån += tot_dag
        tot_h += tot_h_dag
    return tot_mån

def beräknaeffektdag(dag, n):
    global tot_eng, tot_h  
    tot_dag = 0
    tot_h_dag = 0
    n = n
    for hour in hours:
        power = beräknaeffekt(phi, n + dag, hour, I_0, solalt_p, azimut_p, epsilon, A)
        #print("Effekt för timme {}: {}".format(hour, power))
        tot_dag += power
        if power > 0:
            tot_h_dag += 1
            #tot_h += 1  
        tot_eng += power  
    return tot_dag, tot_h_dag

def calculate_energy(month):
    # Här kan du lägga till din beräkningslogik beroende på vilken månad som valts
    if month == "Januari":
        energy = beräknaeffektmånad(1)
        print("Total energi i Januari: ", int(energy / 1000), "kWh")
        print("Total energi i Januari med soltimmar: ", int((energy / 1000)*(32/210)), "kWh")
    elif month == "Juli":
        energy = beräknaeffektmånad(182)
        print("Total energi i Juni: ", int(energy / 1000), "kWh")
        print("Total energi i Juni med soltimmar: ", int((energy / 1000)*(244/330)), "kWh")
    else:
        tot_h_år = 0
        energy = beräknaeffektmånad(1)
        print("tot_h :", tot_h)
        tot_h_år += tot_h
        #energy = int((energy/1000) * (32/210))
        energy += beräknaeffektmånad(32)
        print("tot_h :", tot_h)
        energy += beräknaeffektmånad(60)
        energy += beräknaeffektmånad(91)
        energy += beräknaeffektmånad(121)
        energy += beräknaeffektmånad(152)
        energy += beräknaeffektmånad(182)
        energy += beräknaeffektmånad(213)
        energy += beräknaeffektmånad(243)
        energy += beräknaeffektmånad(274)
        energy += beräknaeffektmånad(204)
        energy += beräknaeffektmånad(335)

        tot_h_år += tot_h

        print("Antal effekt hela året:", int(energy/1000), "kWh")

        sum = 0
        #print("sum:", sum)
        for i in soltimmelista:
            sum = sum + i
            #print("i:", i)
            #print("sum:", sum)
        #print("tot_h_år", tot_h_år)
        #print("sum:", sum)
        #print("Antal timmar per år:", tot_h_år)
        print("Antal timmar per år soltimmar:", int((energy * (sum/tot_h_år))/1000), "kWh")
        #print("solalt_p: ", solalt_p)
        optimal_solalt = optimize_solalt_p()
        print("Optimalt solalt_p", optimal_solalt)

        return energy

def check_best_vinkel():
    global optimalvinkel  # Declare optimalvinkel as global
    phi = 55.6
    I_0 = 1360
    azimut_p = 180
    epsilon = 0.15
    A = 30
    optimaleffekt = 0
    optimalvinkel = 0  # Initialize optimalvinkel
    
    for k in range(1, 91, 1):
        powerdaily = 0
        tot_h = 0  # Reset tot_h for each iteration of the outer loop
        
        for i in range(1, 366):
            for j in range(0, 24, 1):
                powerhourly = beräknaeffekt(phi, i, j, I_0, k, azimut_p, epsilon, A)
                powerdaily += powerhourly
                
                if powerhourly > 0:
                    tot_h = tot_h + 1

            if powerdaily > optimaleffekt:
                optimaleffekt = powerdaily
                optimalvinkel = k
                
    roi = 200000 / (2 * (optimaleffekt / 1000))
    print("Bästa effekt om alltid soligt:", (optimaleffekt / 1000))
    print("ROI om alltid soligt:", roi, "år")
    
    optimaleffekt = optimaleffekt * (1640 / 3632) / 1000
    roi = 200000 / (2 * optimaleffekt)
    print("Bästa effekt med hänsyn till faktiska soltimmar är:", optimaleffekt, "kWh, det får vi vid vinkeln:", optimalvinkel)
    print("ROI med hänsyn till faktiska soltimmar:", roi, "år")

        
def check_best_vinkelvinkel():
    #global optimalvinkel  # Declare optimalvinkel as global before using it
    phi = 55.6
    I_0 = 1360
    azimut_p = 0
    epsilon = 0.15
    A = 30
    #solalt = optimalvinkel
    optimaleffekt = 0
    optimalvinkel = 0
    
    for k in range(1, 91, 1):
        powerdaily = 0
        tot_h = 0  # Reset tot_h for each iteration of the outer loop
        
        for i in range(1, 366):
            for j in range(0, 24,1):
                powerhourly = beräknaeffekt(phi, i, j, I_0, 45, k, epsilon, A)
                powerdaily += powerhourly
                
                if powerhourly > 0:
                    tot_h = tot_h + 1

            if powerdaily > optimaleffekt:
                optimaleffekt = powerdaily
                optimalvinkel = k
                
    print("")
    print("Bästa effekt om alltid soligt:", (optimaleffekt / 1000))
    roi = 400000/(2*(optimaleffekt/1000))
    print("ROI om alltid soligt:", roi, "år")
    
    optimaleffekt = optimaleffekt * (1640 / 3632) / 1000
    print("Bästa effekt med hänseende på faktiska soltimmar är:", optimaleffekt, "kWh, det får vi vid vinkeln:", optimalvinkel)
    roi = 200000/(2*optimaleffekt)
    print("ROI med hänseende till faktiska soltimmar:", roi, "år")

#antalsoltimmartabell / antal timmar solen upp = solengereffekt
#check_best_vinkel()
check_best_vinkelvinkel()
#1640 / 