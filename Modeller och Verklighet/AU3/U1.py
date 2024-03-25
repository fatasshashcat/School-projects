import matplotlib.pyplot as plt
import numpy as np

def uppg1(): 
    def beräknaeffekt(phi, n, t, I_0, solalt_p, azimut_p, epsilon, A):
        deklination = -23.44 * np.cos(2 * np.pi / 365 * n)
        timvinkel = 15 * t - 180
        timvinkel = timvinkel if timvinkel != 0 else 0.01
        solalt = np.degrees(np.arcsin(np.sin(np.radians(phi)) * np.sin(np.radians(deklination)) + np.cos(np.radians(phi)) * np.cos(np.radians(deklination)) * np.cos(np.radians(timvinkel))))
        soltimmelista  = [32,74,142,196,237,243,244,199,132,81,38,22]

        if timvinkel > 0:
            az_east = np.degrees(np.arccos(((np.sin(np.radians(phi)) * np.sin(np.radians(solalt))) - (np.sin(np.radians(deklination)))) / (np.cos(np.radians(phi)) * np.cos(np.radians(solalt)))))
            azimut = 180 - az_east
            #print("East:", az_east)
        else: 
            az_west = np.degrees(np.arccos(((np.sin(np.radians(phi)) * np.sin(np.radians(solalt))) - (np.sin(np.radians(deklination)))) / ((np.cos(np.radians(phi)) * np.cos(np.radians(solalt))))))
            azimut = 180 + az_west
            #print("West:", az_west)

        I = 1.1 * I_0 * 0.7**(1/np.sin(np.radians(solalt))**0.678)
        I_p = I * (np.cos(np.radians(solalt_p - solalt)) * np.cos(np.radians(azimut_p - azimut)) + (1 - np.cos(np.radians(azimut_p - azimut))) * np.sin(np.radians(solalt_p)) * np.sin(np.radians(solalt)))
        P = epsilon * I_p * A 
        P = P if P >= 0 else 0

        return P

    phi = 55.6
    n = 183
    t = 14.5
    I_0 = 1360
    solalt_p = 20
    azimut_p = 180
    epsilon = 0.15
    A = 30

    hours = np.arange(0, 24, 1)
    powers = []

    solalt = beräknaeffekt(phi, n, t, I_0, solalt_p, azimut_p, epsilon, A)
    print("Solens altitud:", solalt)

    for hour in hours:
        power = beräknaeffekt(phi, n, hour, I_0, solalt_p, azimut_p, epsilon, A)
        powers.append(power)

    # Plotta effekten över timmarna
    plt.plot(hours, powers)
    plt.title('Effekt över en dag')
    plt.xlabel('Timme på dagen')
    plt.ylabel('Effekt (W)')
    plt.grid(True)
    plt.show()
    
uppg1()
