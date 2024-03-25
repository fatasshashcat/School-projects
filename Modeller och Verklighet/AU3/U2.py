import tkinter as tk
import numpy as np

tot_eng = 0
tot_h = 0

def beräknaeffekt(phi, n, t, I_0, solalt_p, azimut_p, epsilon, A):
    deklination = -23.44 * np.cos(2 * np.pi / 365 * n)
    timvinkel = 15 * t - 180
    timvinkel = timvinkel if timvinkel != 0 else 0.01
    solalt = np.degrees(np.arcsin(np.sin(np.radians(phi)) * np.sin(np.radians(deklination)) + np.cos(np.radians(phi)) * np.cos(np.radians(deklination)) * np.cos(np.radians(timvinkel))))
    
    if timvinkel > 0:
        az_east = np.degrees(np.arccos(((np.sin(np.radians(phi)) * np.sin(np.radians(solalt))) - (np.sin(np.radians(deklination)))) / (np.cos(np.radians(phi)) * np.cos(np.radians(solalt)))))
        azimut = 180 - az_east
    else: 
        az_west = np.degrees(np.arccos(((np.sin(np.radians(phi)) * np.sin(np.radians(solalt))) - (np.sin(np.radians(deklination)))) / ((np.cos(np.radians(phi)) * np.cos(np.radians(solalt))))))
        azimut = 180 + az_west

    I = 1.1 * I_0 * 0.7**(1/np.sin(np.radians(solalt))**(0.678))
    I_p = I * (np.cos(np.radians(solalt_p - solalt)) * np.cos(np.radians(azimut_p - azimut)) + (1 - np.cos(np.radians(azimut_p - azimut))) * np.sin(np.radians(solalt_p)) * np.sin(np.radians(solalt)))
    P = epsilon * I_p * A 
    P = P if P >= 0 else 0
    return P

def beräknaeffektår(n):
    global tot_eng,tot_h
    tot_eng = 0
    tot_år = 0
    n=n
    for dag in range(1,31):
        tot_dag, tot_h_dag = beräknaeffektdag(dag, n)
        tot_mån += tot_dag
        tot_h += tot_h_dag
    n+31
    for dag in range(32,60):
        tot_dag, tot_h_dag = beräknaeffektdag(dag, n)
        tot_mån += tot_dag
        tot_h += tot_h_dag
    n+31 

    return tot_år


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
        tot_dag += power
        if power > 0:
            tot_h_dag += 1
        tot_eng += power  
    return tot_dag, tot_h_dag

def calculate_energy(month):
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
        energy += beräknaeffektmånad(32)
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

        print("Effekt hela året:", int(energy/1000), "kWh")

        sum = 0
        for i in soltimmelista:
            sum = sum + i
        print("Effekt per år soltimmar:", int((energy * (sum/tot_h_år))/1000), "kWh")
    # Visa beräknad energi i en messagebox
    #tk.messagebox.showinfo("Beräknad Energi", energy)

phi = 55.6
n = 4
t = 14.5
I_0 = 1360
solalt_p = 20
azimut_p = 180
epsilon = 0.15
A = 30

soltimmelista  = [32,74,142,196,237,243,244,199,132,81,38,22]

hours = np.arange(0, 24, 1)
powers = []


# Skapa ett fönster
root = tk.Tk()
root.title("Energi Beräknare")

# Skapa en etikett för titeln
title_label = tk.Label(root, text="Beräkna energi för")
title_label.pack()

# Funktion för att skapa knappar för varje månad
def create_button(month):
    return tk.Button(root, text=month, command=lambda: calculate_energy(month))

# Skapa knappar för varje månad och hela året
january_button = create_button("Januari")
january_button.pack()

july_button = create_button("Juli")
july_button.pack()

whole_year_button = create_button("Hela året")
whole_year_button.pack()

# Starta huvudloopen för GUI:et
root.mainloop()
