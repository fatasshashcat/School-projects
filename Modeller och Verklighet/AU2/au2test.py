import numpy as np
import matplotlib.pyplot as plt
from tkinter import *

# Tidssteg och maxtid
dt = 0.000001
tmax = 0.1
t = np.arange(0, tmax, dt)
dim = len(t)



def trans_func_spole(w, L, C, R):
    H_L = (w * L) / np.sqrt(R ** 2 + ((1 / (w * C)) - w * L) ** 2)
    return H_L

def trans_func_induc(w, L, C, R):
    H_C = (1/w*C) / np.sqrt(R ** 2 + ((1 / (w * C)) - w * L) ** 2)
    return H_C

def trans_func_res(w, L, C, R):
    H_R = R/np.sqrt(R ** 2 + ((1 / (w * C)) - w * L) ** 2)
    return H_R

# Funktion för att plotta lösningar med val av parametrar
def plot_solutions(L, C, R, U0):
    # Uppdatera lösningar med nya parametrar
    #U0 = 1 # Spänning vid t=0
    f = 127
    global w 
    w = 2 * np.pi * f
    fi = np.arctan(((1 / (w * C)) - (w * L)) / R)
    Z = np.sqrt((R ** 2) + ((1 / (w * C)) - (w * L)) ** 2)
    U = U0 * np.sin(w * t)
    Up = w * U0 * np.cos(w * t)
    xa = U0 * np.sin(w * t + fi)

    Ip = np.zeros(dim)
    I = np.zeros(dim)
    a = np.zeros(dim)
    vz = np.zeros(dim)
    u = np.zeros(dim)
    up = np.zeros(dim)
    I[0] = U0 / Z # I(t)
    Ip[0] = -100 # I'(t)
    a[0] = (Up[0] - (1 / C) * Ip[0] - R * I[0]) / L # I''(t)

    for i in range(dim - 1):
        I[i + 1] = I[i] + Ip[i] * dt # I(t)
        Ip[i + 1] = Ip[i] + a[i] * dt # I'(t)
        a[i + 1] = (Up[i + 1] - (1 / C) * I[i + 1] - R * Ip[i + 1]) / L # I''(t) 
        vz[i] = I[i] * Z 
        u[i] = U0 * np.sin(w * dt * i) 
        up[i] = w * U0 * np.cos(w * dt * (i + 1)) 

    # Plotta lösningarna
    plt.figure(figsize=(10, 6))
    plt.plot(t, vz, label='Numerisk lösning', linestyle='--', color='black')
    plt.plot(t, xa, label='Analytisk lösning', linestyle='-', color='red')
    plt.plot(t, u, label='U(t)', linestyle='-', color='blue')
    plt.title('Numerisk och analytisk lösning av differentialekvation')
    plt.xlabel('Tid')
    plt.ylabel('Position')
    plt.legend()
    plt.grid(True)
    plt.show()

# Funktion för att plotta överföringsfunktionerna
def plot_transfer_functions(L, C, R):
    w_values = np.linspace(0, 1000, 1000) # Frekvensvärden att beräkna överföringsfunktionen för
    H_L_values = [trans_func_spole(w, L, C, R) for w in w_values]
    H_C_values = [trans_func_induc(w, L, C, R) for w in w_values]
    H_R_values = [trans_func_res(w, L, C, R) for w in w_values]

    plt.figure(figsize=(10, 6))
    plt.plot(w_values/(1/np.sqrt(L*C)), H_L_values, label='Överföringsfunktion för spolen', linestyle='--', color='green')
    plt.title('Överföringsfunktion för spolen')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(w_values/(1/np.sqrt(L*C)), H_C_values, label='Överföringsfunktion för kondensatorn', linestyle='--', color='blue')
    plt.title('Överföringsfunktion för kondensatorn')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(w_values/(1/np.sqrt(L*C)), H_R_values, label='Överföringsfunktion för resistorn', linestyle='--', color='blue')
    plt.title('Överföringsfunktion för Resistorn')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()

#FUNKTION FÖR FILTER
def plot_filter_function(L, R, C):
    #R = 30   
    #C = 0.000065
    #L = 0.053
    F = 127
    w = float(2*np.pi*F)    
    w0 = 1 / np.sqrt(L*C)
    w0 = int(w0)
    dim = 3*w0
    H_R = np.zeros(dim)
    H_C = np.zeros(dim)
    H_L = np.zeros(dim)

    w_array = np.arange(1,dim+1,1)

    compare_value = 0
    optimerat_R = 0
    for r in range(1,200,1):
        print_me = 2
        R = r
        färdig = True
        H_L[0] = 0
        for i in range(1, dim, 1):
            H_L[i] = (i * L) / np.sqrt(R ** 2 + ((1 / (i * C)) - i * L) ** 2)

        for x in range(1,dim,10):
            if (H_L[x] > 1.05) or (H_L[x-1] > H_L[x]):
                färdig = False
                break
        
        if färdig and H_L[w0] > compare_value:
            optimerat_HR = r 
            compare_value = H_L[w0]

    print(f"RESISTOR VÄRDE: {optimerat_HR}")

    R = 35#optimerat_HR
    H_L[0] = 0
    for i in range(1, dim, 1):
        H_L[i] = (i * L) / np.sqrt(R ** 2 + ((1 / (i * C)) - i * L) ** 2)

    plt.figure(figsize=(10, 6))
    plt.plot(w_array/(1/np.sqrt(L*C)), H_L, label='Spänning över spolen med 35ohm resistans', linestyle='--', color='green')
    plt.title('Högpass-filter')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_LowPass_filter_function(L, R, C):
    F = 127
    w = float(2*np.pi*F)    
    w0 = 1 / np.sqrt(L*C)
    w0 = int(w0)
    dim = 3*w0
    print(L,R,C)

    #H_R = np.zeros(dim)
    H_C = np.zeros(dim)
    #H_L = np.zeros(dim)

    w_array = np.arange(1,dim+1,1)

    compare_value = 0
    optimerat_LR = 0

    print(L,R,C)
    for r in range(1,300,1):
        print(f"VI ÄR I LOOPEN PÅ: {r}")
        R = r
        färdig = True
        H_C[0] = 1
        for i in range(1, dim, 1):
            H_C[i] = (1/(i*C))/np.sqrt(R**2 + ((1/(i*C)) - i*L)**2)

        for a in range(1,dim,10):
            if (H_C[a] > 1.05) or (H_C[a-1] < H_C[a]):
                print(f"{R} klarade sig inte")
                färdig = False
                break
        
        if färdig and H_C[w0] > compare_value:
            optimerat_LR = r 
            compare_value = H_C[w0]

    print(L,R,C)
    print(f"RESISTOR VÄRDE: {optimerat_LR}")

    R = optimerat_LR
    H_C[0] = 1
    for i in range(1, dim, 1):
        H_C[i] = (1/(i*C))/np.sqrt(R**2 + ((1/(i*C)) - i*L)**2)
    

    plt.figure(figsize=(10, 6))
    plt.plot(w_array/(1/np.sqrt(L*C)), H_C, label='Spänning över kondensatorn med 40ohms resistans', linestyle='--', color='blue')
    plt.title('LågPass-Filter')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()


def Calc_R():
    f_egen = 1 / (2*np.pi*np.sqrt(L*C))

    R_new = 1 / (2*np.pi * f_egen * np.sqrt(L*C))

    R_tot = 1 / (1 / R + 1 / R_new)

    R_tot = 1 / (1/R + 1/R_new)

    return R_tot





# Funktion för att visa dialogruta och plotta baserat på användarens val
def show_plot():
    # Skapa en Tkinter-fönster
    root = Tk()
    root.title("Plotta grafer")

    # Skapa en etikett för parametrar
    param_label = Label(root, text="Ange värden för L, C och R:")
    param_label.pack()

    # Skapa inmatningsfält för parametrar med förinställda standardvärden
    L_label = Label(root, text="L:")
    L_label.pack()
    L_entry = Entry(root)
    L_entry.insert(0, '0.053') # Förinställt standardvärde för L
    L_entry.pack()

    C_label = Label(root, text="C:")
    C_label.pack()
    C_entry = Entry(root)
    C_entry.insert(0, '0.000065') # Förinställt standardvärde för C
    C_entry.pack()

    R_label = Label(root, text="R:")
    R_label.pack()
    R_entry = Entry(root)
    R_entry.insert(0, '10') # Förinställt standardvärde för R
    R_entry.pack()

    U0_label = Label(root, text="U:")
    U0_label.pack()
    U0_entry = Entry(root)
    U0_entry.insert(0, '1') # Förinställt standardvärde för L
    U0_entry.pack()


    # Skapa en knapp för att plotta lösningar med angivna parametrar
    plot_button = Button(root, text="Plotta Lösningar", command=lambda: plot_solutions(float(L_entry.get()), float(C_entry.get()), float(R_entry.get()), float(U0_entry.get())))
    plot_button.pack()

    # Skapa en knapp för att plotta överföringsfunktioner
    transfer_functions_button = Button(root, text="Plotta överföringsfunktioner", command=lambda: plot_transfer_functions(float(L_entry.get()), float(C_entry.get()), float(R_entry.get())))
    transfer_functions_button.pack()

    filter_button = Button(root, text ="Visa Högpass-filter", command=lambda: plot_filter_function(float(L_entry.get()), 35 ,float(C_entry.get())))
    filter_button.pack()

    filter_button = Button(root, text ="Visa lågpass-Filter", command=lambda: plot_LowPass_filter_function(float(L_entry.get()), 35 ,float(C_entry.get())))
    filter_button.pack()

    # Visa fönstret
    root.mainloop()

# Visa dialogruta och plotta baserat på användarens val
show_plot()

'''
def plot_överföringsfunktion(print_me):
    fig, ax = plt.subplots()
    ax.set_ylim(0, 2)
    ax.set_xlim(0,3)

    if print_me == 1:
        ax.plot(w_array/w0, H_C, label = 'H_C')
    if print_me == 2:
        ax.plot(w_array/w0, H_L, label = 'H_L')
    if print_me == 3:
        ax.plot(w_array/w0, H_R, label = 'H_R')

    ax.set_ylabel(r'$H(w)$',fontsize=14)
    ax.set_xlabel(r'$w/w0$',fontsize=14)
    ax.legend()
    ax.tick_params(labelsize=14)
    
    plt.show()

#H_C = (1/w*C) / np.sqrt(R ** 2 + ((1 / (w * C)) - w * L) ** 2)
R = 30   
C = 0.000065
L = 0.053
F = 127
w = float(2*np.pi*F)    
w0 = 1 / np.sqrt(L*C)
w0 = int(w0)
dim = 3*w0

H_R = np.zeros(dim)
H_C = np.zeros(dim)
H_L = np.zeros(dim)

w_array = np.arange(1,dim+1,1)

compare_value = 0
optimerat_R = 0

#--------------------------------LÅGPASS--------------------------------

print_me = 1
for r in range(1,200,1):
    R = r
    färdig = True
    H_C[0] = 1
    for i in range(1, dim, 1):
        H_C[i] = (1/(i*C))/np.sqrt(R**2 + ((1/(i*C)) - i*L)**2)

    for x in range(1,dim,10):
        if (H_C[x] > 1.05) or (H_C[x-1] < H_C[x]):
            färdig = False
            break
    
    if färdig and H_C[w0] > compare_value:
        optimerat_LR = r 
        compare_value = H_C[w0]

print(f"RESISTOR VÄRDE: {optimerat_LR}")

R = optimerat_LR
H_C[0] = 1
for i in range(1, dim, 1):
    H_C[i] = (1/(i*C))/np.sqrt(R**2 + ((1/(i*C)) - i*L)**2)

#--------------------------------HÖGPASS--------------------------------

for r in range(1,200,1):
    print_me = 2
    R = r
    färdig = True
    H_L[0] = 0
    for i in range(1, dim, 1):
        H_L[i] = (i * L) / np.sqrt(R ** 2 + ((1 / (i * C)) - i * L) ** 2)

    for x in range(1,dim,10):
        if (H_L[x] > 1.05) or (H_L[x-1] > H_L[x]):
            färdig = False
            break
    
    if färdig and H_L[w0] > compare_value:
        optimerat_HR = r 
        compare_value = H_L[w0]

print(f"RESISTOR VÄRDE: {optimerat_HR}")

R = optimerat_HR
H_L[0] = 0
for i in range(1, dim, 1):
    H_L[i] = (i * L) / np.sqrt(R ** 2 + ((1 / (i * C)) - i * L) ** 2)
'''
'''
#--------------------------------BANDPASS--------------------------------
print_me = 3

for r in range(1,200,1):
    print_me = 2
    R = r
    färdig = False
    H_R[0] = 0

    for i in range(1, dim, 1):
        H_R[i] = R / np.sqrt(R ** 2 + ((1 / (i * C)) - i * L) ** 2)
    
    if(0.95 < H_R[w0] < 1.05) and (H_R[2*w0] <= 0.5):
        färdig = True

    if färdig:
        if(H_R[2*w0] > compare_value):
            compare_value = H_R[2*w0]
            optimerat_R = R
 
print(f"RESISTOR VÄRDE: {optimerat_R}")

R = optimerat_R
H_R[0] = 0
for i in range(1, dim, 1):
    H_R[i] = R / np.sqrt(R ** 2 + ((1 / (i * C)) - i * L) ** 2)
    

plot_överföringsfunktion(print_me)
'''