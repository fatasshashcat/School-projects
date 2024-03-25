import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import messagebox

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
    H_R =  R/np.sqrt(R ** 2 + ((1 / (w * C)) - w * L) ** 2)
    return H_R

# Funktion för att plotta lösningar med val av parametrar
def plot_solutions(L, C, R):
    # Uppdatera lösningar med nya parametrar
    U0 = 1  # Spänning vid t=0
    f = 127
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
    I[0] = U0 / Z  # I(t)
    Ip[0] = -100   # I'(t)
    a[0] = (Up[0] - (1 / C) * Ip[0] - R * I[0]) / L  # I''(t)

    for i in range(dim - 1):
        I[i + 1] = I[i] + Ip[i] * dt             # I(t)
        Ip[i + 1] = Ip[i] + a[i] * dt             # I'(t)
        a[i + 1] = (Up[i + 1] - (1 / C) * I[i + 1] - R * Ip[i + 1]) / L  # I''(t) 
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
def plot_transfer_functions():
    # Standardvärden för parametrarna
    L = 0.053
    C = 0.000065
    R = 10

    w_values = np.linspace(1, 1000, 1000)  # Frekvensvärden att beräkna överföringsfunktionen för
    H_L_values = [trans_func_spole(w, L, C, R) for w in w_values]
    H_C_values = [trans_func_induc(w, L, C, R) for w in w_values]
    H_R_values = [trans_func_res(w, L, C, R) for w in w_values]

    plt.figure(figsize=(10, 6))
    plt.plot(w_values, H_L_values, label='Överföringsfunktion för spolen', linestyle='--', color='green')
    plt.title('Överföringsfunktion för spolen')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(w_values, H_C_values, label='Överföringsfunktion för kondensatorn', linestyle='--', color='blue')
    plt.title('Överföringsfunktion för kondensatorn')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(w_values, H_R_values, label='Överföringsfunktion för resistorn', linestyle='--', color='blue')
    plt.title('Överföringsfunktion för Resistorn')
    plt.xlabel('Frekvens (Hz)')
    plt.ylabel('Amplitud')
    plt.legend()
    plt.grid(True)
    plt.show()

# Funktion för att visa dialogruta och plotta baserat på användarens val
def show_plot():
    # Skapa en Tkinter-fönster
    root = Tk()
    root.title("Plotta grafer")

    # Funktion för att hantera valda parametrar och plotta lösningar
    def plot_solutions_with_params():
        L = float(L_entry.get())
        C = float(C_entry.get())
        R = float(R_entry.get())
        plot_solutions(L, C, R)

    # Skapa en etikett för parametrar
    param_label = Label(root, text="Ange värden för L, C och R:")
    param_label.pack()

    # Skapa inmatningsfält för parametrar med förinställda standardvärden
    L_label = Label(root, text="L:")
    L_label.pack()
    L_entry = Entry(root)
    L_entry.insert(0, '0.053')  # Förinställt standardvärde för L
    L_entry.pack()

    C_label = Label(root, text="C:")
    C_label.pack()
    C_entry = Entry(root)
    C_entry.insert(0, '0.000065')  # Förinställt standardvärde för C
    C_entry.pack()

    R_label = Label(root, text="R:")
    R_label.pack()
    R_entry = Entry(root)
    R_entry.insert(0, '10')  # Förinställt standardvärde för R
    R_entry.pack()

    # Skapa en knapp för att plotta lösningar med angivna parametrar
    plot_button = Button(root, text="Plotta lösningar", command=plot_solutions_with_params)
    plot_button.pack()

    # Skapa en knapp för att plotta överföringsfunktioner
    transfer_functions_button = Button(root, text="Plotta överföringsfunktioner", command=plot_transfer_functions)
    transfer_functions_button.pack()

    # Visa fönstret
    root.mainloop()

# Visa dialogruta och plotta baserat på användarens val
show_plot()
