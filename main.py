import tkinter as tk
from tkinter import messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import csv

# Constantes físicas
h = 6.62607015e-34
c = 2.99792458e8
e = 1.60217662e-19
m_e = 9.1093837e-31
epsilon_0 = 8.854187817e-12
pi = np.pi
Ry = 13.59844

# Límites de validación
MAX_Z = 118
MAX_N = 50

class BohrModel:
    def __init__(self, Z, n):
        self.Z = Z
        self.n = n

    def radio(self):
        return (self.n**2 * h**2 * epsilon_0) / (pi * m_e * e**2 * self.Z)

    def energia(self, relativista=False):
        E = -Ry * (self.Z**2) / (self.n**2)
        if relativista:
            alpha = 1/137.035999
            E *= (1 + (alpha**2 * self.Z**2) / (self.n**2))
        return E

    def velocidad(self):
        return (self.Z * e**2) / (2 * epsilon_0 * h * self.n)

    @staticmethod
    def transicion(Z, n1, n2):
        if n1 == n2:
            raise ValueError("n1 y n2 deben ser diferentes")
        delta_E = Ry * Z**2 * (1/n2**2 - 1/n1**2) if n1 > n2 else -Ry * Z**2 * (1/n1**2 - 1/n2**2)
        freq = abs(delta_E) * e / h
        lambda_onda = c / freq if freq != 0 else float('inf')
        serie = BohrModel.clasificar_serie(n1, n2)
        return delta_E, freq, lambda_onda * 1e9, serie

    @staticmethod
    def clasificar_serie(n1, n2):
        min_n, max_n = sorted([n1, n2])
        if min_n == 1: return "Lyman (UV)"
        elif min_n == 2: return "Balmer (Visible)"
        elif min_n == 3: return "Paschen (IR)"
        elif min_n == 4: return "Brackett (IR)"
        else: return "Serie general"

class App:
    def __init__(self, root):
        self.root = root
        root.title("Calculadora Modelo Atómico de Bohr - Kevin Julián Gómez Oliveros y Alejandro Diaz")

        tk.Label(root, text="Número atómico Z (1-118):").grid(row=0, column=0)
        self.Z_entry = tk.Entry(root)
        self.Z_entry.grid(row=0, column=1)

        tk.Label(root, text="Nivel n inicial (1-50):").grid(row=1, column=0)
        self.n1_entry = tk.Entry(root)
        self.n1_entry.grid(row=1, column=1)

        tk.Label(root, text="Nivel n final (opcional):").grid(row=2, column=0)
        self.n2_entry = tk.Entry(root)
        self.n2_entry.grid(row=2, column=1)

        # self.relativista_var = tk.BooleanVar()
        # tk.Checkbutton(root, text="Corrección relativista", variable=self.relativista_var).grid(row=3, column=0, columnspan=2)

        tk.Button(root, text="Calcular 🧮", command=self.calcular).grid(row=4, column=0)
        tk.Button(root, text="Graficar 📈", command=self.graficar_niveles).grid(row=4, column=1)
        tk.Button(root, text="Exportar CSV 📥", command=self.exportar).grid(row=4, column=2)
        tk.Button(root, text="Limpiar ⚡", command=self.limpiar).grid(row=4, column=3)
        tk.Button(root, text="Info. ℹ️", command=self.mostrar_info).grid(row=4, column=4)

        self.result_text = tk.Text(root, height=10, width=60)
        self.result_text.grid(row=5, column=0, columnspan=5)

        self.figure = plt.Figure(figsize=(5, 4))
        self.canvas = FigureCanvasTkAgg(self.figure, root)
        self.canvas.get_tk_widget().grid(row=6, column=0, columnspan=5)

    def validar(self):
        Z = int(self.Z_entry.get())
        n1 = int(self.n1_entry.get())
        if not (1 <= Z <= MAX_Z):
            raise ValueError(f"Z debe estar entre 1 y {MAX_Z}")
        if not (1 <= n1 <= MAX_N):
            raise ValueError(f"n debe estar entre 1 y {MAX_N}")
        n2 = self.n2_entry.get()
        n2 = int(n2) if n2 else None
        if n2 and not (1 <= n2 <= MAX_N):
            raise ValueError(f"n final debe estar entre 1 y {MAX_N}")
        return Z, n1, n2

    def calcular(self):
        # Validación básica
        try:
            Z = int(self.Z_entry.get())
            n1 = int(self.n1_entry.get())
            if Z < 1 or n1 < 1:
                raise ValueError("Z y n deben ser enteros positivos")
        except:
            messagebox.showerror("Error", "Ingresa valores válidos y positivos para Z y n")
            return

        try:
            Z, n1, n2 = self.validar()
            model1 = BohrModel(Z, n1)
            resultados = f"Para n={n1}:\nRadio: {model1.radio():.2e} m\nEnergía: {model1.energia(
                #self.relativista_var.get()
            ):.2f} eV\nVelocidad: {model1.velocidad():.2e} m/s\n"
            if n2:
                delta_E, freq, lambda_nm, serie = BohrModel.transicion(Z, n1, n2)
                resultados += f"\nTransición {n1}->{n2}:\nΔE: {delta_E:.2f} eV\nFrecuencia: {freq:.2e} Hz\nLongitud de onda: {lambda_nm:.2f} nm ({serie})\n"
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, resultados)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def graficar_niveles(self):
        try:
            Z, n1, n2 = self.validar()
            max_n = max(n1, n2 or 0) + 5
            ns = np.arange(1, max_n + 1)
            Es = [-Ry * Z**2 / n**2 for n in ns]
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.plot(ns, Es, 'bo-')
            ax.set_xlabel('Nivel n')
            ax.set_ylabel('Energía (eV)')
            ax.set_title('Niveles de Energía Atómico')
            self.canvas.draw()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def exportar(self):
        try:
            Z, n1, _ = self.validar()
            model = BohrModel(Z, n1)
            data = [["Propiedad", "Valor"], ["Radio (m)", model.radio()], ["Energía (eV)", model.energia()], ["Velocidad (m/s)", model.velocidad()]]
            file = filedialog.asksaveasfilename(defaultextension=".csv")
            if file:
                with open(file, 'w', newline='') as f:
                    csv.writer(f).writerows(data)
                messagebox.showinfo("Éxito", "Archivo CSV guardado")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def limpiar(self):
        self.Z_entry.delete(0, tk.END)
        self.n1_entry.delete(0, tk.END)
        self.n2_entry.delete(0, tk.END)
        self.result_text.delete(1.0, tk.END)
        self.figure.clear()
        self.canvas.draw()

    def mostrar_info(self):
        txt = "Este programa utiliza el modelo atómico de Bohr para calcular el radio, la energía y la velocidad de un electrón en un átomo según el nivel n y número atómico Z. También puede mostrar transiciones electrónicas y su longitud de onda."
        messagebox.showinfo("Acerca del programa", txt)

root = tk.Tk()
root.resizable(False, False)
app = App(root)
root.mainloop()
