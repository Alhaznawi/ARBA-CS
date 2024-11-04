""""
File: main.py
Author: Abdulaziz S. Alghamdi
Date of modification: 03 November 2024
Description: ARPA-CS source code
"""
from cycler import cycler
import time
import serial
import copy
import math
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from tkmacosx import Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
plt.rcParams.update({'font.size': 5})
plt.rcParams.update({'axes.facecolor': 'black'})
cases = [None]
colors = ['#FFFF00']
plt.rcParams['axes.prop_cycle'] = cycler(markevery=cases, color=colors)


# this function performs Walsh-Hadamard Transform on received/input signal
def fwowht(x):
    n = len(x)
    for i in range(0, n-1, 2):
        x[i] = x[i] + x[i+1]
        x[i+1] = x[i] - 2 * x[i+1]
    f = 1
    y1 = np.zeros(n)
    for nStage in range(2, int(math.log2(n))+1):
        m = 2 ** f
        p = 0
        k = 0
        while k < n:
            for j in range(p, p+m, 2):
                y1[k] = x[j] + x[j+m]
                y1[k+1] = x[j] - x[j+m]
                y1[k+2] = x[j+1] - x[j+1+m]
                y1[k+3] = x[j+1] + x[j+1+m]
                k = k + 4
            p = p + 2 * m
        x = copy.deepcopy(y1)
        f = f + 1
    y1 = np.divide(x, n)
    y = y1
    return y


# this function generates Walsh function (reference signal)
def walsh_generator(order, k):
    n = 2 ** k
    if order == 0:
        return np.ones(n)
    else:
        j = np.arange(1, n+1)
        m = 1 + np.floor(math.log(order, 2))
        r = (-1) ** np.floor(2 ** m * ((j - 1)/n))
        return r * walsh_generator(2 ** m-1-order, k)


def myprogram1(w1, select, seq):
    mcu_data1 = serial.Serial('/dev/cu.usbmodem1101', 115200, timeout=2)
    time.sleep(1)

    # we need to re-write this again to prevent program crashing
    if mcu_data1.is_open:

        my_data1 = str(select) + '\r'
        mcu_data1.write(my_data1.encode())  # encode(): encodes string using UTF-8

        for i in range(len(w1)):
            if w1[i] == 1:
                w1[i] = 4096  # maximum value of DAC (12-bit)
                my_data1 = str(w1[i]) + '\r'
                # print(my_data.encode())
                mcu_data1.write(my_data1.encode())  # encode(): encodes string using UTF-8
            else:
                w1[i] = 0
                my_data1 = str(w1[i]) + '\r'
                mcu_data1.write(my_data1.encode())  # encode(): encodes string using UTF-8
    else:
        # print("check COM3")
        summary_label1.config(text="Error: Check COM port and try again")

    # this while loop to waiting
    while mcu_data1.inWaiting() == 0:
        summary_label1.config(text="Status: Waiting...")
        pass

    # this code to read received data
    counter = 0  # do we need this?
    p = 0
    dx = np.zeros(len(w1))
    while counter < 10:
        data1 = mcu_data1.readline().decode('utf-8')
        if data1 == "":
            break
        else:
            dx[p] = int(data1)
            print(dx[p])
            p = p + 1
    # ################################ phase alignment algorithm #######################################################
    mv = np.zeros(len(w1))
    for q in range(len(w1)-1):
        dv = copy.deepcopy(dx)
        dv = np.roll(dv, -1*q)
        fv = fwowht(dv)
        mv[q] = fv[seq]
    index = np.argmax(mv)
    s1 = copy.deepcopy(index)
    # ################################ ######################### #######################################################
    dt = np.roll(dx, -1*index)
    df = np.roll(dx, -1*index)
    # ################################ ######################### #######################################################
    dz = copy.deepcopy(dt)
    fx = fwowht(df)
    fz = copy.deepcopy(fx)
    # print(f)
    return dz, fz, s1


def myprogram2(w2, select, seq):
    mcu_data2 = serial.Serial('/dev/cu.usbmodem1101', 115200, timeout=2)
    time.sleep(1)

    # we need to re-write this again to prevent program crashing
    if mcu_data2.is_open:

        my_data2 = str(select) + '\r'
        mcu_data2.write(my_data2.encode())  # encode(): encodes string using UTF-8

        for i in range(len(w2)):
            if w2[i] == 1:
                w2[i] = 4096  # maximum value of DAC (12-bit)
                my_data2 = str(w2[i]) + '\r'
                mcu_data2.write(my_data2.encode())  # encode(): encodes string using UTF-8
            else:
                w2[i] = 0
                my_data = str(w2[i]) + '\r'
                mcu_data2.write(my_data.encode())  # encode(): encodes string using UTF-8
    else:
        # print("check COM3")
        summary_label1.config(text="Error: Check COM port and try again")

    # this while loop to waiting
    while mcu_data2.inWaiting() == 0:
        summary_label1.config(text="Status: Waiting...")
        pass

    # this code to read received data
    counter = 0  # do we need this?
    p = 0
    dy = np.zeros(len(w2))
    while counter < 10:
        data2 = mcu_data2.readline().decode('utf-8')
        if data2 == "":
            break
        else:
            dy[p] = int(data2)
            p = p + 1
    # ################################ phase alignment algorithm #######################################################
    mr = np.zeros(len(w2))
    for k in range(len(w2) - 1):
        dr = copy.deepcopy(dy)
        dr = np.roll(dr, -1 * k)
        fr = fwowht(dr)
        mr[k] = fr[seq]
    index2 = np.argmax(mr)
    s2 = copy.deepcopy(index2)
    # ################################ ######################### #######################################################
    du = np.roll(dy, -1 * index2)
    dm = np.roll(dy, -1 * index2)
    # ################################ ######################### #######################################################
    dq = copy.deepcopy(du)
    fy = fwowht(dm)
    fq = copy.deepcopy(fy)
    # print(f)
    return dq, fq, s2

#####################################################
# ##### program main window #######


window = tk.Tk()
window.title("ARPA-CS User Interface")
# window.geometry("1000x700")

# frame
frame = tk.Frame(window)
frame.pack()

# info frame
info_frame = tk.LabelFrame(frame, text="Excitation Signal Control & System Management", font=("Helvetica", 12, "bold"))
info_frame.grid(row=0, column=0, sticky="news", padx=10, pady=5)

length_label = tk.Label(info_frame, text="Length (2^K)")
length_label.grid(row=0, column=0)
length_entry = ttk.Entry(info_frame)
length_entry.grid(row=0, column=1)
seq_label = tk.Label(info_frame, text="Sequency")
seq_label.grid(row=0, column=2)
seq_entry = tk.Entry(info_frame)
seq_entry.grid(row=0, column=3)

# dashboard_frame
dashboard_frame = tk.LabelFrame(frame, text="Observed Signal Analysis", font=("Helvetica", 12, "bold"))
dashboard_frame.grid(row=1, column=0, sticky="news", padx=10, pady=5)

#da1 = [0, 0, 0]
da1 = [0] * 1024
da2 = [0] * 1024
da3 = [0] * 1024
da4 = [0] * 1024
fig, axs = plt.subplots(2, 2, figsize=(6, 3))
axs[0, 0].plot(da1)
axs[0, 0].grid(False)
axs[0, 0].set_title("Received Signal at 590 nm Photodetector")
axs[0, 0].set_xlabel("Sample Index")
axs[0, 0].set_ylabel("Fluorescence Intensity (a.u.)")
axs[1, 0].plot(da2)
axs[1, 0].grid(False)
axs[1, 0].set_title("Walsh-ordered Walsh-Hadamard Transform")
axs[1, 0].set_xlabel("Sequency Index")
axs[1, 0].set_ylabel("Magnitude (a.u)")
axs[0, 1].plot(da3)
axs[0, 1].grid(False)
axs[0, 1].set_title("Received Signal at 630 nm Photodetector")
axs[0, 1].set_xlabel("Sample Index")
axs[0, 1].set_ylabel("Fluorescence Intensity (a.u.)")
axs[1, 1].plot(da4)
axs[1, 1].grid(False)
axs[1, 1].set_title("Walsh-ordered Walsh-Hadamard Transform")
axs[1, 1].set_xlabel("Sequency Index")
axs[1, 1].set_ylabel("Magnitude (a.u.)")
plt.tight_layout()
# plt.show(block=False)
canvas = FigureCanvasTkAgg(fig, dashboard_frame)
canvas.draw()
canvas.get_tk_widget().pack(expand=True, fill="both", padx=5, pady=5)


# summary frame
summary_frame = tk.LabelFrame(frame, text="Result Summary", font=("Helvetica", 12, "bold"))
summary_frame.grid(row=2, column=0, sticky="news", padx=10, pady=5)
#####################################################
summary_label1 = tk.Label(summary_frame, text=" ", fg='#00008B')
summary_label1.grid(row=1, column=1, sticky="w")

summary_label2 = tk.Label(summary_frame, text=" ", fg='#00008B')
summary_label2.grid(row=1, column=2, sticky="w")

summary_label3 = tk.Label(summary_frame, text=" ", fg='#00008B')
summary_label3.grid(row=1, column=3, sticky="w")

summary_label4 = tk.Label(summary_frame, text=" ", fg='#f00')
summary_label4.grid(row=1, column=4, sticky="w")

summary_label5 = tk.Label(summary_frame, text=" ", fg='#00008B')
summary_label5.grid(row=2, column=1, sticky="w")

summary_label6 = tk.Label(summary_frame, text=" ", fg='#00008B')
summary_label6.grid(row=2, column=2, sticky="w")

summary_label7 = tk.Label(summary_frame, text=" ", fg='#00008B')
summary_label7.grid(row=2, column=3, sticky="w")

summary_label8 = tk.Label(summary_frame, text=" ", fg='#f00')
summary_label8.grid(row=2, column=4, sticky="w")


def press():

    #  check if entered values are correct
    if seq_entry.get() != '' and length_entry.get() != '':

        #  check if sequency value is greater than the length of reference signal
        if int(seq_entry.get()) >= 2**int(length_entry.get()):
            summary_label1.config(text="Error: Sequency value is greater than or equal signal length.")
        else:
            summary_label1.config(text="Status: In Progress...\n ")
            w1 = walsh_generator(int(seq_entry.get()), int(length_entry.get()))
            w2 = walsh_generator(int(seq_entry.get()), int(length_entry.get()))
            # #########
            seq1 = copy.deepcopy(int(seq_entry.get()))
            seq2 = copy.deepcopy(int(seq_entry.get()))
            d1, f1, s1 = myprogram1(w1, 590, seq1)
            # print(f1)
            d2, f2, s2 = myprogram2(w2, 630, seq2)
            # print(f2)
            fig.clear()
            fig.add_subplot(221).plot(d1, linewidth=1)
            plt.title("Received signal at 590 nm photodetector")
            plt.xlabel("Sample Index")
            plt.ylabel("Fluorescence Intensity (a.u.)")

            fig.add_subplot(223).plot(f1, '-o',   markevery=[seq1], c='yellow', mfc='red', mec='k',
                                      markersize=4, linewidth=1)
            #plt.title("Walsh-Hadamard Transform (WHT)")
            plt.title("Walsh-ordered Walsh-Hadamard Transform")
            plt.xlabel("Sequency Index")
            plt.ylabel("Magnitude (a.u.)")
            fig.add_subplot(222).plot(d2, linewidth=1)
            plt.title("Received Signal at 630 nm Photodetector")
            plt.xlabel("Sample index")
            plt.ylabel("Fluorescence Intensity (a.u.)")
            fig.add_subplot(224).plot(f2, '-o',   markevery=[seq1], c='yellow', mfc='red', mec='k',
                                      markersize=4, linewidth=1)
            plt.title("Walsh-ordered Walsh-Hadamard Transform")
            plt.xlabel("Sequency Index")
            plt.ylabel("Magnitude (a.u.)")
            canvas.draw_idle()
            dp1 = f1[int(seq_entry.get())]
            dp2 = f2[int(seq_entry.get())]
            r12 = dp2/dp1
            summary_label1.config(text="Photodetector (590 nm):         ", font=("Helvetica", 12, "bold"))
            summary_label2.config(text="Fluorescence intensity = {:05.3f} (a.u)         ".format(dp1), font=("Helvetica", 12, "bold"))
            summary_label3.config(text="Cyclic shift count = {}     ".format(s1), font=("Helvetica", 12, "bold"))
            summary_label4.config(text="Ratio (630/590):", font=("Helvetica", 12, "bold"))
            summary_label5.config(text="Photodetector (630 nm):         ", font=("Helvetica", 12, "bold"))
            summary_label6.config(text="Fluorescence intensity = {:05.3f} (a.u)         ".format(dp2), font=("Helvetica", 12, "bold"))
            summary_label7.config(text="Cyclic shift count = {}     ".format(s2), font=("Helvetica", 12, "bold"))
            summary_label8.config(text="{:.3f}".format(r12), font=("Helvetica", 12, "bold"))

    else:
        summary_label1.config(text="Error: Please enter missing value(s)")
        fig.clear()
        c = [0, 0, 0]
        fig.add_subplot(221).plot(c, color='green', marker='o', linestyle='dashed',
                                  linewidth=0.1, markersize=1, label='line 2')
        plt.grid(color='w', linestyle='-', linewidth=0.2)
        plt.title("Received Signal at 590 nm Photodetector")
        plt.xlabel("Sample Index")
        plt.ylabel("Fluorescence Intensity (a.u.)")

        fig.add_subplot(223).plot(c, linewidth=0.1)
        plt.title("Walsh-ordered Walsh-Hadamard Transform")
        plt.xlabel("Sequency Index")
        plt.ylabel("Magnitude (a.u.)")

        fig.add_subplot(222).plot(c)
        plt.title("Received signal at 630 nm photodetector")
        plt.xlabel("Sample Index")
        plt.ylabel("Fluorescence Intensity (a.u.)")

        fig.add_subplot(224).plot(c)
        plt.title("Walsh-ordered Walsh-Hadamard Transform")
        plt.xlabel("Sequency Index")
        plt.ylabel("Magnitude (a.u.)")
        canvas.draw_idle()

# Buttons
#####################################################


button = Button(info_frame, text="Start Data Acquisition", fg='#fff', bg='#548B54', command=press)
button.grid(row=0, column=4)

button2 = Button(info_frame, text="   Start Calibration   ", fg='#fff', bg='#FF6A6A', command=press)
button2.grid(row=0, column=5)

button3 = Button(info_frame, text="   Check Connection   ", fg='#fff', bg='#7B68EE', command=press)
button3.grid(row=0, column=6)

for widget in info_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)
#####################################################

window.mainloop()
