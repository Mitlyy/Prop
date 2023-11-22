import tkinter as tk
from tkinter import ttk
from tkinter import * 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from propp import *
plt.style.use(['science', 'notebook', 'grid'])

# this is the function called when the button is clicked
def start_calc():
    # w_0, c, l0 = initial(float(lenght.get()))
    ZnTe = material(n2 = float(n_cube.get())*1e-12, n1 = float(n_quad.get())*1e-12, N0 = 3)
    f = Pulse(tau = float(tau.get()), n2 = ZnTe.n2, n1 = ZnTe.n1, N0 = ZnTe.N0, I = float(intensivity.get()), l = float(lenght.get()))
    if Crank_var.get():
        f.crank_nik(steps)
    if split_var.get():
        f.split_step(steps)
        
    plot1.set_xlim([0, 10])
    try:
        plot1.plot(f.time[len(f.time)-len(f.E_out):]/abs(f.norm()[0]), abs(f.E_out))
    except: 
        plot1.plot(f.time[len(f.time)-len(f.E_out):], abs(f.E_out))
    
    # plot1.plot(f.time[len(f.time)-len(f.E_out):]/abs(f.time[len(f.time)-len(f.E_out):][abs(f.E0_k) == max(abs(f.E0_k))][0])
    #            , abs(f.E_out)/max(abs(f.E_out))[0])
    # plot1.set_xlim(0, 1e-14)
    # plot1.plot(f.time[len(f.time)-len(f.E_out):]/16.53, abs(f.E_out))
    fig.canvas.draw()

def clear():
    plot1.cla()
    plot1.set_xlabel(r'$\nu/\nu_{max}$')
    plot1.set_ylabel(r'$|G|$')
    fig.canvas.draw()
    
def E_in():
    ZnTe = material(n2 = float(n_cube.get())*1e-12, n1 = float(n_quad.get())*1e-12, N0 = 2.48)
    f = Pulse(tau = float(tau.get()), n2 = ZnTe.n2, n1 = ZnTe.n1, N0 = ZnTe.N0, I = float(intensivity.get()), l = float(lenght.get()))
    plt.plot(f.time/f.tau, f.gauss(f.time))
    # print(f.w_0)
    plt.xlabel(r'$t/\tau$')
    plt.ylabel(r'$E/E_0$')
    # plt.plot(f.time/f.tau, f.B1)
    # plt.plot(f.time/f.tau, f.B2)
    
    plt.show()
    

# this is a function to check the status of the checkbox (1 means checked, and 0 means unchecked)
def getCheckboxValue():
    checkedOrNot = Crank_var.get()
    return checkedOrNot


# this is a function to check the status of the checkbox (1 means checked, and 0 means unchecked)
def getCheckboxValue():
	checkedOrNot = split_var.get()
	return checkedOrNot


# this is a function to get the user input from the text input box
def getInputBoxValue():
	userInput = n_quad.get()
	return userInput


# this is a function to get the user input from the text input box
def getInputBoxValue():
	userInput = n_cube.get()
	return userInput


# this is a function to get the user input from the text input box
def getInputBoxValue():
	userInput = intensivity.get()
	return userInput

def getInputBoxValue():
	userInput = lenght.get()
	return userInput


root = Tk()
#this is the declaration of the variable associated with the checkbox
Crank_var = tk.IntVar()
#this is the declaration of the variable associated with the checkbox
split_var = tk.IntVar()


w, h = 1280, 640

# This is the section of code which creates the main window
root.geometry('1280x640')
# root.configure(background='#426A8C')
root.title('Prop')

# This is the section of code which creates a button
Button(root, text='Calc', font=('arial', 12, 'normal'), command=start_calc).place(x=777, y = h - 70 - 8)
# root["bg"] = "white"
Button(root, text='Clear', font=('arial', 12, 'normal'), command=clear).place(x=830, y = h - 70 - 8)
Button(root, text='E_in', font=('arial', 12, 'normal'), command=E_in).place(x=900, y = h - 70 - 8)
# This is the section of code which creates a checkbox
Crank_Nicolson=Checkbutton(root, text='Crank-Nicolson', bg = "white", variable=Crank_var, font=('arial', 12, 'normal'))
Crank_Nicolson.place(x=17, y = h - 70 - 5 )



# This is the section of code which creates a checkbox
Split_step=Checkbutton(root, text='Split-step', bg = "white", variable=split_var, font=('arial', 12, 'normal'))
Split_step.place(x=17, y= h - 70 - 25 - 8)


# This is the section of code which creates the a label
Label(root, text='chi_2(1e-12):', font=('arial', 12, 'normal'),  bg = "white").place(x=170, y =  h - 70 - 25-2)


# This is the section of code which creates the a label
Label(root, text='n_2(1e-12):', font=('arial', 12, 'normal'),  bg = "white").place(x=170, y = h - 70-2)
root["bg"] = "white"

# This is the section of code which creates a text input box
n_quad=Entry(root)
n_quad.place(x=267, y=h - 70 - 25)



# This is the section of code which creates a text input box
n_cube=Entry(root)
n_cube.place(x=267, y=h - 70)

# This is the section of code which creates a text input box
intensivity=Entry(root)
intensivity.place(x=470, y = h - 70)

# This is the section of code which creates the a label
Label(root, text='I:', font=('arial', 12, 'normal'), bg = "white").place(x=417, y=h - 70 -2)
fig = Figure(figsize = (13, 5.4),
                dpi = 100)


lenght=Entry(root)
lenght.place(x=470, y = h - 76 -25)

Label(root, text='lmbd:', font=('arial', 12, 'normal'), bg = "white").place(x=417, y=h - 76 -25)


Label(root, text='tau:', font=('arial', 12, 'normal'), bg = "white").place(x=607, y=h - 76 -25)
tau=Entry(root)
tau.place(x=640, y=h - 74 - 25)
#n_quad = 2.5e-12




# adding the subplot
plot1 = fig.add_subplot(111)
plot1.set_xlabel(r'$\nu/\nu_{max}$')
plot1.set_ylabel(r'$|G|$')
# plotting the graph

# creating the Tkinter canvas
# containing the Matplotlib figure
canvas = FigureCanvasTkAgg(fig,
                            master = root)  
canvas.draw()

# placing the canvas on the Tkinter window
canvas.get_tk_widget().pack()

# creating the Matplotlib toolbar
toolbar = NavigationToolbar2Tk(canvas,
                                root)
toolbar.update()





root.mainloop()
