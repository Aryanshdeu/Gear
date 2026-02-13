import tkinter as tk
from tkinter import ttk, messagebox
import math
import sqlite3
import os



db_path = 'Bevel_Gear_Results.db'
if os.path.exists(db_path):
    os.remove(db_path)


def create_database(db_path='Bevel_Gear_Results.db'):

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS results")

    cursor.execute('''
        CREATE TABLE results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            input_Nump INTEGER,
            input_Ng INTEGER,
            input_Pd REAL,
            input_phi REAL,
            input_np REAL,
            input_P REAL,
            output_Pitch_Diameter_Pinion REAL,
            output_Pitch_Diameter_Gear REAL,
            output_Pitch_Cone_Angle_Pinion REAL,
            output_Pitch_Cone_Angle_Gear REAL,
            output_Angle_Check REAL,
            output_Whole_depth REAL,
            output_Working_depth REAL,
            output_Clearance REAL,
            output_Outside_Diameter_Pinion REAL,
            output_Outside_Diameter_Gear REAL,
            output_Face_Cone_Distance REAL,
            output_Face_Width REAL,
            output_Rotational_Speed_Pinion REAL,
            output_Rotational_Speed_Gear REAL,
            output_Torque_Pinion REAL,
            output_Torque_Gear REAL,
            output_Tangential_Force_Pinion REAL,
            output_Tangential_Force_Gear REAL,
            output_Radial_Force_Pinion REAL,
            output_Radial_Force_Gear REAL,
            output_Axial_Force_Pinion REAL,
            output_Axial_Force_Gear REAL
        )
    ''')
    conn.commit()
    return conn


def save_to_database(conn, inputs, results):
    cursor = conn.cursor()
    try:
        cursor.execute('''
            INSERT INTO results(
                input_Nump, input_Ng, input_Pd, input_phi, input_np, input_P,
                output_Pitch_Diameter_Pinion, output_Pitch_Diameter_Gear,
                output_Pitch_Cone_Angle_Pinion, output_Pitch_Cone_Angle_Gear,
                output_Angle_Check, output_Whole_depth, output_Working_depth,
                output_Clearance, output_Outside_Diameter_Pinion,
                output_Outside_Diameter_Gear, output_Face_Cone_Distance,
                output_Face_Width, output_Rotational_Speed_Pinion,
                output_Rotational_Speed_Gear, output_Torque_Pinion,
                output_Torque_Gear, output_Tangential_Force_Pinion,
                output_Tangential_Force_Gear, output_Radial_Force_Pinion,
                output_Radial_Force_Gear, output_Axial_Force_Pinion,
                output_Axial_Force_Gear
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
            (
                inputs['Num_p'], inputs['N_g'], inputs['p_d'], inputs['phi'], inputs['n_p'], inputs['P'],
                results['Pitch Diameter Pinion (in)'], results['Pitch Diameter Gear (in)'],
                results['Pitch Cone Angle Pinion (deg)'], results['Pitch Cone Angle Gear (deg)'],
                results['Angle Check (deg)'], results['Whole depth (in)'],
                results['Working depth (in)'], results['Clearance (in)'],
                results['Outside Diameter Pinion (in)'], results['Outside Diameter Gear (in)'],
                results['Face Cone Distance (in)'], results['Face Width (in)'],
                results['Rotational Speed Pinion (rpm)'], results['Rotational Speed Gear (rpm)'],
                results['Torque Pinion (lb-in)'], results['Torque Gear (lb-in)'],
                results['Tangential Force Pinion(lb)'], results['Tangential Force Gear(lb)'],
                results['Radial Force Pinion(lb)'], results['Radial Force Gear(lb)'],
                results['Axial Force Pinion(lb)'], results['Axial Force Gear(lb)']
            ))
        conn.commit()
    except Exception as e:
        messagebox.showerror("Database Error", f"Failed to save data: {e}")


def calculate_bevel_gear(Num_p, N_g, p_d, phi, n_p, P):
    try:
        Num_p = int(Num_p)
        N_g = int(N_g)
        p_d = float(p_d)
        phi_deg = float(phi)
        phi_rad = math.radians(phi_deg)
        n_p = float(n_p)
        P = float(P)

        d = Num_p / p_d
        D = N_g / p_d

        gamma = math.degrees(math.atan(Num_p / N_g))
        Gamma = math.degrees(math.atan(N_g / Num_p))
        angle_check = gamma + Gamma

        h_t = (2.188 / p_d) + 0.002
        h_k = 2.00 / p_d
        c = (0.188 / p_d) + 0.002

        a_g = (0.54 / p_d) + (0.46 / p_d * (N_g / Num_p) ** 2)
        a_p = h_k - a_g

        D_o = D + 2 * a_g * math.cos(math.radians(Gamma))
        d_o = d + 2 * a_p * math.cos(math.radians(gamma))

        A_o = math.sqrt((d / 2) ** 2 + (D / 2) ** 2)
        F = A_o / 3

        r_m = (d / 2) - (F / 2) * math.sin(phi_rad) * math.cos(math.radians(gamma))
        R_m = (D / 2) - (F / 2) * math.sin(phi_rad) * math.cos(math.radians(Gamma))

        n_g = n_p * (Num_p / N_g)

        T_p = 63000 * (P / Num_p)
        T_g = 63000 * (P / n_g)

        w_t_p = T_p / r_m
        w_t_g = T_g / R_m
        w_r_p = w_t_p * math.tan(phi_rad) * math.cos(math.radians(gamma))
        w_r_g = w_t_g * math.tan(phi_rad) * math.cos(math.radians(Gamma))
        w_x_p = w_t_p * math.tan(phi_rad) * math.sin(math.radians(gamma))
        w_x_g = w_t_g * math.tan(phi_rad) * math.sin(math.radians(Gamma))

        results = {
            'Pitch Diameter Pinion (in)': d,
            'Pitch Diameter Gear (in)': D,
            'Pitch Cone Angle Pinion (deg)': gamma,
            'Pitch Cone Angle Gear (deg)': Gamma,
            'Angle Check (deg)': angle_check,
            'Whole depth (in)': h_t,
            'Working depth (in)': h_k,
            'Clearance (in)': c,
            'Outside Diameter Pinion (in)': d_o,
            'Outside Diameter Gear (in)': D_o,
            'Face Cone Distance (in)': A_o,
            'Face Width (in)': F,
            'Rotational Speed Pinion (rpm)': n_p,
            'Rotational Speed Gear (rpm)': n_g,
            'Torque Pinion (lb-in)': T_p,
            'Torque Gear (lb-in)': T_g,
            'Tangential Force Pinion(lb)': w_t_p,
            'Tangential Force Gear(lb)': w_t_g,
            'Radial Force Pinion(lb)': w_r_p,
            'Radial Force Gear(lb)': w_r_g,
            'Axial Force Pinion(lb)': w_x_p,
            'Axial Force Gear(lb)': w_x_g
        }

        inputs = {'Num_p': Num_p, 'N_g': N_g, 'p_d': p_d, 'phi': phi_deg, 'n_p': n_p, 'P': P}
        return inputs, results
    except (ValueError, ZeroDivisionError) as e:
        return None, {'Error': f'Invalid input: {e}'}


def run_calculation():
    N_p_val = entry_Nump.get()
    N_g_val = entry_Ng.get()
    p_d_val = entry_pd.get()
    phi_val = entry_phi.get()
    n_p_val = entry_np.get()
    P_val = entry_P.get()

    inputs, results = calculate_bevel_gear(N_p_val, N_g_val, p_d_val, phi_val, n_p_val, P_val)

    if 'Error' in results:
        messagebox.showerror('Input Error', results['Error'])
        return

    conn = create_database()
    save_to_database(conn, inputs, results)
    conn.close()

    show_results(results)


def show_results(results):
    new_window = tk.Toplevel(root)
    new_window.title('Calculation Results')
    new_window.geometry('400x550')
    new_window.config(bg="#f0f0f0")
    result_frame = ttk.Frame(new_window, padding='10')
    result_frame.pack(pady=10, fill='x')

    for key, value in results.items():
        label_text = f'{key}: {value:.2f}'
        ttk.Label(result_frame, text=label_text).pack(anchor='w', pady=2)

    ttk.Button(new_window, text="Close", command=new_window.destroy).pack(pady=10)


def Delete():
    entry_Nump.delete(0, tk.END)
    entry_Ng.delete(0, tk.END)
    entry_pd.delete(0, tk.END)
    entry_phi.delete(0, tk.END)
    entry_np.delete(0, tk.END)
    entry_P.delete(0, tk.END)


# GUI
root = tk.Tk()
root.title('Bevel Gear Design Calculator')
root.geometry('450x650')
root.config(bg='#f0f0f0')

input_frame = ttk.Frame(root, padding='10')
input_frame.pack(pady=10)

ttk.Label(input_frame, text='Pinion Teeth (Np):').grid(row=0, column=0, sticky=tk.W, pady=5)
entry_Nump = ttk.Entry(input_frame)
entry_Nump.grid(row=0, column=1)

ttk.Label(input_frame, text='Gear Teeth (Ng):').grid(row=1, column=0, sticky=tk.W, pady=5)
entry_Ng = ttk.Entry(input_frame)
entry_Ng.grid(row=1, column=1)

ttk.Label(input_frame, text='Diametral Pitch (Pd):').grid(row=2, column=0, sticky=tk.W, pady=5)
entry_pd = ttk.Entry(input_frame)
entry_pd.grid(row=2, column=1)

ttk.Label(input_frame, text='Pressure Angle(phi,deg):').grid(row=3, column=0, sticky=tk.W, pady=5)
entry_phi = ttk.Entry(input_frame)
entry_phi.grid(row=3, column=1)

ttk.Label(input_frame, text='Pinion Speed (np,rpm):').grid(row=4, column=0, sticky=tk.W, pady=5)
entry_np = ttk.Entry(input_frame)
entry_np.grid(row=4, column=1)

ttk.Label(input_frame, text='Power (P,hp):').grid(row=5, column=0, sticky=tk.W, pady=5)
entry_P = ttk.Entry(input_frame)
entry_P.grid(row=5, column=1)

ttk.Button(root, text='Calculate', command=run_calculation).pack(pady=10)
ttk.Button(root, text='Delete', command=Delete).pack(pady=10)

root.mainloop()
