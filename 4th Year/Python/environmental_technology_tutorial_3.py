import math
import matplotlib.pyplot as plt

# Constants
theta_k1_less_20 = 1.135  # Temperature correction factor for k1 when T < 20°C
theta_k1_more_20 = 1.047  # Temperature correction factor for k1 when T > 20°C
U_to_mps = 1 / 60  # Conversion factor from m/min to m/s
seconds_per_day = 86400  # Number of seconds in a day

def calculate_DO_saturation(T):
    # Calculation of DO at saturation
    DO_sat = 1 / (0.0686 + 0.00188432 * T + 0.00000612 * T ** 2)
    return DO_sat

def calculate_k1(k1_20, T):
    # Calculate k1 at temperature T
    if T < 20:
        theta = theta_k1_less_20
    else:
        theta = theta_k1_more_20
    k1 = k1_20 * theta ** (T - 20)
    return k1

def calculate_k2(U_mps, H_mm, T):
    # Calculate k2 at temperature T
    f_20 = 78200 * U_mps ** 0.67 * H_mm ** -0.85
    f_T = f_20 * 0.98 ** (20 - T)  # Temperature correction
    k2 = f_T / H_mm  # h⁻¹
    k2 = k2 * 24  # Convert to d⁻¹
    return k2

def calculate_BODu(BOD5, k1):
    # Calculate ultimate BOD
    BODu = BOD5 / (1 - math.exp(-k1 * 5))
    return BODu

# Part (a)
def part_a():
    print("Part (a): Calculating DOmin and distance downstream")
    # Given data
    Qr = 15  # m³/min
    Qw = 0.25 * Qr
    Qc = Qr + Qw
    Lr = 3  # mg/L
    Lw = 25  # mg/L
    DO_r = 9  # mg/L
    DO_w = 2  # mg/L
    U = 1.0  # m/min
    H = 1.5  # m
    T = 17  # °C
    k1_20 = 0.230  # d⁻¹

    # Step 1: Combined BOD5 and DO
    Lc = (Qr * Lr + Qw * Lw) / Qc
    DO_c = (Qr * DO_r + Qw * DO_w) / Qc

    # Step 2: Calculate k1 at T°C
    k1 = calculate_k1(k1_20, T)

    # Step 3: Calculate ultimate BOD (BODu_c)
    BODu_c = calculate_BODu(Lc, k1)

    # Step 4: Calculate k2 at T°C
    U_mps = U * U_to_mps
    H_mm = H * 1000
    k2 = calculate_k2(U_mps, H_mm, T)

    # Step 5: Calculate DO saturation
    DO_sat = calculate_DO_saturation(T)

    # Step 6: Initial oxygen deficit D0
    D0 = DO_sat - DO_c

    # Step 7: Calculate tmax
    k_diff = abs(k2 - k1)
    numerator = k2 / k1
    denominator = 1 - ((k2 - k1) * D0) / (k1 * BODu_c)
    if denominator <= 0:
        print("Error: Denominator in tmax calculation is non-positive.")
        return
    tmax = (1 / k_diff) * math.log(numerator * denominator)

    # Step 8: Calculate Dmax
    Dmax = (k1 * BODu_c / k2) * math.exp(-k1 * tmax)

    # Step 9: Calculate DOmin
    DOmin = DO_sat - Dmax

    # Step 10: Distance downstream at DOmin
    distance = U_mps * tmax * seconds_per_day / 1000  # km

    # Print results
    print(f"Combined BOD5 (Lc): {Lc:.2f} mg/L")
    print(f"Combined DO (DO_c): {DO_c:.2f} mg/L")
    print(f"k1 at {T}°C: {k1:.4f} d⁻¹")
    print(f"Combined ultimate BOD (BODu_c): {BODu_c:.2f} mg/L")
    print(f"k2 at {T}°C: {k2:.4f} d⁻¹")
    print(f"DO saturation at {T}°C (DO*): {DO_sat:.2f} mg/L")
    print(f"Initial oxygen deficit (D0): {D0:.2f} mg/L")
    print(f"Time to reach maximum deficit (tmax): {tmax:.2f} days")
    print(f"Maximum oxygen deficit (Dmax): {Dmax:.2f} mg/L")
    print(f"Minimum dissolved oxygen concentration (DOmin): {DOmin:.2f} mg/L")
    print(f"Distance downstream at DOmin: {distance:.2f} km")

    return {
        'BODu_c': BODu_c,
        'tmax': tmax,
        'Dmax': Dmax,
        'DOmin': DOmin,
        'distance': distance,
        'k1': k1,
        'k2': k2,
        'D0': D0,
        'DO_sat': DO_sat,
        'U_mps': U_mps
    }

# Execute Part (a)
results_a = part_a()

# Part (b)
def part_b(results):
    print("\nPart (b): Plotting the oxygen sag curve")
    distances_km = [i for i in range(0, 26, 2)]
    times = [distance * 1000 / (results['U_mps'] * seconds_per_day) for distance in distances_km]
    D_t = []
    for t in times:
        D = ((results['k1'] * results['BODu_c']) / (results['k2'] - results['k1'])) * \
            (math.exp(-results['k1'] * t) - math.exp(-results['k2'] * t)) + \
            results['D0'] * math.exp(-results['k2'] * t)
        D_t.append(D)
        print(f"At distance {distances_km[times.index(t)]} km (time {t:.2f} days): Deficit D(t) = {D:.2f} mg/L")

    # Plotting the oxygen sag curve
    plt.figure()
    plt.plot(distances_km, D_t, marker='o')
    plt.title('Oxygen Sag Curve')
    plt.xlabel('Distance downstream (km)')
    plt.ylabel('Oxygen Deficit D(t) (mg/L)')
    plt.grid(True)
    plt.show()

# Execute Part (b)
part_b(results_a)

# Part (c)
def part_c():
    print("\nPart (c): Effects of drought conditions")
    # Updated velocity and depth
    U = 0.5  # m/min
    H = 1.0  # m

    # Flow rates proportional to velocity and depth
    Qr = 15 * (U / 1.0) * (H / 1.5)
    Qw = 0.25 * Qr
    Qc = Qr + Qw

    # Same BOD5 and DO concentrations
    Lr = 3
    Lw = 25
    DO_r = 9
    DO_w = 2

    # Temperature remains the same
    T = 17
    k1_20 = 0.230

    # Recalculate combined concentrations
    Lc = (Qr * Lr + Qw * Lw) / Qc
    DO_c = (Qr * DO_r + Qw * DO_w) / Qc

    # Recalculate k1 and BODu_c
    k1 = calculate_k1(k1_20, T)
    BODu_c = calculate_BODu(Lc, k1)

    # Recalculate k2
    U_mps = U * U_to_mps
    H_mm = H * 1000
    k2 = calculate_k2(U_mps, H_mm, T)

    # Recalculate DO saturation and D0
    DO_sat = calculate_DO_saturation(T)
    D0 = DO_sat - DO_c

    # Recalculate tmax and Dmax
    k_diff = abs(k2 - k1)
    numerator = k2 / k1
    denominator = 1 - ((k2 - k1) * D0) / (k1 * BODu_c)
    if denominator <= 0:
        print("Error: Denominator in tmax calculation is non-positive.")
        return
    tmax = (1 / k_diff) * math.log(numerator * denominator)
    Dmax = (k1 * BODu_c / k2) * math.exp(-k1 * tmax)
    DOmin = DO_sat - Dmax

    # Distance downstream at DOmin
    distance = U_mps * tmax * seconds_per_day / 1000  # km

    # Print results
    print(f"Combined BOD5 (Lc): {Lc:.2f} mg/L")
    print(f"Combined DO (DO_c): {DO_c:.2f} mg/L")
    print(f"k1 at {T}°C: {k1:.4f} d⁻¹")
    print(f"Combined ultimate BOD (BODu_c): {BODu_c:.2f} mg/L")
    print(f"k2 at {T}°C: {k2:.4f} d⁻¹")
    print(f"DO saturation at {T}°C (DO*): {DO_sat:.2f} mg/L")
    print(f"Initial oxygen deficit (D0): {D0:.2f} mg/L")
    print(f"Time to reach maximum deficit (tmax): {tmax:.2f} days")
    print(f"Maximum oxygen deficit (Dmax): {Dmax:.2f} mg/L")
    print(f"Minimum dissolved oxygen concentration (DOmin): {DOmin:.2f} mg/L")
    print(f"Distance downstream at DOmin: {distance:.2f} km")

    # Plot new oxygen sag curve
    distances_km = [i for i in range(0, 26, 2)]
    times = [distance * 1000 / (U_mps * seconds_per_day) for distance in distances_km]
    D_t = []
    for t in times:
        D = ((k1 * BODu_c) / (k2 - k1)) * \
            (math.exp(-k1 * t) - math.exp(-k2 * t)) + \
            D0 * math.exp(-k2 * t)
        D_t.append(D)
        print(f"At distance {distances_km[times.index(t)]} km (time {t:.2f} days): Deficit D(t) = {D:.2f} mg/L")

    plt.figure()
    plt.plot(distances_km, D_t, marker='o', color='red')
    plt.title('Oxygen Sag Curve under Drought Conditions')
    plt.xlabel('Distance downstream (km)')
    plt.ylabel('Oxygen Deficit D(t) (mg/L)')
    plt.grid(True)
    plt.show()

    return {
        'BODu_c': BODu_c,
        'tmax': tmax,
        'Dmax': Dmax,
        'DOmin': DOmin,
        'distance': distance
    }

# Execute Part (c)
results_c = part_c()

# Part (d)
def part_d():
    print("\nPart (d): Effects of increased temperature to 20°C")
    # Use same flow rates and concentrations as in part (c)
    U = 0.5
    H = 1.0
    Qr = 15 * (U / 1.0) * (H / 1.5)
    Qw = 0.25 * Qr
    Qc = Qr + Qw
    Lr = 3
    Lw = 25
    DO_r = 9
    DO_w = 2
    T = 20  # Increased temperature
    k1_20 = 0.230

    Lc = (Qr * Lr + Qw * Lw) / Qc
    DO_c = (Qr * DO_r + Qw * DO_w) / Qc

    # Recalculate k1
    k1 = calculate_k1(k1_20, T)
    BODu_c = calculate_BODu(Lc, k1)

    # Recalculate k2
    U_mps = U * U_to_mps
    H_mm = H * 1000
    # At 20°C, f_T = f_20 (no temperature correction)
    f_20 = 78200 * U_mps ** 0.67 * H_mm ** -0.85
    k2 = (f_20 / H_mm) * 24  # Convert to d⁻¹

    # Recalculate DO saturation and D0
    DO_sat = calculate_DO_saturation(T)
    D0 = DO_sat - DO_c

    # Recalculate tmax and Dmax
    k_diff = abs(k2 - k1)
    numerator = k2 / k1
    denominator = 1 - ((k2 - k1) * D0) / (k1 * BODu_c)
    if denominator <= 0:
        print("Error: Denominator in tmax calculation is non-positive.")
        return
    tmax = (1 / k_diff) * math.log(numerator * denominator)
    Dmax = (k1 * BODu_c / k2) * math.exp(-k1 * tmax)
    DOmin = DO_sat - Dmax

    # Distance downstream at DOmin
    distance = U_mps * tmax * seconds_per_day / 1000  # km

    # Print results
    print(f"Combined BOD5 (Lc): {Lc:.2f} mg/L")
    print(f"Combined DO (DO_c): {DO_c:.2f} mg/L")
    print(f"k1 at {T}°C: {k1:.4f} d⁻¹")
    print(f"Combined ultimate BOD (BODu_c): {BODu_c:.2f} mg/L")
    print(f"k2 at {T}°C: {k2:.4f} d⁻¹")
    print(f"DO saturation at {T}°C (DO*): {DO_sat:.2f} mg/L")
    print(f"Initial oxygen deficit (D0): {D0:.2f} mg/L")
    print(f"Time to reach maximum deficit (tmax): {tmax:.2f} days")
    print(f"Maximum oxygen deficit (Dmax): {Dmax:.2f} mg/L")
    print(f"Minimum dissolved oxygen concentration (DOmin): {DOmin:.2f} mg/L")
    print(f"Distance downstream at DOmin: {distance:.2f} km")

    return {
        'BODu_c': BODu_c,
        'tmax': tmax,
        'Dmax': Dmax,
        'DOmin': DOmin,
        'distance': distance
    }

# Execute Part (d)
results_d = part_d()

# Part (e)
def part_e():
    print("\nPart (e): Tabulating results")
    print("{:<15} {:<10} {:<10} {:<10} {:<15}".format('Condition', 'BODu_c', 'tmax (d)', 'Dmax', 'DOmin', 'Distance (km)'))
    print("{:<15} {:<10.2f} {:<10.2f} {:<10.2f} {:<15.2f}".format(
        'Part (a)', results_a['BODu_c'], results_a['tmax'], results_a['Dmax'], results_a['DOmin'], results_a['distance']))
    print("{:<15} {:<10.2f} {:<10.2f} {:<10.2f} {:<15.2f}".format(
        'Part (c)', results_c['BODu_c'], results_c['tmax'], results_c['Dmax'], results_c['DOmin'], results_c['distance']))
    print("{:<15} {:<10.2f} {:<10.2f} {:<10.2f} {:<15.2f}".format(
        'Part (d)', results_d['BODu_c'], results_d['tmax'], results_d['Dmax'], results_d['DOmin'], results_d['distance']))

# Execute Part (e)
part_e()
