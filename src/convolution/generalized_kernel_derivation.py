#%%
import sympy as sp
from sympy import solve
from sympy import *
import itertools
from itertools import combinations
import numpy as np
import math
import time
import threading
import matplotlib.pyplot as plt

def kernel_func(s, coefs_final, coefficients):
    # evaluate the derived kernel
    s_abs = abs(s)
    kernel_val = 0.0
    for i in range(M_eqs):
        if s_abs < i+1:
            for j in range(p_deg+1):
                try:
                    kernel_val += float(coefs_final[coefficients[i, j]]) * s_abs**j
                except:
                    kernel_val += 0.0
            return kernel_val
        else:
            continue
    return kernel_val

def ideal(x):
    # sinc function frequency response (simplified)
    y = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] <= 0.5 and x[i] >= -0.5:
            y[i] = 1.0
        else:
            y[i] = 0.0
    return y

def find_closest(numbers, target):
    # find the index of the number closest to target
    # this function is used for calculating slope of frequency response at cuttoff
    closest_index = 0
    smallest_difference = abs(numbers[0] - target)
    
    for i, num in enumerate(numbers[1:], 1):
        difference = abs(num - target)
        if difference < smallest_difference:
            closest_index = i
            smallest_difference = difference
    
    return closest_index

def generate_symbolic_bvp_equations(M_eqs, p_deg, d_cont):
    """
    Generate symbolic equations for the kernel boundary value problem.
    
    Parameters:
    M_eqs (int): Number of polynomials (equations)
    p_deg (int): Degree of polynomials
    d_cont (int): Highest derivative continuity
    
    Returns:
    tuple: (equations, coefficients)
    """
    
    # Define symbolic coefficients
    A = sp.Matrix([[sp.Symbol(f'A_{{{i+1},{j}}}') for j in range(p_deg + 1)] for i in range(M_eqs)])
    
    equations = []
    
    def polynomial_term(i, j, x, n):
        """Generate term of polynomial and its derivatives"""
        # if j - n < 0:
            # return 0
        coeff = sp.prod([(j-k) for k in range(n)])
        return A[i, j] * coeff * x**(j-n)
    
    # Zeroth derivative (values) - write equation for each end of each polynomial
    for i in range(M_eqs):
        # Left end of the polynomial
        x_left = i
        if i == 0:
            # First boundary condition: u(0) = 1
            equations.append(sum(polynomial_term(i, j, x_left, 0) for j in range(p_deg + 1)) - 1)
            # # extra boundary condition: u(-1) = 0
            # x_left = i - 1
            # equations.append(sum(polynomial_term(i, j, x_left, 0) for j in range(p_deg + 1)) + 1)
        else:
            # Interior point: u(i) = 0
            equations.append(sum(polynomial_term(i, j, x_left, 0) for j in range(p_deg + 1)))
            
        # Right end of the polynomial
        # if i == M_eqs-1: # extra addition (from paper)
        #     pass
        # else:
        x_right = i + 1
        equations.append(sum(polynomial_term(i, j, x_right, 0) for j in range(p_deg + 1)))
    
    # Higher derivatives - loop over nodes
    for n in range(1, d_cont + 1):
        for k in range(M_eqs + 1):
            x = k
            if k == 0:
                # if n == 1:
                if n % 2 == 0:
                    pass
                    # new boundary condition: du(-1) = -du(+1)
                    # equations.append(sum(polynomial_term(0, j, -1, n) for j in range(n, p_deg + 1)) + sum(polynomial_term(0, j, 1, n) for j in range(n, p_deg + 1))) # if p_deg-j-n >= 0))
                else:
                    # Left boundary condition: n-th derivative = 0
                    equations.append(sum(polynomial_term(0, j, x, n) for j in range(n, p_deg + 1))) # if p_deg-j-n >= 0))
            elif k == M_eqs:
                # Right boundary condition: n-th derivative = 0
                equations.append(sum(polynomial_term(M_eqs-1, j, x, n) for j in range(n, p_deg + 1))) # if p_deg-j-n >= 0))
            else:
                # Continuity of n-th derivative between equations
                left_eq = sum(polynomial_term(k-1, j, x, n) for j in range(n,p_deg + 1)) # if p_deg-j-n >= 0)
                right_eq = sum(polynomial_term(k, j, x, n) for j in range(n,p_deg + 1)) # if p_deg-j-n >= 0)
                equations.append(left_eq - right_eq)

    return sp.Matrix(equations), A


class TimeoutException(Exception):
    pass

def solve_with_timeout(A, b, timeout=15):
    result = []
    exception = []

    def solve_function():
        try:
            coefs = sp.solve(A, b)
            result.append(coefs)
        except Exception as e:
            exception.append(e)

    thread = threading.Thread(target=solve_function)
    thread.start()
    thread.join(timeout)

    if thread.is_alive():
        print("Linear solve timed out")
        return None
    elif exception:
        print(f"An error occurred: {exception[0]}")
        return None
    else:
        return result[0]
    
#%%
# Setup the inputs to get the desired kernel
extra_eq = 0 # number of additional equations (improves order of accuracy and frequency response)
M_eqs = 2 + extra_eq # number of polynomials
p_deg = 2*(M_eqs-extra_eq)-1 # degree of polynomials
d_cont = p_deg-2 # highest derivative continuity (can be increased at the expense of order of accuracy and frequency response)
max_taylor = 7 # maximum number of terms in the Taylor expansion (can make derivations computationally demanding, slow)
solver_timeout = 2*60 # seconds till time out of symbolic linear solves
highest_fixed_coef = p_deg # highest degree of fixed coefficient (can be used to limit solution space)
skip_first_equation = False # optional limitation of solution space
print(f"Number of equations: {M_eqs}")
print(f"Degree of polynomials: {p_deg}")
print(f"Highest derivative continuity: {d_cont}")
print(f"Maximum number of terms in the Taylor expansion: {max_taylor}")
print(f"Solver timeout: {solver_timeout} seconds")
print(f"Skip first equation: {skip_first_equation}")
print(f"Highes degree of fixed coefficient: {highest_fixed_coef}")

#%%
# generate the basic equations
equations_raw, coefficients = generate_symbolic_bvp_equations(M_eqs, p_deg, d_cont)

print(f"Number of raw equations: {len(equations_raw)}")
print("\nEquations:")
for i, eq in enumerate(equations_raw):
    print(f"Equation {i+1}: {eq} = 0")

print("\nCoefficients:")
sp.pprint(coefficients[:, ::-1])

#%%
# test different boundary conditions
# until we find some that work

# Calculate the number of coefficients and existing equations
num_coefficients = M_eqs * (p_deg + 1)
num_existing_eqs = len(equations_raw)

# Calculate how many additional equations are needed
num_additional_eqs = num_coefficients - num_existing_eqs
lowest_fixed_coef = math.ceil((M_eqs*p_deg - num_additional_eqs)/M_eqs)
print(f"Number of additional equations needed: {num_additional_eqs}")

# Find coefficients set to specific values
set_coeffs = {coefficients[0,j] for j in range(d_cont+1)}
print(set_coeffs)

#%%
if num_additional_eqs == 0:
    fft_grads = {0: 100} # arbitrary large number
    zeros_vals = {0:'N/A'}
    coefs_all = {0: solve(equations_raw)}
else:
    # Generate additional equations for the lowest power coefficients

    # Generate all possible indices of combinations for the additional equations
    rows, cols = M_eqs, highest_fixed_coef
    all_indices = [(i, j) for i in range(1 if skip_first_equation else 0, rows) for j in range(lowest_fixed_coef,highest_fixed_coef+1)]
    combinations = list(itertools.combinations(all_indices, num_additional_eqs))
    print(f"Number of combinations: {len(combinations)}")
    print(combinations)
    equations = sp.zeros(M_eqs*(p_deg+1), 1)
    for i in range((p_deg+1)*(M_eqs)-num_additional_eqs):
        equations[i] = equations_raw[i]

    # loop through the combinations to find the best solution
    fft_grads = dict()
    coefs_all = dict()
    zeros_vals = dict()
    best_derivative = 0.0
    combi_count = 0
    for combination in combinations:
        print(f"\nCombination {combi_count+1}/{len(combinations)}: {combination}")
        combi_count += 1

        break_flag = False
        for index in combination:
            i, j = index
            if coefficients[i,j] in set_coeffs:
                print(f"Skipping {coefficients[i,j]}, it is already fixed.")
                break_flag = True
                break
        if break_flag:
            continue # skip the rest, go to the next combination
        
        # set the coefficients
        k = 0
        unknowns_range = num_additional_eqs
        a_k = sp.symbols([f'a_{i}' for i in range(unknowns_range)])
        for index in combination:
            i, j = index
            equations[k+M_eqs*(p_deg+1)-num_additional_eqs] = coefficients[i,j] - a_k[k]
            k += 1
            if k == num_additional_eqs:
                break

        # Solve the equations
        sol = sp.solve(equations, coefficients)
        
        ## now set up the kernel equations
        new_eqs = sp.zeros(M_eqs, 1)
        s = sp.symbols('s')
        for i in range(M_eqs):
            for j in range(p_deg+1):
                try:
                    new_eqs[i] += sol[coefficients[i, j]] * s**j
                except:
                    print(f"Empty solution for coefficient[{i},{j}], setting this to zero.")
                    pass

        # evaluate the kernel

        # negative direction
        g_neg = sp.zeros(M_eqs,1)
        for i in range(M_eqs): # go from s and backwards
            g_neg[i] = new_eqs[i].subs(s, s+i).simplify()

        # positive direction
        g_pos = sp.zeros(M_eqs,1)
        for i in range(M_eqs):
            g_pos[i] = sum([(-1)**j * new_eqs[i].coeff(s, j)*s**j for j in range(p_deg+1)]).subs(s, s-i-1).simplify()

        coef_neg = sp.Matrix([sp.Symbol(f'coef_{{{i}}}') for i in range(M_eqs)])
        coef_pos = sp.Matrix([sp.Symbol(f'coef_{{{i}}}') for i in range(M_eqs)])

        g_neg = [coef_neg[i] * g_neg[i] for i in range(M_eqs)]
        g_pos = [coef_pos[i] * g_pos[i] for i in range(M_eqs)]

        # build the taylor series and match coefficients with the kernel

        fp = sp.symbols(f'fp:{max_taylor+1}') # derivatives (index zero is fx)
        h = sp.symbols('h') # step

        taylor_step = sp.zeros(max_taylor, 1)
        print(f"taylor_step: {taylor_step}")
        for i in range(max_taylor):
            taylor_step[i] = fp[i]*h**i/factorial(i)

        for i in range(M_eqs):
            g_neg[i] = g_neg[i].subs(coef_neg[i], sum(taylor_step[0:max_taylor]).subs({h:-i*h}))
            g_pos[i] = g_pos[i].subs(coef_pos[i], sum(taylor_step[0:max_taylor]).subs({h:(i+1)*h}))

        f_taylor = sum(taylor_step[0:max_taylor]).subs({h:s*h})

        g = sum( g_neg[i] + g_pos[i] for i in range(M_eqs) )

        subtract = f_taylor - g
        subtract = sp.collect(sp.expand(subtract),s)

        if subtract.coeff(s, 0) != 0:
            print(f"Skipping combination {combination} due to non-zero constant term.")
            continue # skip the rest, go to the next combination

        # solve the symbolic equations
        eqs_sol = sp.zeros(p_deg, 1)
        unknowns_range = num_additional_eqs if num_additional_eqs < p_deg else p_deg
        for i in range(unknowns_range):
            eqs_sol[i] = subtract.coeff(s, i+1)
        solve_for = sp.symbols([f'a_{i}' for i in range(unknowns_range)])
        print(f"solving equations {eqs_sol}")
        print(f"solving for {solve_for}")

        # Attempt to solve with a timeout
        start_time = time.time()
        coefs = solve_with_timeout(eqs_sol, solve_for, solver_timeout)
        end_time = time.time()

        if coefs is not None:
            print(f"Solved iteration {combi_count} in {end_time - start_time:.2f} seconds")
            if coefs == [] or nan in coefs:
                print(f"Empty solution for {p_deg-j} for combination {combination}. Skipping equations.")
                continue
        else:
            print(f"Moving to next iteration due to timeout or error.")
            continue

        # find the order of accuracy of the derived kernel
        # by setting the highest derivatives to zero one after the other
        # and checking if there are any free symbols in the solution
        set_to_zero = {i:0 for i in a_k}
        break_flag = false
        for min_taylor in range(max_taylor,0,-1):
            taylor_count = 0
            deriv_zeros = {fp[i]:0 for i in range(min_taylor, max_taylor+1)}
            coefs_final = {key:sol[key].subs(coefs).subs(set_to_zero).subs(deriv_zeros) for key in sol}
            for index in coefs_final.keys():
                expr = coefs_final[index]
                symbols = expr.free_symbols
                if not symbols:
                    taylor_count += 1
                    if taylor_count == M_eqs*(p_deg+1):
                        print(f"Eliminated symbols by setting Taylor terms >= {min_taylor} to zero.")
                        break_flag = True # all good, no symbols
                        min_taylor_save = min_taylor
                        break
                    else:
                        pass
                else:
                    print(f'Free symbols for {index}: {symbols} for Taylor term >= {min_taylor} set to zero.')
                    break # some symbols still in expression, go to next
            if break_flag:
                break

        print(f"The non-simplified coeffs are:")
        print(coefs)
        print(f"The simplified coefficients are:")
        print(coefs_final)

        # calculate FFT of the derived kernel
        # to check if it is better than the best up till now
        tstart = -5
        tend = 5
        steps = 1000
        ttot = tend - tstart
        tstep = ttot/steps
        t = np.linspace(tstart, tend, steps, endpoint=False)
        fft_result = np.fft.fft([kernel_func(s, coefs_final, coefficients) for s in t])
        freq = np.fft.fftfreq(len(t), tstep)
        # shift
        fft_result_shifted = np.fft.fftshift(fft_result)
        freq_shifted = np.fft.fftshift(freq)
        fft_result_normalized = np.real(fft_result_shifted)/(np.max(np.abs(fft_result_shifted)))
        fft_result_normalized_abs = np.abs(fft_result_shifted)/(np.max(np.abs(fft_result_shifted)))
        # calculate indices and values at important points to ensure that the kernel is good
        index_zero = find_closest(freq_shifted, 0.0)
        index_near_zero = find_closest(freq_shifted, -0.25)
        index_near_near_zero = find_closest(freq_shifted, -0.25/2)
        index_half = find_closest(freq_shifted, -0.5)
        zero_val = fft_result_normalized_abs[index_zero]
        near_zero_val = fft_result_normalized_abs[index_near_zero]
        near_near_zero_val = fft_result_normalized_abs[index_near_near_zero]
        # calculate approximate derivative at cutoff frequency (do not take absolute value, since this derivative may be negative)
        derivative = (fft_result_normalized[index_half+1]-fft_result_normalized[index_half-1])/(freq_shifted[index_half+1]-freq_shifted[index_half-1])
        print(f"Current derivative: {derivative:.3f}, best derivative: {best_derivative:.3f}")

        # update best solution if criteria are met
        if derivative > best_derivative and math.isclose(zero_val, 1.0, rel_tol=0.05) and \
                    math.isclose(near_near_zero_val, 1.0, rel_tol=0.075) and \
                    math.isclose(near_zero_val, 1.0, rel_tol=0.1):
            print(f"New best combination at {combination}!")
            abs_derivative = np.abs(fft_result_normalized[index_half+1]-fft_result_normalized[index_half-1])/np.abs(freq_shifted[index_half+1]-freq_shifted[index_half-1])
            print(f"Derivative approximately: {abs_derivative:.3f}")
            best_derivative = derivative
            fft_grads[combi_count-1] = derivative
            zeros_vals[combi_count-1] = zero_val
            coefs_all[combi_count-1] = coefs_final
        
#%%
### Lastly, unpack the best combination and plot the kernel and its frequency response
large_grad = 0.0
for (key, grad) in fft_grads.items():
    if grad > large_grad:
        large_grad = grad
        fft_grad_best = fft_grads[key]
        coefs_best = coefs_all[key]
        key_best = key

print(f"Best derivative: {fft_grad_best} at {key_best}")
print(f'Final coefs: {coefs_best}')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
r = np.linspace(-4, 4, 10_000)
quint = [kernel_func(s, coefs_best, coefficients) for s in r]

# plot kernel in time domain (relative to t=0)
ax1.plot(r, quint)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Amplitude')
ax1.set_title('Kernel in time domain')
ax1.set_xlim(-4, 4)
ax1.legend([f"{p_deg}rd degree kernel"])
ax1.grid(True)

# plot kernel in frequency domain
fs = 10_000 # sampling frequency
tmax = 100.0
t = np.arange(-tmax, tmax, 1/fs)
f = [kernel_func(s, coefs_best, coefficients) for s in t]
k_shifted = np.fft.fftshift(np.fft.fftfreq(len(t), 1/fs))
fft_shifted = np.fft.fftshift(np.fft.fft(f))
fft_result_normalized = np.abs(fft_shifted)/np.max(np.abs(fft_shifted))
index_half = find_closest(k_shifted, -0.5)
derivative = -np.abs(fft_result_normalized[index_half+1]-fft_result_normalized[index_half-1])/(k_shifted[index_half+1]-k_shifted[index_half-1])
ax2.plot(k_shifted, fft_result_normalized)
ax2.set_title(label=f"Frequency response, derivative at f=0.5: {derivative:.3f}")
ax2.set_xlim(0, 5)
freq_ideal = np.linspace(-5.0, 5.0, 1000)
ideal_signal = ideal(freq_ideal)
ax2.plot(freq_ideal, ideal_signal)
ax2.set_yscale('log')
ax2.set_ylim(1e-7, 2)
ax2.set_xlim(-5,5)
ax2.set_xlabel("Frequency [Hz]")
ax2.set_ylabel("Magnitude")
ax2.grid(True)
ax2.legend([f"FFT of {p_deg}rd degree kernel", "Ideal (sinc-function)"])
plt.tight_layout()
plt.show()
# %%
