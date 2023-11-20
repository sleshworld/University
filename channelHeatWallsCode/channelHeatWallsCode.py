import numpy as np

# Решение трехдиагональной матричной системы
def solve_tridiagonal_system(num_points, diag_below, diag_main, diag_above, rhs):
    # Прямой проход
    diag_above[0] /= diag_main[0]
    rhs[0] /= diag_main[0]
    diag_main[0] = 1.0

    for i in range(1, num_points):
        factor = diag_main[i] - diag_below[i] * diag_above[i-1]
        diag_above[i] /= factor
        rhs[i] = (rhs[i] - diag_below[i] * rhs[i-1]) / factor
        diag_main[i] = 1.0

    # Обратный проход
    for i in range(num_points-2, -1, -1):
        rhs[i] -= diag_above[i] * rhs[i+1]

    return rhs

# Правило трапеций для интегрирования
def trapezoidal_integration(num_intervals, x_values, y_values, step):
    integral_sum = 0.0
    for i in range(1, num_intervals):
        integral_sum += (y_values[i-1] + y_values[i]) * (x_values[i] - x_values[i-1]) / step
    return integral_sum / 2.0

# Функция интерполяции для вязкости
def viscosity_interpolation(temperature):
    # Коэффициенты получены в Excel
    return 3.381e-11 * temperature**4 - 9.304e-9 * temperature**3 + 9.963e-7 * temperature**2 - 5.531e-5 * temperature + 1.781e-3

# Вывод результатов в файл
def output_to_file(filename, num_points, x, velocities, analytical_velocities, temperatures, interpolated_viscosities):
    with open(filename, 'w') as file:
        file.write('VARIABLES = "X","U","U_AN","T","MU_Interpolation"\n')
        file.write(f'ZONE I={num_points}, DATAPACKING=BLOCK\n')
        np.savetxt(file, x)
        np.savetxt(file, velocities)
        np.savetxt(file, analytical_velocities)
        np.savetxt(file, temperatures)
        np.savetxt(file, interpolated_viscosities)

# Основная функция программы
def main():
    num_points = 21
    flow_scheme = 2
    mu_def = 2 # способ определения вязкости для неизотермич случая 1 или 2
    H = 0.02
    temp_dif = 40.0
    T1 = 45.0 - temp_dif
    T2 = 45.0 + temp_dif

    # Инициализация массивов
    velocities = np.zeros(num_points) 
    x = np.zeros(num_points)
    diag_below = np.zeros(num_points)
    diag_main = np.zeros(num_points)
    diag_above = np.zeros(num_points)
    rhs = np.zeros(num_points)
    analytical_velocities = np.zeros(num_points)
    temperatures = np.zeros(num_points)
    interpolated_viscosities = np.zeros(num_points)

    # Расчет шага по позициям и определение позиций
    position_step = H / (num_points - 1)
    x = np.linspace(0, H, num_points)

    
    if flow_scheme == 1: # изотермическая
        mu_default = 0.599e-3 # Выбор коэффициента вязкости при температуре 45 градусов
        velocity_boundary = 1000.0 * mu_default / (1000.0 * H) # Расчет аналитической среднерасходной скорости
        constant = 12.0 * 1000.0 * velocity_boundary**2 / (1000.0 * H) # расчет градиента давления
        # коэффициенты для метода прогонки
        diag_below[1:-1] = mu_default / position_step**2 
        diag_main[1:-1] = -2 * mu_default / position_step**2
        diag_above[1:-1] = mu_default / position_step**2
        rhs[1:-1] = -constant

        # Граничные условия
        diag_below[0], diag_main[0], diag_above[0], rhs[0] = 0.0, 1.0, 0.0, 0.0
        diag_below[-1], diag_main[-1], diag_above[-1], rhs[-1] = 0.0, 1.0, 0.0, 0.0

        velocities = solve_tridiagonal_system(num_points, diag_below, diag_main, diag_above, rhs)
        velocity_average_analytical = trapezoidal_integration(num_points, x, velocities, H)
        print(f'Аналитическая средневзвешенная скорость: {velocity_boundary}, Средневзвешенная скорость: {velocity_average_analytical}')
        analytical_velocities = -constant / (2 * mu_default) * x**2 + constant * H / (2 * mu_default) * x


    
    if flow_scheme == 2: # неизотермическая
        mu_default = 0.599e-3
        velocity_boundary = 1000.0 * mu_default / (1000.0 * H)
        constant = 12.0 * 1000.0 * velocity_boundary**2 / (1000.0 * H)

        temperatures = T1 + (T2 - T1) * x / H
        interpolated_viscosities = np.vectorize(viscosity_interpolation)(temperatures)
        
        for j in range(1, num_points - 1):
            if mu_def == 1:
                mu_plus = (interpolated_viscosities[j + 1] + interpolated_viscosities[j]) / 2.0
                mu_minus = (interpolated_viscosities[j] + interpolated_viscosities[j - 1]) / 2.0
            else:
                x_plus = (x[j+1]+x[j])/2.0
                x_minus = (x[j]+x[j-1])/2.0
                t_plus = T1 + (T2 - T1) * x_plus/H
                t_minus = T1 + (T2 - T1) * x_minus/H
                mu_plus = viscosity_interpolation(t_plus)
                mu_minus = viscosity_interpolation(t_minus)
            diag_below[j] = mu_minus / position_step**2
            diag_main[j] = -(mu_minus + mu_plus) / position_step**2
            diag_above[j] = mu_plus / position_step**2
            rhs[j] = -constant

        # Граничные условия
        diag_below[0], diag_main[0], diag_above[0], rhs[0] = 0.0, 1.0, 0.0, 0.0
        diag_below[-1], diag_main[-1], diag_above[-1], rhs[-1] = 0.0, 1.0, 0.0, 0.0

        velocities = solve_tridiagonal_system(num_points, diag_below, diag_main, diag_above, rhs)
        velocity_average_analytical = trapezoidal_integration(num_points, x, velocities, H)
        print('Средневзвешенная скорость:', velocity_average_analytical)

        # Расчет коэффициента гидравлического сопротивления
        hydraulic_resistance_coefficient = 2 * constant * H / (1000 * velocity_average_analytical**2)
        print('Коэффициент гидравлического сопротивления:', hydraulic_resistance_coefficient)

    # Вывод результатов в файл
    output_to_file(f'Results_{flow_scheme}_{num_points}_deltaT_{temp_dif}.plt', num_points, x, velocities, analytical_velocities, temperatures, interpolated_viscosities)

if __name__ == "__main__":
    main()