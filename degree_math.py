import math
import copy

# ==========================================
# РАЗДЕЛ 1, 2, 3 & 4: БАЗОВЫЕ ГРАДУСНЫЕ ЧИСЛА
# ==========================================

class DegreeNumber:
    """
    Класс градусного числа (a*).
    Автоматически поддерживает Закон Единичного Градуса (a* = a^1*),
    Мнимый Градус 'i' для извлечения корней из отрицательных пространств
    и Закон Ориентированного Нуля.
    """
    def __init__(self, value, is_degree=False, is_imaginary=False):
        # Если передано классическое отрицательное число, конвертируем в градус
        if isinstance(value, (int, float)) and value < 0:
            value = abs(value)
            is_degree = not is_degree
            
        self.value = value
        self.is_degree = is_degree
        self.is_imaginary = is_imaginary

    def to_classical(self):
        """Внутренняя конвертация в комплексные/вещественные числа Python"""
        if self.is_imaginary:
            mult = -1 if self.is_degree else 1
            return complex(0, self.value * mult)
        return -self.value if self.is_degree else self.value

    def __add__(self, other):
        if not isinstance(other, DegreeNumber): other = DegreeNumber(other)
        return DegreeNumber(self.to_classical() + other.to_classical())

    def __sub__(self, other):
        if not isinstance(other, DegreeNumber): other = DegreeNumber(other)
        return DegreeNumber(self.to_classical() - other.to_classical())

    def __mul__(self, other):
        if not isinstance(other, DegreeNumber): other = DegreeNumber(other)
        new_val = self.value * other.value
        
        # Цикличность маркера: i * i = *
        if self.is_imaginary and other.is_imaginary:
            return DegreeNumber(new_val, is_degree=not (self.is_degree ^ other.is_degree ^ True))
        
        return DegreeNumber(new_val, is_degree=self.is_degree ^ other.is_degree, 
                            is_imaginary=self.is_imaginary ^ other.is_imaginary)

    def __truediv__(self, other):
        if not isinstance(other, DegreeNumber): other = DegreeNumber(other)
        
        # Закон Ориентированного Нуля (Раздел 5)
        if other.value == 0:
            return DegreeNumber(float('inf'), is_degree=self.is_degree ^ other.is_degree)
            
        new_val = self.value / other.value
        return DegreeNumber(new_val, is_degree=self.is_degree ^ other.is_degree,
                            is_imaginary=self.is_imaginary ^ other.is_imaginary)

    def __pow__(self, other):
        if not isinstance(other, DegreeNumber): other = DegreeNumber(other)
        # Синтез через Закон Единичного Градуса
        base_p = DegreeNumber(1, is_degree=self.is_degree)
        exp_p = DegreeNumber(other.value, is_degree=other.is_degree)
        total_exp = base_p * exp_p
        
        return DegreeNumber(self.value ** total_exp.value, is_degree=total_exp.is_degree)

    def sqrt(self):
        """Закон Радикала: √a* = √a * i"""
        res_val = math.sqrt(self.value)
        if self.is_degree:
            return DegreeNumber(res_val, is_degree=False, is_imaginary=True)
        return DegreeNumber(res_val)

    def sin(self):
        return DegreeNumber(math.sin(self.value), is_degree=self.is_degree)

    def cos(self):
        return DegreeNumber(math.cos(self.value), is_degree=False) # Поглощает градус

    def __repr__(self):
        if self.value == float('inf'): return "∞*" if self.is_degree else "∞"
        marker = f"{'i' if self.is_imaginary else ''}{'*' if self.is_degree else ''}"
        return f"{self.value}{marker}"


# ==========================================
# РАЗДЕЛ 5: ЛИНЕЙНАЯ АЛГЕБРА И ПРОСТРАНСТВА
# ==========================================

class DegreeMatrix:
    """
    Градусные Матрицы. Вычисление детерминанта (det*) с автоматическим Шахматным Градусом.
    """
    def __init__(self, matrix_data):
        self.data = [[item if isinstance(item, DegreeNumber) else DegreeNumber(item) for item in row] for row in matrix_data]
        self.rows = len(self.data)
        self.cols = len(self.data[0]) if self.rows > 0 else 0

    def det(self):
        if self.rows != self.cols: raise ValueError("Матрица должна быть квадратной!")
        if self.rows == 1: return self.data[0][0]
        
        # Формула 2x2: (a*d) + (b*c)*
        if self.rows == 2:
            a, b, c, d = self.data[0][0], self.data[0][1], self.data[1][0], self.data[1][1]
            return (a * d) + (b * c * DegreeNumber(1, is_degree=True))

        # Развертывание 3x3 по первой строке
        if self.rows == 3:
            total_det = DegreeNumber(0)
            for j in range(3):
                element = self.data[0][j]
                minor = DegreeMatrix([[self.data[r][c] for c in range(3) if c != j] for r in range(1, 3)])
                term = element * minor.det()
                if j % 2 != 0:  # Шахматный Флип
                    term = term * DegreeNumber(1, is_degree=True)
                total_det = total_det + term
            return total_det
        raise NotImplementedError("Размерности выше 3x3 пока не задекларированы.")

    def __repr__(self):
        return "\n".join(["\t".join([str(item) for item in row]) for row in self.data])


class DegreeVector:
    """Векторное Пространство градусных векторов (Degree Vectors)."""
    def __init__(self, x, y, z=0):
        self.x = x if isinstance(x, DegreeNumber) else DegreeNumber(x)
        self.y = y if isinstance(y, DegreeNumber) else DegreeNumber(y)
        self.z = z if isinstance(z, DegreeNumber) else DegreeNumber(z)

    def dot(self, other):
        """Скалярное произведение"""
        return (self.x * other.x) + (self.y * other.y) + (self.z * other.z)

    def cross(self, other):
        """Векторное произведение с градусным разворотом"""
        rx = (self.y * other.z) + (self.z * other.y * DegreeNumber(1, is_degree=True))
        ry = ((self.x * other.z) + (self.z * other.x * DegreeNumber(1, is_degree=True))) * DegreeNumber(1, is_degree=True)
        rz = (self.x * other.y) + (self.y * other.x * DegreeNumber(1, is_degree=True))
        return DegreeVector(rx, ry, rz)

    def __repr__(self):
        return f"({self.x}, {self.y}, {self.z})"


def solve_cramer(matrix_A, vector_B):
    """Решение СЛАУ методом Крамера в градусном пространстве"""
    det_main = matrix_A.det()
    results = []
    for j in range(matrix_A.cols):
        modified_data = copy.deepcopy(matrix_A.data)
        for i in range(matrix_A.rows):
            modified_data[i][j] = vector_B[i] if isinstance(vector_B[i], DegreeNumber) else DegreeNumber(vector_B[i])
        results.append(DegreeMatrix(modified_data).det() / det_main)
    return results


# ==========================================
# РАЗДЕЛ 5.7 & 5.8: МАТЕМАТИЧЕСКИЙ АНАЛИЗ
# ==========================================

def diff_degree_pow(base_str, power_degree_num):
    """Символьное дифференцирование: (x^(n*))' = (n * x^(n-1))*"""
    return f"({power_degree_num.value} * {base_str}^{power_degree_num.value - 1})*"


class DegreeIntegral:
    """Градусное Интегрирование (аналитическое и численное)"""
    @staticmethod
    def integrate_pow(base_str, power_degree_num):
        """Символьное первообразное: ∫(x^(n*))dx = ((x^(n+1))/(n+1))*"""
        new_exp = power_degree_num.value + 1
        if new_exp == 0:
            return f"(ln|{base_str}|)*"
        return f"({base_str}^{new_exp} / {new_exp})*"

    @staticmethod
    def definite_integral_vector(vector_func, limit_a, limit_b, steps=1000):
        """Определенный интеграл для траектории градусного вектора"""
        a = limit_a.to_classical() if isinstance(limit_a, DegreeNumber) else limit_a
        b = limit_b.to_classical() if isinstance(limit_b, DegreeNumber) else limit_b
        
        dx = (b - a) / steps
        total_vector = DegreeVector(0, 0, 0)
        
        for i in range(steps):
            x = a + i * dx
            v_current = vector_func(x)
            dx_deg = DegreeNumber(dx)
            
            step_vector = DegreeVector(v_current.x * dx_deg, v_current.y * dx_deg, v_current.z * dx_deg)
            total_vector = DegreeVector(
                total_vector.x + step_vector.x,
                total_vector.y + step_vector.y,
                total_vector.z + step_vector.z
            )
        return total_vector


# ==========================================
# РАЗДЕЛ 6: ЛОГАРИФМЫ И АВТОРЕШАТЕЛЬ
# ==========================================

def log_degree(argument, base=None):
    """Градусный Логарифм. Раскрывает логарифмы от градусных (*) пространств."""
    if base is None: base = DegreeNumber(math.e)
    if not isinstance(argument, DegreeNumber): argument = DegreeNumber(argument)
    if not isinstance(base, DegreeNumber): base = DegreeNumber(base)

    if argument.is_degree and not argument.is_imaginary:
        real_part = math.log(argument.value) / math.log(base.value)
        return DegreeNumber(real_part, is_degree=True, is_imaginary=True)
        
    return DegreeNumber(math.log(argument.value) / math.log(base.value))


class DegreeSolver:
    """Универсальный автоматический решатель градусных показательных уравнений (Base^x = Result)"""
    @staticmethod
    def solve_exponential(base, result):
        if not isinstance(base, DegreeNumber): base = DegreeNumber(base)
        if not isinstance(result, DegreeNumber): result = DegreeNumber(result)

        if base.value == 1 and result.value != 1: raise ValueError("Нет решений")
        if base.value == 0 and result.value == 0: return "x — любое число"

        log_res = log_degree(result, base)
        
        if log_res.is_imaginary and log_res.is_degree:
            return DegreeNumber(log_res.value, is_degree=base.is_degree ^ result.is_degree)
            
        if base.is_degree:
            return DegreeNumber(log_res.value, is_degree=base.is_degree ^ result.is_degree)
            
        return log_res

# ==========================================
# ПРОВЕРОЧНЫЙ ТЕСТ-ДРАЙВ СИСТЕМЫ
