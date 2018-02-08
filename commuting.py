class Polynomial:
    def __init__(self, qcwa, monomials):
        assert all([m.qcwa == qcwa for m in monomials])
        self.qcwa = qcwa
        self.monomials = [Monomial.copy_mono(m) for m in monomials]
        
    def simplify (self):
        for m1 in self.monomials:
            for m2 in self.monomials:
                if m1.indices == m2.indices and m1.q_power == m2.q_power and m1 != m2:
                    m1.coefficient += m2.coefficient
                    m2.coefficient = 0
        self.monomials = [m for m in self.monomials if m.coefficient != 0]
    
    def __str__(self):
        toret = ""
        for m in self.monomials:
            toret += str(m)
            if m != self.monomials[-1]:
                toret += " + "
        if toret == "":
            toret = "0"
        return toret

class Monomial:
    def __init__(self, qcwa, indices, coefficient, q_power):
        self.qcwa = qcwa
        self.indices = indices
        self.coefficient = coefficient
        self.q_power = q_power

    @staticmethod
    def copy_mono(monomial):
        return Monomial(monomial.qcwa, list(monomial.indices), monomial.coefficient, monomial.q_power)

    def __str__(self):
        toret = ""
        if self.coefficient != 1:
            if self.coefficient == -1:
                toret += "-"
            else:
                toret += str(self.coefficient)
        if self.q_power != 0:
            if self.q_power != 1:
                toret += "q^"
                if self.q_power < 0:
                    toret += "{"
                toret += str(self.q_power)
                if self.q_power < 0:
                    toret += "}"
            else:
                toret += "q"
        for i in range(self.qcwa.num_variables):
            if self.indices[i] != 0:
                toret += "x_"
                toret += str(i+1)
                if self.indices[i] != 1:
                    toret += "^"
                    toret += str(self.indices[i])
        if toret == "" or toret == "-":
            toret += "1"
        return toret

    def to_poly(self):
        return Polynomial(self.qcwa, [self])

class QuantisedConnectedWeylAlgebra:
    def __init__(self, num_variables, q_matrix, r_matrix):
        assert all([q_matrix[i][j] == -q_matrix[j][i] for i in range(num_variables) for j in range(num_variables)])
#        assert all([
        self.num_variables = num_variables
        self.q_matrix = q_matrix
        self.r_matrix = r_matrix

    def sum(self, polynomials):
        assert all([p.qcwa == self for p in polynomials])
        summed = Polynomial(self, [m for p in polynomials for m in p.monomials])
        summed.simplify()
        return summed
    
    def multiply(self, p1, p2):
        assert p1.qcwa == self
        assert p2.qcwa == self
        expanded_brackets = []
        for m1 in p1.monomials:
            for m2 in p2.monomials:
                expanded_brackets.append(self.multiply_monomials(m1, m2))
        return self.sum(expanded_brackets)

    def multiply_list(self, polynomials, coefficient=1, q_power=0):
        toret = Polynomial(self, [Monomial(self, [0]*self.num_variables, coefficient, q_power)])
        for p in polynomials:
            toret = self.multiply(toret, p)
        return toret

    def x_i(self, i, power=1):
        xilist = [0]*self.num_variables
        xilist[i] = power
        return Monomial(qcwa, xilist, 1, 0)

    def multiply_monomials(self, m1, m2):
        if all([m1.indices[x] == 0 for x in range(self.num_variables)]) or all ([m2.indices[x] == 0 for x in range(self.num_variables)]): 
            indices = []
            for x in range(self.num_variables):
                indices.append(m1.indices[x] + m2.indices[x])
            return Monomial(self, indices, m1.coefficient * m2.coefficient, m1.q_power + m2.q_power).to_poly()
        else:
            i = max([x for x in range(self.num_variables) if m1.indices[x] != 0])
            j = min([x for x in range(self.num_variables) if m2.indices[x] != 0])
            if i > j:
                left = Monomial.copy_mono(m1)
                left.indices[i] -= 1
                right = Monomial.copy_mono(m2)
                right.indices[j] -= 1
                xi = self.x_i(i)
                xj = self.x_i(j)
                
                term1 = self.multiply_list([left.to_poly(), xj.to_poly(), xi.to_poly(), right.to_poly()], q_power=self.q_matrix[i][j])
                if self.r_matrix[i][j] != 0:
                    term2 = self.multiply_list([left.to_poly(), right.to_poly()])
                    term3 = self.multiply_list([left.to_poly(), right.to_poly()], -1, self.r_matrix[i][j])
                    return self.sum([term1, term2, term3])
                else:
                    return term1
                
            else:
                indices = []
                for x in range(self.num_variables):
                    indices.append(m1.indices[x] + m2.indices[x])
                return Monomial(self, indices, m1.coefficient * m2.coefficient, m1.q_power + m2.q_power).to_poly()

    def commutator(self, p1, p2, q_power=0):
        return self.sum([self.multiply_list([p1, p2]), self.multiply_list([p2, p1], -1, q_power)])

class CyclicQCWA(QuantisedConnectedWeylAlgebra):
    def __init__(self, n):
        assert n % 2 == 1
        q_matrix = [[0 for i in range(n)] for j in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    if ((i + j) % 2 == 0 and i < j) or ((i + j) % 2 == 1 and i > j):
                        q_matrix[i][j] = -1
                    else:
                        q_matrix[i][j] = 1
        
        r_matrix = [[0 for i in range(n)] for j in range(n)]
        for i in range(n-1):
            r_matrix[i][i+1] = 1
            r_matrix[i+1][i] = -1
        r_matrix[n-1][0] = 1
        r_matrix[0][n-1] = -1
        super().__init__(n, q_matrix, r_matrix)

    def z_i(self, i):
        assert i < self.num_variables
        if i == 0:
            return Monomial(self, [0]*self.num_variables, 1, 0).to_poly()
        elif i == 1:
            return self.x_i(0).to_poly()
        else:
            zim1xi = self.multiply(self.z_i(i-1), self.x_i(i-1).to_poly())
            zim2 = self.multiply_list([self.z_i(i-2)], -1)
            return self.sum([zim1xi, zim2])

    def theta(self, poly):
        monos = [Monomial.copy_mono(m) for m in poly.monomials]
        newpolys = []
        for m in monos:
            m.indices.insert(0, 0)
            x1power = m.indices.pop()
            x1 = self.x_i(0, x1power)
            newpolys.append(self.multiply(m.to_poly(), x1.to_poly()))
        return self.sum(newpolys)

    def theta_j(self, poly, j):
        while j < 0:
            j += self.num_variables
        p = poly
        for r in range(j):
            p = self.theta(p)
        return p

    def t_i(self, i):
        polys = []
        for li in self.t_i_iterator(i):
            polys.append(self.list_to_poly(li))
        return self.sum(polys)

    def t_i_iterator(self, i):
        assert i % 2 == 1
        if i == self.num_variables:
            yield range(self.num_variables)
        else:
            for li in self.t_i_iterator(i+2):
                for j in range(self.num_variables):
                    if j in li and ((j == self.num_variables - 1 and 0 in li) or (j != self.num_variables - 1 and j+1 in li)):
                        toret = list(li)
                        toret.remove(j)
                        if j+1 in toret:
                            toret.remove(j+1)
                        elif 0 in toret:
                            toret.remove(0)
                        yield toret

    def list_to_poly(self, li):
        i = 0
        looped = False
        if li[-1] == self.num_variables - 1:
            while i < len(li) and li[i] == i:
                i += 1
        polys = []
        while not looped:
            base = i
            while i < len(li) and li[i] == i - base + li[base]:
                i += 1
            if i == len(li):
                looped = True
                if li[-1] == self.num_variables - 1:
                    i = 0
                    while i == li[i]:
                        i += 1
            length = li[i-1] - li[base] + 1
            if length < 0:
                length += self.num_variables
            polys.append(self.theta_j(self.z_i(length), li[base]))
        return self.multiply_list(polys)

for n in range(7, 15, 2):
    qcwa = CyclicQCWA(n)
    print ("n = " + str(n) + ":")
    for i in range(1, n, 2):
        for j in range(i + 2, n, 2):
            print ("commutator of t_" + str(i) + " with t_" + str(j) + ": ")
            print (qcwa.commutator(qcwa.t_i(i), qcwa.t_i(j)))

