import numpy as np

def A_express(x):
    
    return 1 + (2.2*((x-1.5)**2))

def p_express(x):
    
    return 1 - (0.3146*x)

def T_express(x):
    
    return 1 - (0.2314*x)

def V_express(x, T):
    
    return (0.1 + 1.09*x)*(np.sqrt(T))

def dp_dt(p, V, A, n, dx):
    
    result = np.zeros((n))
    for i in range(1, n-1):
        
        result[i] = (-p[i]*((V[i+1] - V[i])/dx))\
                      - (p[i]*V[i]*((np.log(A[i+1]) - np.log(A[i]))/dx))\
                      -(V[i]*((p[i+1] - p[i])/dx))

    return result

def dV_dt(p, V, T, n, dx, lam):
    
    result = np.zeros((n))
    for i in range(1, n-1):
        
        result[i] = (-V[i]*((V[i+1] - V[i])/dx))\
                     - ( (1/lam) * ( ((T[i+1] - T[i])/dx) + ( (T[i]/p[i]) * ((p[i+1]-p[i])/dx) ) ) )

    return result        
              
def dT_dt(V, T, A, n, dx, lam):
    
    result = np.zeros((n))
    for i in range(1, n-1):
        
        result[i] = (-V[i]*((T[i+1]-T[i])/dx))\
                     -( ((lam - 1) * T[i]) * ( ((V[i+1] - V[i])/dx) + (V[i]*(np.log(A[i+1]) - np.log(A[i]))/dx) ) )
    
    return result

def passo_preditor(p, V, T, A, n, dx, dt, lam):
    
    dp_dt_ant = dp_dt(p, V, A, n, dx)
    dV_dt_ant = dV_dt(p, V, T, n, dx, lam)
    dT_dt_ant = dT_dt(V, T, A, n, dx, lam)
    
    p_ = np.zeros((n))
    V_ = np.zeros((n))
    T_ = np.zeros((n))
    
    p_[0] = p[0]
    V_[0] = V[0]
    T_[0] = T[0]
    
    p_[-1] = p[-1]
    V_[-1] = V[-1]
    T_[-1] = T[-1]
    
    for i in range(1, n-1):
        
        p_[i] = p[i] + (dp_dt_ant[i] * dt)
        V_[i] = V[i] + (dV_dt_ant[i] * dt)
        T_[i] = T[i] + (dT_dt_ant[i] * dt)
        
    return dp_dt_ant, dV_dt_ant, dT_dt_ant, p_, V_, T_

def passo_corretor(p_ant, V_ant, T_ant, dp_dt_ant, dV_dt_ant, dT_dt_ant, p_, V_, T_, A, n, dx, dt, lam):
    
    dp_dt_atual = dp_dt(p_, V_, A, n, dx)
    dV_dt_atual = dV_dt(p_, V_, T_, n, dx, lam)
    dT_dt_atual = dT_dt(V_, T_, A, n, dx, lam)
    
    p = np.zeros((n))
    V = np.zeros((n))
    T = np.zeros((n))
    
    for i in range(1, n-1):
        
        dp_dt_avg_i = 0.5*(dp_dt_ant[i] + dp_dt_atual[i])
        dV_dt_avg_i = 0.5*(dV_dt_ant[i] + dV_dt_atual[i])
        dT_dt_avg_i = 0.5*(dT_dt_ant[i] + dT_dt_atual[i])
        
        #p[t+dt], V[t+dt], T[t+dt]
        p[i] = p_ant[i] + (dp_dt_avg_i * dt)
        V[i] = V_ant[i] + (dV_dt_avg_i * dt)
        T[i] = T_ant[i] + (dT_dt_avg_i * dt)
    
    return p, V, T
    
def solver(dominiox, n, nt, lam, courant, dx, dt):
    
    A = np.array([A_express(dominiox[0] + i*dx) for i in range(n)])
    
    p_atual = np.array([p_express(dominiox[0] + i*dx) for i in range(n)])
    p_ant = np.array([p_express(dominiox[0] + i*dx) for i in range(n)])
    
    T_atual = np.array([T_express(dominiox[0] + i*dx) for i in range(n)])
    T_ant = np.array([T_express(dominiox[0] + i*dx) for i in range(n)])
    
    
    V_atual = np.array([V_express(dominiox[0] + i*dx, T_atual[i]) for i in range(n)])
    V_ant = np.array([V_express(dominiox[0] + i*dx, T_atual[i]) for i in range(n)])
    
    print("it: 0")
    printTabela(dominiox, n, dx, A, p_atual, V_atual, T_atual)
    
    for t in range(nt):
        
        p_atual, p_ant = p_ant, p_atual
        
        #passo preditor para calcular dp_dt, dV_dt, dT_dt em t
        #e os valores previstos p_, V_ e T_ em t + dt
        dp_dt_ant, dV_dt_ant, dT_dt_ant, p_, V_, T_  = \
                            passo_preditor(p_ant, V_ant, T_ant, A, n, dx, dt, lam)
    
        #passo corretor para calcular p, V e T em t+dt
        #com p_, V_ e T_ calculados no passo preditor
        p, V, T = passo_corretor(p_ant, V_ant, T_ant, dp_dt_ant, dV_dt_ant, dT_dt_ant,\
                       p_, V_, T_, A, n, dx, dt, lam)
        
        #contorno no ponto 1 (indice 0)
        p[0] = 1
        T[0] = 1
        V[0] = (2*V[1]) - V[2]
        
        #contorno no ponto supersonico N (indice -1)
        V[-1] = (2*V[-2] - V[-3])
        p[-1] = (2*p[-2] - p[-3])
        T[-1] = (2*T[-2] - T[-3])
        
        p_atual[:] = p[:]
        V_atual[:] = V[:]
        T_atual[:] = T[:]
    
        print("it: " + str(t+1))
        printTabela(dominiox, n, dx, A, p_atual, V_atual, T_atual)
    
    
    return p_atual, V_atual, T_atual

def printTabela(dominiox, n, dx, A, p, V, T):
    
    print("\n\n\nRESULTADOS COMPARACAO COM TABELA 7.3 LIVRO DO Anderson\n")
    print("____________________________________________________________________")
    print("|    x/L    |    A/A*    |    p/p_0    |    V/a_0    |    T/T_0    |")
    print("____________________________________________________________________")
    for i in range(n):
        
        print("|    {:.3f}".format(dominiox[0] + i*dx) + "  |"\
              + "    {:.3f}".format(A[i]) + "   |"\
              + "    {:.3f}".format(p[i]) + "    |"\
              + "    {:.3f}".format(V[i]) + "    |"\
              + "    {:.3f}".format(T[i]) + "    |")
    print("____________________________________________________________________\n\n\n")

def main():
    
    n = 31
        
    dominiox = [0, 3]

    courant = 0.5

    lam = 1.4
    
    nt = 1400
    
    dx = ((dominiox[1] - dominiox[0])/(n-1))
    
    dt = courant*dx
    dt = 2.94e-3
    
    solver(dominiox, n, nt, lam, courant, dx, dt)
    
main()