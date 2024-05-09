#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from gurobipy import *
import numpy as np
import math
import time
from numpy import linalg as LA


# In[2]:


# f = open("1-4-multi-100-1.txt", "r")
# print(f.read())


# In[3]:


# print(f.read())


# In[4]:


# f = open("1-4-multi-100-1.txt", "r")
# print(f.read())

# Set initial values
epsilon = 0.20
theta = 0.10


# In[5]:


filename = '1-7-1-500-2.txt'
f = open("1-7-1-500-2.txt", "r")
a = []
for line in f:
    a.append(line)

# # The data is of the list type.  The Python list type is actually
# # a dynamic array. The lines contain also the \n; hence the .rstrip()
# for n, line in enumerate(data, 1):
#     print '{:2}.'.format(n), line.rstrip()

# print '-----------------'


# In[6]:


# a


# In[7]:


c = a[0]
c = c.replace("[", "")
c = c.replace("]\n", "")
c = c.split(",")
c = np.array(c)
c = c.astype(np.float64)


# In[8]:


b = a[-1]
b = b.replace("[", "")
b = b.replace("]\n", "")
b = b.split(",")
b = np.array(b)
b = b.astype(np.float64)


# In[9]:


xi = [None]*(len(a)-2)


# In[10]:


for i in range(1,len(a)-1):
    dd = a[i]
    dd = dd.replace("[", "")
    dd = dd.replace("],\n", "")
    dd = dd.replace("]]\n", "")
    dd = dd.split(",")
    dd = np.array(dd)
    dd = dd.astype(np.float64)
    xi[i-1] = dd
#     xi[i-1] = np.around(dd, decimals=2)
    
    


# In[11]:


xi


# In[12]:


np.size(c)


# In[13]:


np.size(b)


# In[14]:


np.size(xi)


# In[15]:


np.size(xi[0])


# In[16]:


random_size = np.size(b)
x_random_size = np.size(c)


# In[17]:


k=math.floor(random_size*epsilon)


# In[18]:


c = -c


# In[ ]:


# CVaR Method

def Model4():
    #Model
    print ("Begin to solve Model CVaR")

    m = Model()
    #m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
    s = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="s")

    gamma = m.addVar(lb=0,name="gamma")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
    m.update()
        
    
    # Add objective    
    m.setObjective(sum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    # Add constraints
#     m.addConstrs(  s[i] >= beta for i in range(random_size))
    m.update()
    m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()

   
    m.addConstrs( gamma + s[i] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i]   for i in range(random_size))
    m.update()
    m.addConstr( theta*slack_t -(epsilon)*gamma - s.sum()/float(random_size)  <=0)
    m.update()
    # Solve the problem
    m.update()
        
    m.optimize()
    

#     Store solutions
#     ppp = m.getAttr('x', beta)
    CVAR_s = m.getAttr('x', s)
    CVAR_x = m.getAttr('x', x)
    

    aaaa=m.objVal
    for v in m.getVars():
        print('%s %g' % (v.varName, v.x))
    return m.objVal,CVAR_s,CVAR_x
    
    


# In[ ]:


start=time.time()
CVaR=Model4()
CVAR=CVaR[0]
CVAR_s=CVaR[1]
CVAR_x=CVaR[2]
modeltime_cvar= time.time() - start


# In[ ]:


# quantitle bounds

def Model_q():

    print ("Begin to solve model quantitle relaxtion")

    m = Model()
    
    m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1, name="x")
#     slack_t = m.addVar(lb=0,name="slack_t")
    m.update()
    
    # Set functions
    
    #obj = 0
    g_x_xi = [None] * random_size
    f_x = 0
    g_x = 0
    
#     m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
#     m.addConstr(slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()
    # Set functions
    for i in range(x_random_size):
      
        f_x += c[i]*x[i]
        m.update()
    
#     for i in range(random_size):
#         g_x_si[i] = g_x +       SSS   
    for i in range(random_size):
     #   g_x_si[i] = g_x + si[i] 
        g_x = 0
        for j in range(x_random_size):
            g_x += xi[i][j]*x[j]
            m.update()
        g_x_xi[i] = g_x  
    
    # Add objective    
    m.setObjective(f_x, GRB.MINIMIZE)
    m.update()
    # Add constraints
    
#     con1=m.addConstr(theta/epsilon*slack_t + g_x_xi[0] <= b[0] )
    con1=m.addConstr( g_x_xi[0] <= b[0] )
    # Solve the problem
    m.update()
    m.params.OutputFlag=0    
    m.optimize()
    kkk = m.getAttr('x', x)
    result = [0]*random_size
    obj = 0
    for i in range(x_random_size):
        obj += c[i]*kkk[i]
    result[0]=obj
    
    
    for i in range(1,random_size):   
        m.remove(con1)
#         con1=m.addConstr(theta/epsilon*slack_t + g_x_xi[i]  <= b[i] )
        con1=m.addConstr(g_x_xi[i]  <= b[i] )
        m.update()
        m.params.OutputFlag=0    
        m.optimize()
        kkk = m.getAttr('x', x)
        obj = 0
        for j in range(x_random_size):
            obj += c[j]*kkk[j]
        result[i]=obj
    return result
        



# In[ ]:


start=time.time()
quantile_result=Model_q()
sorted_array=np.sort(quantile_result,-1)
v_q = sorted_array[random_size-k]
print(v_q)
modeltime_quantile = time.time() - start


# In[ ]:


# scenario_value_alsox_sharp = []


# In[ ]:


# for i in range(random_size):
#     scenario_value_alsox_sharp.append(sum(xi[i][j]*solution_alsox_sharp[j] for j in range(x_random_size)))


# In[ ]:


v_q = sorted_array[random_size-k]


# In[ ]:


v_q


# In[ ]:


CVAR


# In[ ]:


violation = math.floor(random_size*epsilon)


# In[ ]:


# DRCC set
def Model_checking(V_x):
    
    m = Model()
    #m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()

#     x = m.addVars(x_random_size,vtype=GRB.BINARY,name="x")
#     x = m.addVars(x_random_size,name="x")

    s = m.addVars(random_size,lb=0,name="s")

    gamma = m.addVar(lb=-GRB.INFINITY,ub=0,name="gamma")
    lambda_1 = m.addVar(lb=0,name="lambda_1")
    m.update()
    

    # Add objective    
    m.setObjective(0, GRB.MINIMIZE)
    m.update()
    
    # Add constraints
   
    m.addConstrs( -max(b[i]-sum(xi[i][j]*V_x[j] for j in range(x_random_size)),0) - gamma  <=   s[i]  for i in range(random_size))
    m.update()
#     m.addConstrs(s[j] >= gamma for j in range(random_size))

    m.addConstr(theta*lambda_1 + s.sum()/float(random_size) + epsilon*gamma <= 0)
    m.update()
    
#     m.addConstr(sum(V_x[j] for j in range(x_random_size))<=lambda_1)
    m.update()
     
    m.addConstrs(lambda_1 >= V_x[j] for j in range(x_random_size))
    
    # Solve the problem
    m.update()
    m.params.OutputFlag=0
        
    m.optimize()
    result = 1
    if  m.status == 4 or m.status == 3 or m.status ==12:
        result = 0
    
   
    return result


# In[ ]:


# CVaR Method

def Model_alsox_sharp():
    #Model
    print ("Begin to solve model ALSOX Sharp")
    delta_1 = 1e-1
    v_lower_upper = v_q
    v_upper_bound = CVAR
    current_bound = (v_lower_upper+v_upper_bound)/2.0
    
    delta_t = v_upper_bound - v_lower_upper
    
    m = Model()
    #m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    s = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="s")

    gamma = m.addVar(lb=0,name="gamma")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
    m.update()
    

    # Add objective    
    m.setObjective(theta*slack_t -(epsilon)*gamma - s.sum()/float(random_size), GRB.MINIMIZE)
    m.update()
    # Add constraints
#     m.addConstrs(  s[i] >= beta for i in range(random_size))
    m.update()
#     m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
    m.addConstrs(slack_t >= x[i] for i in range(x_random_size))
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()

   
    m.addConstrs( gamma + s[i] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i]   for i in range(random_size))
    m.update()

    
    con1=m.addConstr(sum(c[i]*x[i] for i in range(x_random_size)) <= current_bound)
    m.update()
    # Solve the problem
    m.update()
    
    m.params.OutputFlag=0
        
    m.optimize()
        
#     Store solutions
#     ppp = m.getAttr('x', beta)
    current_s = m.getAttr('x', s)
    current_x = m.getAttr('x', x)
    
    
#     support=0
#     for j in current_s:
#         if current_s[j]>1e-6:
#             support +=1
#     if support > violation:
#         v_lower_upper = current_bound
# #         t1=con1.SARHSUp
#     else:
#         v_upper_bound = current_bound
        
    while delta_t >= delta_1:      
        checking = Model_checking(current_x)
        if checking == 1:
            v_upper_bound = current_bound
        else:
            v_lower_upper = current_bound

        current_bound = (v_lower_upper + v_upper_bound)/2.0
        delta_t = v_upper_bound - v_lower_upper
        con1.RHS = current_bound
        for i in range(random_size):
            s[i].start = current_s[i]
        for i in range(x_random_size):
            x[i].start = current_x[i]  
        m.optimize()
        current_s = m.getAttr('x', s)
        current_x = m.getAttr('x', x)

    con1.RHS = v_upper_bound

    for i in range(random_size):
        s[i].start = current_s[i]
    for i in range(x_random_size):
        x[i].start = current_x[i]  

    m.optimize()

    current_s = m.getAttr('x', s)

    current_x = m.getAttr('x', x)
    
    
#     aaaa=m.objVal
    
#     for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
    return v_upper_bound,current_x
    
    


# In[ ]:


start=time.time()    
alsox_sharp = Model_alsox_sharp()
modeltime_alsox_sharp = time.time() - start   
value_alsox_sharp = alsox_sharp[0]
solution_alsox_sharp = alsox_sharp[1]


# In[ ]:


value_alsox_sharp


# In[ ]:





# In[ ]:


#  lower bound improvment
start=time.time()   
delta_1 = 1e-1
current_bound = v_q
# v_upper_bound = CVAR
# current_bound = -20000

delta_t = 10


m = Model()
#m.setParam('TimeLimit', 60*60)
# Create variables
m.update()
# s1 = m.addVars(random_size,name="s1")
x = m.addVars(x_random_size,lb=0,ub=1,name="x")

y = m.addVar(lb=-GRB.INFINITY,name="y")

s = m.addVars(random_size,lb=0,name="s")

beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")


slack_t = m.addVar(lb=0,name="slack_t")

mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")

w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")

# z = m.addVars(random_size,vtype=GRB.BINARY,name="z")

z = m.addVars(random_size,lb=0,ub=1,name="z")

m.update()

m.setObjective(y, GRB.MINIMIZE)
# Add objective 

# m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
m.update()


m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <= beta + b[i]*z[i] + s[i]  for i in range(random_size))
m.update()
m.addConstr( theta*slack_t + s.sum()/float(random_size) + epsilon*beta <=0)
m.update()
# Solve the problem
# m.update()

m.addConstrs( x[i] == mu[i,j] + w[i,j] for i in range(x_random_size) for j in range(random_size) )
m.update()
# m.addConstrs( mu[j,i]<=1-z[j]  for j in range(x_random_size) for i in range(random_size) )
m.addConstrs( w[j,i] <= 1 - z[i]  for j in range(x_random_size) for i in range(random_size) )
m.update()
m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in range(random_size))
m.update()
m.addConstr(z.sum() >= random_size - math.ceil(random_size*epsilon))
m.update()

constr1 = m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + current_bound * (1-z[j])  for j in  range(random_size))
m.update()
constr2 = m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + current_bound * (z[j])  for j in  range(random_size) )
m.update()

m.params.OutputFlag=0
m.params.BarHomogeneous=1 
m.optimize()

# current_s = m.getAttr('x', s)
current_x = m.getAttr('x', x)

dual_bound = m.objVal
print('bound is ', dual_bound)

ite = 0 

while delta_t >= delta_1:
    
    m.remove(m.getConstrs()[-(2*random_size)-1:-1])
#     m.remove(m.getConstrs()[-2])
    delta_t = abs(dual_bound - current_bound)
    
    current_bound = dual_bound
    
#     if m.status == 2:
#         current_bound = dual_bound
#     else:
#         v_lower_upper = current_bound

#     current_bound = (v_lower_upper + v_upper_bound)/2.0
    
#     print('current bound is', current_bound)
    constr1 = m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + current_bound * (1-z[j])  for j in  range(random_size))
    m.update()
    constr2 = m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + current_bound * (z[j])  for j in  range(random_size) )
    m.update()
#     constr1.RHS = current_bound
#     constr2.RHS = current_bound
#     for i in range(random_size):
#         s[i].start = current_s[i]
    for i in range(x_random_size):
        x[i].start = current_x[i]  
    m.params.OutputFlag=0
    print('current bound is', current_bound)
    m.optimize()
    dual_bound = m.objVal
    
    ite = ite + 1 
    
current_bound = dual_bound
print('with iterations', ite)
#     Store solutions
#     ppp = m.getAttr('x', beta)
# first_dual_bound_s = m.getAttr('x', s)
# first_dual_bound_x = m.getAttr('x', x)

# dual_bound = m.objVal

modeltime_second_dual_bound = time.time() - start   


# In[ ]:


lower_bound_second_dual = current_bound


# In[ ]:


lower_bound_second_dual


# In[ ]:


v_q


# In[ ]:


(CVAR-v_q)/abs(v_q)


# In[ ]:


value_alsox_sharp


# In[ ]:


lower_bound_second_dual


# In[ ]:


(value_alsox_sharp-lower_bound_second_dual)/abs(lower_bound_second_dual)


# In[ ]:


value_approximation_bound = value_alsox_sharp


# In[ ]:


sorted_array_indice = [i for i, value in enumerate(quantile_result) if value > value_approximation_bound]


# In[ ]:


cut_z_indice = sorted_array_indice


# In[ ]:


stengthen_difference = [element for element in range(random_size) if element not in cut_z_indice]


# In[ ]:


len(cut_z_indice)


# In[ ]:


subset_indice = np.argsort(quantile_result)[random_size-2*k:random_size]
# subset_indice = np.argsort(scenario_value_alsox_sharp)[random_size-2*k:random_size]
stengthen_difference_update = [element for element in subset_indice if element not in cut_z_indice]
stengthen_difference_update=stengthen_difference_update[::-1]


# In[ ]:


start=time.time()
cut_third_one = []

m = Model()
m.setParam('Seed', 2)
m.setParam('TimeLimit', 60*60)
# Create variables
m.update()
# s1 = m.addVars(random_size,name="s1")
x = m.addVars(x_random_size,lb=0,ub=1,name="x")

y = m.addVar(lb=-GRB.INFINITY,name="y")

s = m.addVars(random_size,lb=0,name="s")

beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")


slack_t = m.addVar(lb=0,name="slack_t")

mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")

w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")

# z = m.addVars(random_size,vtype=GRB.BINARY,name="z")

z = m.addVars(random_size,lb=0,ub=1,name="z")

m.update()


# Add objective    
m.setObjective(y, GRB.MINIMIZE)
m.update()

# m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
# m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
m.update()

m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <= beta + b[i]*z[i] + s[i]  for i in range(random_size))
m.update()

m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <=  b[i]*z[i]  for i in range(random_size))
m.update()

m.addConstr( theta*slack_t + s.sum()/float(random_size) + epsilon*beta <=0)
m.update()

m.addConstrs( x[i] == mu[i,j] + w[i,j] for i in range(x_random_size) for j in range(random_size) )
m.update()

m.addConstrs( w[j,i] <= 1 - z[i]  for j in range(x_random_size) for i in range(random_size) )
m.update()
m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in range(random_size))
m.update()
m.addConstr(z.sum() >= math.ceil(random_size*(1-epsilon))+1)
m.update()

m.addConstr(y >= sum(c[i] * x[i] for i in range(x_random_size)))

m.update()


m.addConstrs(z[i] == 0 for i in cut_z_indice)

m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + value_approximation_bound * (1-z[j])  for j in stengthen_difference)
m.update()
m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + value_approximation_bound * (z[j])  for j in stengthen_difference)
m.update()

m.update()

m.params.OutputFlag=0 

m.optimize()

z_solution_store_1 = m.getAttr('x', z)
mu_solution_store_1 = m.getAttr('x', mu)
w_solution_store_1 = m.getAttr('x', w)
x_solution_store_1 = m.getAttr('x', x)

print('current obj is', m.ObjVal)


m.remove(m.getConstrs()[-(2*len(stengthen_difference))-1:-1])

iteration = 0

for i in range(0,math.ceil(k*0.2)):       


        print(stengthen_difference_update[i])
        
        m.addConstr(z[stengthen_difference_update[i]]==1)
        m.update()

        stengthen_difference_new = [element for element in range(random_size) if element not in cut_z_indice and element!=stengthen_difference_update[i]]

        m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + value_approximation_bound * (1-z[j])  for j in stengthen_difference_new)
        m.update()
        m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + value_approximation_bound * (z[j])  for j in stengthen_difference_new)
        m.update()


        x.start = x_solution_store_1 
        mu.start = mu_solution_store_1
        w.start = w_solution_store_1
        z.start = z_solution_store_1

        m.params.OutputFlag=0 

        m.optimize()

        print('Updated obj is', m.ObjVal)
        
        if m.objVal >= value_approximation_bound:
            
    #     print(m.status)
#         if m.status == 4:
#         if m.status == 12 or m.status == GRB.INFEASIBLE:
            print("Cut Found")
            cut_third_one.append([i])
        
        iteration += 1
        
        m.remove(m.getConstrs()[-(2*len(stengthen_difference_new))-2:-1])


modeltime_fixing_1 = time.time() - start


# In[ ]:


for i in range(len(cut_third_one)):
    cut_z_indice.append(stengthen_difference_update[cut_third_one[i][0]])
    


# In[ ]:


subset_indice = np.argsort(quantile_result)[random_size-2*k:random_size]
stengthen_difference_update = [element for element in subset_indice if element not in cut_z_indice]
stengthen_difference_update=stengthen_difference_update[::-1]


# In[ ]:


start=time.time()
cut_third = []

m = Model()
m.setParam('Seed', 2)
m.setParam('TimeLimit', 60*60)
# Create variables
m.update()
# s1 = m.addVars(random_size,name="s1")
x = m.addVars(x_random_size,lb=0,ub=1,name="x")

y = m.addVar(lb=-GRB.INFINITY,name="y")

s = m.addVars(random_size,lb=0,name="s")

beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")


slack_t = m.addVar(lb=0,name="slack_t")

mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")

w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")

# z = m.addVars(random_size,vtype=GRB.BINARY,name="z")

z = m.addVars(random_size,lb=0,ub=1,name="z")

m.update()


# Add objective    
m.setObjective(y, GRB.MINIMIZE)
m.update()

# m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
# m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
m.update()

m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <= beta + b[i]*z[i] + s[i]  for i in range(random_size))
m.update()

m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <=  b[i]*z[i]  for i in range(random_size))
m.update()

m.addConstr( theta*slack_t + s.sum()/float(random_size) + epsilon*beta <=0)
m.update()

m.addConstrs( x[i] == mu[i,j] + w[i,j] for i in range(x_random_size) for j in range(random_size) )
m.update()

m.addConstrs( w[j,i] <= 1 - z[i]  for j in range(x_random_size) for i in range(random_size) )
m.update()
m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in range(random_size))
m.update()
m.addConstr(z.sum() >= math.ceil(random_size*(1-epsilon))+1)
m.update()

m.addConstr(y >= sum(c[i] * x[i] for i in range(x_random_size)))

m.update()


m.addConstrs(z[i] == 0 for i in cut_z_indice)

m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + value_approximation_bound * (1-z[j])  for j in stengthen_difference)
m.update()
m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + value_approximation_bound * (z[j])  for j in stengthen_difference)
m.update()

m.update()

m.params.OutputFlag=0 

m.optimize()

z_solution_store_1 = m.getAttr('x', z)
mu_solution_store_1 = m.getAttr('x', mu)
w_solution_store_1 = m.getAttr('x', w)
x_solution_store_1 = m.getAttr('x', x)

print('current obj is', m.ObjVal)


m.remove(m.getConstrs()[-(2*len(stengthen_difference))-1:-1])

iteration = 0

for i in range(0,math.ceil(k*0.3)):       

    for j in range(1,2):
        
        print(stengthen_difference_update[i],stengthen_difference_update[i+j])
        
        m.addConstr(z[stengthen_difference_update[i]]+z[stengthen_difference_update[i+j]]==2)
        m.update()

        stengthen_difference_new = [element for element in range(random_size) if element not in cut_z_indice and element!=stengthen_difference_update[i] and element!=stengthen_difference_update[i+j]]

        m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + value_approximation_bound * (1-z[j])  for j in stengthen_difference_new)
        m.update()
        m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + value_approximation_bound * (z[j])  for j in stengthen_difference_new)
        m.update()


        x.start = x_solution_store_1 
        mu.start = mu_solution_store_1
        w.start = w_solution_store_1
        z.start = z_solution_store_1

        m.params.OutputFlag=0 

        m.optimize()

        print('Updated obj is', m.ObjVal)
        
        if m.objVal >= value_approximation_bound:
            
    #     print(m.status)
#         if m.status == 4:
#         if m.status == 12 or m.status == GRB.INFEASIBLE:
            print("Cut Found")
            cut_third.append([i,i+j])
        
        iteration += 1
        
        m.remove(m.getConstrs()[-(2*len(stengthen_difference_new))-2:-1])


modeltime_fixing_2 = time.time() - start


# In[ ]:


value_alsox_sharp


# In[ ]:


cut_third 


# In[ ]:


cut_z_indice 


# In[ ]:


modeltime_fixing = modeltime_fixing_1 + modeltime_fixing_2


# In[19]:


Big_M = []
for i in range(random_size):
    
#     Big_M.append(abs(math.ceil(sum(x for x in xi[i][:] if x > 0) -b[i]))) 
    Big_M.append(max(abs(math.ceil(sum(x for x in xi[i][:] if x > 0) -b[i])),b[i])) 
    
Big_M_coefficient = Big_M 


# In[ ]:


def big_m_find_1(xi_1,xi_2,b_1,b_2):
    
    m = Model()
    
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    
    y = m.addVar(lb=0,name="y")
    
    m.setObjective(sum(xi_1[j]*x[j] for j in range(x_random_size)) - b_1, GRB.MAXIMIZE)
    
    # Add objective 
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()

    m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)
    
    m.update()
    
#     m.addConstrs( theta/epsilon*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in stengthen_difference)
#     m.update()
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon)))
#     m.update()
    
    m.params.OutputFlag=0
    
    m.optimize()
    
    return m.objVal
    
    
    
    
    
    
    
    
    


# In[ ]:


def big_m_find_2(xi_1,xi_2,b_1,b_2):
    
    m = Model()
    
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
    y = m.addVar(lb=0,name="y")
    
    
    m.setObjective(b_1 - sum(xi_1[j]*x[j] for j in range(x_random_size)) , GRB.MAXIMIZE)
    
    # Add objective
#     m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
    m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()

    m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)
    m.update()
    
    
    m.params.OutputFlag=0
    
    m.optimize()
    
    return m.objVal
    


# In[20]:


start=time.time()
Big_M_strength_only = []
for j in range(random_size):
    eta_value=[]
    for i in range(random_size):
        if sum(xi[i]) <= b[i]:
            eta_value.append(max(sum(xi[j]) - b[j],0))
#             print("Special:",i)
        else:  
            division_sort = xi[i][np.argsort(np.divide(xi[j],xi[i]))][::-1]

            division_index = np.argsort(np.divide(xi[j],xi[i]))[x_random_size-np.where(division_sort.cumsum() > b[i])[0][0]:x_random_size]

            eta_value.append(max(sum(xi[j][division_index]) - b[j],0)+theta/epsilon)
#             print(i)
    Big_M_strength_only.append(theta/epsilon+1*eta_value[np.argsort(eta_value)[random_size-math.floor(random_size*epsilon)+1]])
modeltime_big_m_stengthen_only = time.time() - start
Big_M_coefficient_VaR = Big_M_strength_only


# In[ ]:


def big_m_find_strengthen_1(xi_1,b_1):
    
    m = Model()
    
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
    z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    
    y = m.addVar(lb=0,name="y")
    
    m.setObjective(sum(xi_1[j]*x[j] for j in range(x_random_size)) - b_1, GRB.MAXIMIZE)
    
    # Add objective 
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()
    m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) >= b[i]   for i in cut_z_indice)
    m.update() 

#     m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)
    
    m.update()
    
    m.addConstrs( theta/epsilon*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in stengthen_difference)
    m.update()
    m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon)))
    m.update()
    
    m.params.OutputFlag=0
    
    m.optimize()
    
    return m.objVal
    
    
    
    
    
    
    
    
    


# In[ ]:


start=time.time()
Big_M_updated_further_1 = [0]*random_size
for j in stengthen_difference:
    print('current j is', j)
 
    Big_M_updated_further_1[j]=max(big_m_find_strengthen_1(xi[j],b[j]),0)
    
modeltime_big_m_stengthen_further  = time.time() - start


# In[ ]:


Big_M_updated_further_1


# In[ ]:





# In[ ]:


# def big_m_find_strengthen_1_regular(xi_1,b_1):
    
#     m = Model()
    
#     x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
#     slack_t = m.addVar(lb=0,name="slack_t")
    
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    
# #     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
    
#     z_regular = m.addVars(random_size,lb=0,ub=1,name="z_regular")
    
        
#     s = m.addVars(random_size,lb=0,name="s")

#     beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")

#     mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")

#     w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")

#     z = m.addVars(random_size,lb=0,ub=1,name="z")

    
#     m.setObjective(sum(xi_1[j]*x[j] for j in range(x_random_size)) - b_1, GRB.MAXIMIZE)
    
#     # Add objective 
# #     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
#     m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
# #     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
#     m.update()
#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) >= b[i]   for i in cut_z_indice)
#     m.update() 

# #     m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)
    
#     m.update()
    
#     m.addConstrs( theta/epsilon*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in stengthen_difference)
#     m.update()
#     m.addConstrs( sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_regular[i])  for i in stengthen_difference)
#     m.update()
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon)))
#     m.update()
#     m.addConstr(z_regular.sum() >= random_size - math.ceil(random_size*(epsilon)))
#     m.update()
    
#     m.params.OutputFlag=0
    

    


#     m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <= beta + b[i]*z[i] + s[i]  for i in range(random_size))
#     m.update()

#     m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <=  b[i]*z[i]  for i in range(random_size))
#     m.update()

#     m.addConstr( theta*slack_t + s.sum()/float(random_size) + epsilon*beta <=0)
#     m.update()

#     m.addConstrs( x[i] == mu[i,j] + w[i,j] for i in range(x_random_size) for j in range(random_size) )
#     m.update()

#     m.addConstrs( w[j,i] <= 1 - z[i]  for j in range(x_random_size) for i in range(random_size) )
#     m.update()
#     m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in range(random_size))
#     m.update()
#     m.addConstr(z.sum() >= math.ceil(random_size*(1-epsilon))+1)
#     m.update()

#     m.addConstr(sum(c[i] * x[i] for i in range(x_random_size))<=value_approximation_bound)

#     m.update()

    
#     m.optimize()
    

    
    
#     return m.objVal
    
    
    
    
    
    
    
    
    


# In[ ]:


# start=time.time()
# Big_M_updated_further_regular_1 = [0]*random_size
# for j in stengthen_difference:
#     print('current j is', j)
 
#     Big_M_updated_further_regular_1[j]=max(big_m_find_strengthen_1_regular(xi[j],b[j]),0)
    
# modeltime_big_m_stengthen_further_regular_1  = time.time() - start


# In[ ]:


# Big_M_updated_further_regular_1 


# In[ ]:





# In[ ]:


# m = Model()

# x = m.addVars(x_random_size,lb=0,ub=1,name="x")

# slack_t = m.addVar(lb=0,name="slack_t")

# z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")

# #     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")

# z_regular = m.addVars(random_size,lb=0,ub=1,name="z_regular")


# s = m.addVars(random_size,lb=0,name="s")

# beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")

# mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")

# w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")

# z = m.addVars(random_size,lb=0,ub=1,name="z")


# m.setObjective(sum(xi[0][j]*x[j] for j in range(x_random_size)) - b[0], GRB.MAXIMIZE)

# # Add objective 
# #     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
# m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
# #     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
# m.update()
# m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) >= b[i]   for i in cut_z_indice)
# m.update() 

# #     m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)

# m.update()

# m.addConstrs( theta/epsilon*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in stengthen_difference)
# m.update()
# # m.addConstrs( sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_regular[i])  for i in stengthen_difference)
# m.update()
# m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon)))
# m.update()
# # m.addConstr(z_regular.sum() >= random_size - math.ceil(random_size*(epsilon)))
# m.update()

# m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <= beta + b[i]*z[i] + s[i]  for i in range(random_size))
# m.update()

# m.addConstrs( sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <=  b[i]*z[i]  for i in range(random_size))
# m.update()

# m.addConstr( theta*slack_t + s.sum()/float(random_size) + epsilon*beta <=0)
# m.update()

# m.addConstrs( x[i] == mu[i,j] + w[i,j] for i in range(x_random_size) for j in range(random_size) )
# m.update()

# m.addConstrs( w[j,i] <= 1 - z[i]  for j in range(x_random_size) for i in range(random_size) )
# m.update()
# m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in range(random_size))
# m.update()
# m.addConstr(z.sum() >= math.ceil(random_size*(1-epsilon))+1)
# m.update()

# m.addConstr(sum(c[i] * x[i] for i in range(x_random_size))<=value_approximation_bound)

# m.update()
# m.params.OutputFlag=0

# m.optimize()

# z_solution_store_update = m.getAttr('x', z)
# mu_solution_store_update = m.getAttr('x', mu)
# w_solution_store_update = m.getAttr('x', w)
# x_solution_store_update = m.getAttr('x', x)

# print('current obj is', m.ObjVal)


# for i in range(0,10):       

#     m.setObjective(sum(xi[i+1][j]*x[j] for j in range(x_random_size)) - b[i+1], GRB.MAXIMIZE)
    
#     print('current i is', i)
#     z.start = z_solution_store_update 
#     mu.start = mu_solution_store_update
#     w.start = w_solution_store_update
#     x.start = x_solution_store_update

#     m.params.OutputFlag=0 

#     m.optimize()

#     print('Updated obj is', m.ObjVal)


# In[ ]:


# Big_M_updated_1


# In[ ]:





# In[ ]:


def big_m_find_strengthen_2(xi_1,b_1):
    
    m = Model()
    
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
    z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    
    y = m.addVar(lb=0,name="y")
    
    m.setObjective( b_1-sum(xi_1[j]*x[j] for j in range(x_random_size)), GRB.MAXIMIZE)
    
    # Add objective 
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()
    m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) >= b[i]   for i in cut_z_indice)
    m.update() 

#     m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)
    
    m.update()
    
    m.addConstrs( theta/epsilon*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in stengthen_difference)
    m.update()
    m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon)))
    m.update()
    
    m.params.OutputFlag=0
    
    m.optimize()
    
    return m.objVal
    
    
    
    
    
    
    
    
    


# In[ ]:


start=time.time()
Big_M_updated_further_2 = [0]*random_size
for j in stengthen_difference:
    print('current j is', j)
 
    Big_M_updated_further_2[j]=max(big_m_find_strengthen_2(xi[j],b[j]),0)
    
modeltime_big_m_stengthen_further_2  = time.time() - start


# In[ ]:


Big_M_updated_further_2 


# In[ ]:


modeltime_big_m_stengthen_update = modeltime_big_m_stengthen_further + modeltime_big_m_stengthen_further_2


# In[ ]:


# start=time.time()
# Big_M_updated_2 = [0]*random_size
# for j in stengthen_difference:
#     Big_M_updated_2[j]=b[j]
    
    


# In[ ]:


Big_M_updated = Big_M_coefficient


# In[ ]:


# Big_M_updated_1 = Big_M_coefficient
# Big_M_updated_2 = Big_M_coefficient


# In[ ]:


# zero_positions_1 = [index for index, num in enumerate(Big_M_updated_1) if num == 0]


# In[ ]:


# zero_positions_2 = [index for index, num in enumerate(Big_M_updated_2) if num == 0]


# In[ ]:


# # Big-M Method
# #Big_M = [None] * random_size

# # Big_M = [30]*random_size

# def Model3_orginial_1():
#     print ("Begin to solve Big M")
    
#     m = Model()
#     m.setParam('Seed', 2)
#     m.setParam('TimeLimit', 4*60*60)
#     # Create variables
#     m.update()
#     # s1 = m.addVars(random_size,name="s1")
#     x = m.addVars(x_random_size,lb=0,ub=1,name="x")
#     z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
# #     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
#     y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
#     s = m.addVars(random_size,lb=0,name="s")
#     gamma = m.addVars(1,lb=0,name="gamma")
#     lambda_1 = m.addVars(1,lb=0,name="lambda_1")
# #     z_var_continuous = m.addVars(stengthen_difference,name="z")
   
#     m.update()
    

#     m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)

#     m.update()
    
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound_second_dual)
#     m.update()
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= value_alsox_sharp)
#     m.update()
    
#     m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated[i]*(1-z[i]) for i in stengthen_difference)
#     m.update()
#     m.addConstrs(s[i] <=  Big_M_updated[i]*z[i] for i in stengthen_difference)
#     m.update()
#     m.addConstrs(s[i]==0 for i in cut_z_indice)
#     m.update()
#     m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))

    
#     m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
#     m.update()

#     m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_appro[i])  for i in stengthen_difference)
#     m.update()

#     m.update()
#     m.addConstrs(y[i]+gamma[0] <= s[i] for i in stengthen_difference)
    
#     m.update()
#     m.addConstrs(y[i]+gamma[0] <= 0 for i in cut_z_indice)
    
#     m.update()
#     m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
#     m.update()
    
#     m.addConstrs(z[i] == 0 for i in cut_z_indice)
#     m.update()

#     m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
#     m.update()
    
#     m.addConstr(z.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
#     m.update()
   
#     for i in range(x_random_size):

#         x[i].start = solution_alsox_sharp[i]

#     m.optimize()
    
#     # Store solutions
#     #pp = m.getAttr('x', z)
#     z_VaR = m.getAttr('x', z)
#     x_VaR = m.getAttr('x', x)
    
#     #return pp,ii,kk
    
#     return m.objVal,z_VaR,x_VaR
    


# In[ ]:


# start=time.time()
# Big_M_solution_1 = Model3_orginial_1()
# z_VaR_1 = Big_M_solution[1]
# x_VaR_1 = Big_M_solution[2]
# # value_z = Model3_cut()[0]
# # kk=Model3_cut()[1]
# # obj_cut = 0
# # for i in range(x_random_size):
# #     obj_cut += c[i]*kk[i]
# modeltime_Model3_var_1 = time.time() - start
# # print(f'The value of optimal value t is exactly {obj_cut_stengthen_updated:.5f}.')
# print(f'The value of running time is approximately {modeltime_Model3_var:.3f} s.')
# # print(f'The value of total running time is approximately {modeltime_cut_stengthen_updated+modeltime_big_m_stengthen+modeltime_fixing:.3f} s.')


# In[ ]:


# Big-M Method
#Big_M = [None] * random_size

# Big_M = [30]*random_size

def Model3_orginial_improve_lower_bound_combine_relax_drccp():
    print ("Begin to solve Big M")
    
    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 20*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
#     z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
    z = m.addVars(random_size,lb=0,ub=1,name="z")
    z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
    s = m.addVars(random_size,lb=0,name="s")
    gamma = m.addVars(1,lb=0,name="gamma")
    lambda_1 = m.addVars(1,lb=0,name="lambda_1")
#     z_var_continuous = m.addVars(stengthen_difference,name="z")
   
    m.update()
    

    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)

    m.update()
    
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= v_q)
    m.update()
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= value_alsox_sharp)
    m.update()
    
    m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated[i]*(1-z[i]) for i in stengthen_difference)
    m.update()
    m.addConstrs(s[i] <=  Big_M_updated[i]*z[i] for i in stengthen_difference)
    m.update()
    m.addConstrs(s[i] == 0 for i in cut_z_indice)
    m.update()
#     m.addConstr(lambda_1[0]>=sum(x[i] for i in range(x_random_size)))
    m.update()
    m.addConstrs(lambda_1[0]>= x[i] for i in range(x_random_size))
    m.update()

    
    m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
    m.update()

    m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in stengthen_difference)
    m.update()

    m.update()
    m.addConstrs(y[i]+gamma[0] <= s[i] for i in stengthen_difference)
    
    m.update()
    m.addConstrs(y[i]+gamma[0] <= 0 for i in cut_z_indice)
    
    m.update()
    m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
    m.update()
    
    m.addConstrs(z[i] == 0 for i in cut_z_indice)
    m.update()

    m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
   
    m.update()
    
#     m.addConstrs(z[stengthen_difference_update[cut_third_three[i][0]]]+z[stengthen_difference_update[cut_third_three[i][1]]]+z[stengthen_difference_update[cut_third_three[i][2]]] <= 2 for i in range(len(cut_third)))
   
    m.update()
    
    
    m.addConstr(z.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    
#     m.addConstrs(z[i]>=z_var_appro_solution[i] for i in range(random_size))
    m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    m.update()
   
#     initial_x = big_m_improve_lower[2]
#     for i in range(x_random_size):

#         x[i].start = initial_x[i]
        
#     initial_z = big_m_improve_lower[1]
#     for i in range(random_size):

#         z[i].start = initial_z[i]


    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_VaR = m.getAttr('x', z_appro)
    x_VaR = m.getAttr('x', x)
    
    lower_bound = m.getAttr('ObjBound')
    #return pp,ii,kk
    
    return m.objVal,z_VaR,x_VaR,lower_bound
    


# In[ ]:


start=time.time() 
lower_var_solution_combine = Model3_orginial_improve_lower_bound_combine_relax_drccp()
modeltime_lower_improve_var = time.time() - start


# In[ ]:


z_VaR_initial = lower_var_solution_combine[1]
x_VaR_initial = lower_var_solution_combine[2]


# In[ ]:


lower_var_solution_combine_value = lower_var_solution_combine[3]


# In[ ]:


lower_var_solution_combine_value


# In[ ]:





# In[ ]:


# satisfy_indice


# In[ ]:


modeltime_var_lower  = modeltime_big_m_stengthen_only + modeltime_lower_improve_var


# In[ ]:


lower_bound = max(lower_var_solution_combine_value,lower_bound_second_dual)


# In[ ]:


z_var_appro_solution_new = lower_var_solution_combine[1]


# In[ ]:


lower_bound


# In[ ]:


value_alsox_sharp


# In[ ]:


value_approximation_bound


# In[ ]:


# DRCC set
def Model_checking_gamma(V_x):
    
    m = Model()
    #m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()

#     x = m.addVars(x_random_size,vtype=GRB.BINARY,name="x")
#     x = m.addVars(x_random_size,name="x")

    s = m.addVars(random_size,lb=0,name="s")

    gamma = m.addVars(1,lb=-GRB.INFINITY,ub=0,name="gamma")
    lambda_1 = m.addVar(lb=0,name="lambda_1")
    m.update()
    

    # Add objective    
    m.setObjective(0, GRB.MINIMIZE)
    m.update()
    
    # Add constraints
   
    m.addConstrs( -max(b[i]-sum(xi[i][j]*V_x[j] for j in range(x_random_size)),0) - gamma[0]  <=   s[i]  for i in range(random_size))
    m.update()
#     m.addConstrs(s[j] >= gamma for j in range(random_size))

    m.addConstr(theta*lambda_1 + s.sum()/float(random_size) + epsilon*gamma[0] <= 0)
    m.update()
    
#     m.addConstr(sum(V_x[j] for j in range(x_random_size))<=lambda_1)
    m.update()
     
    m.addConstrs(lambda_1 >= V_x[j] for j in range(x_random_size))
    
    # Solve the problem
    m.update()
    m.params.OutputFlag=0
        
    m.optimize()
    result = 1
    if  m.status == 4 or m.status == 3 or m.status ==12:
        result = 0
    gamma_Var = m.getAttr('x', gamma)
   
    return result,gamma_Var


# In[ ]:


gamma_alsox_sharp = abs(Model_checking_gamma(solution_alsox_sharp)[1].values()[0])


# In[ ]:


# new_quantile_index_alsox_sharp =[]
# for i in range(random_size):
#     new_quantile_index_alsox_sharp.append(sum(xi[i][j]*solution_alsox_sharp[j]  for j in range(x_random_size)))
# subset_indice_new_alsox_sharp = np.argsort(new_quantile_index_alsox_sharp)
# z_information = [0] * 1000
# for index in subset_indice_new_alsox_sharp[:random_size-k+1]:
#     z_information[index] = 1


# In[ ]:





# In[ ]:


# # AM algorithm 
# v_upper_best = value_alsox_sharp
# v_lower_best = lower_bound
# current_v = (v_upper_best+v_lower_best)/float(2.0)
# delta_v = v_upper_best - v_lower_best
# delta_1 = 1e-1

# m = Model()
# x = m.addVars(x_random_size,lb=0,ub=1,name="x")

# s = m.addVars(random_size,lb=0,name="s")
# gamma = m.addVars(random_size,lb=0,name="gamma")
# beta = m.addVars(1,lb=0,name="beta")
# lambda_1 = m.addVars(1,lb=0,name="lambda_1")

# m.setObjective(theta*lambda_1[0]+epsilon*beta[0]+sum(s[i] for i in range(random_size))/float(random_size), GRB.MINIMIZE)
# m.update()

# m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))
# m.update()

# con1=m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size))<=current_v)
# m.update()

# m.addConstrs(s[i]>=gamma[i]-beta[0] for i in range(random_size))
# m.update()

# m.addConstrs(gamma[i]>=(sum(xi[i][j]*x[j]  for j in range(x_random_size)) - b[i])*z_information[i]  for i in range(random_size))
# m.update()

# m.params.OutputFlag=0

# m.optimize()


# s_solution_store_am = m.getAttr('x', s)
# gamma_solution_store_am = m.getAttr('x', gamma)
# beta_solution_store_am = m.getAttr('x', beta)
# x_solution_store_am = m.getAttr('x', x)

# while delta_v >= delta_1: 
#     iteration = 0
#     z_information_update =  [0] * 1000
#     z_difference = sum(abs(z_information[i]-z_information_update[i]) for i in range(random_size))
    
#     while iteration<=20 or z_difference>=1:
    
#         m.remove(m.getConstrs()[-(random_size)-1:-1])


#         new_quantile_index_new_solution =[]
#         for i in range(random_size):
#             new_quantile_index_new_solution.append(sum(xi[i][j]*x_solution_store_am[j]  for j in range(x_random_size)))
#         subset_indice_new_solution = np.argsort(new_quantile_index_new_solution)
#         z_information_update = [0] * 1000
#         for index in subset_indice_new_solution[:random_size-k+1]:
#             z_information_update[index] = 1

#         z_difference = sum(abs(z_information[i]-z_information_update[i]) for i in range(random_size))

#         z_information = z_information_update

#         m.addConstrs(gamma[i]>=(sum(xi[i][j]*x[j]  for j in range(x_random_size)) - b[i])*z_information[i]  for i in range(random_size))
#         m.update()

#         m.optimize()

#         s_solution_store_am = m.getAttr('x', s)
#         gamma_solution_store_am = m.getAttr('x', gamma)
#         beta_solution_store_am = m.getAttr('x', beta)
#         x_solution_store_am = m.getAttr('x', x)


#         iteration = iteration + 1
    

#     checking = Model_checking(x_solution_store_am)
#     if checking == 1:
#         v_upper_best = current_v
#     else:
#         v_lower_best = current_v

#     current_v = (v_upper_best+v_lower_best)/float(2.0)
#     delta_v = v_upper_best - v_lower_best
#     con1.RHS = current_v
        
# con1.RHS = v_upper_best
# m.optimize()

# current_x_am = m.getAttr('x', x)





# In[ ]:





# In[ ]:





# In[ ]:


# m = Model()
# x = m.addVars(x_random_size,lb=0,ub=1,name="x")

# y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
# s = m.addVars(random_size,lb=0,name="s")
# gamma = m.addVars(1,lb=0,name="gamma")
# lambda_1 = m.addVars(1,lb=0,name="lambda_1")

# m.setObjective(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size), GRB.MINIMIZE)
# m.update()

# m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))
# m.update()

# con1=m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size))<=-16660)
# m.update()

# m.addConstrs(y[i]+gamma[0] <= s[i] for i in range(random_size))

# m.update()
# m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated_1[i]*(1-z_information[i]) for i in range(random_size))
# m.update()
# m.addConstrs(s[i] <=  Big_M_updated_2[i]*z_information[i] for i in range(random_size))
# m.update()


# m.params.OutputFlag=0

# m.optimize()


# s_solution_store_am = m.getAttr('x', s)
# gamma_solution_store_am = m.getAttr('x', gamma)
# x_solution_store_am = m.getAttr('x', x)


# iteration = 0
# z_information_update =  [0] * 1000
# z_difference = sum(abs(z_information[i]-z_information_update[i]) for i in range(random_size))

# while z_difference>=1:

#     m.remove(m.getConstrs()[-2*(random_size)-1:-1])


#     new_quantile_index_new_solution =[]
#     for i in range(random_size):
#         new_quantile_index_new_solution.append(sum(xi[i][j]*x_solution_store_am[j]  for j in range(x_random_size)))
#     subset_indice_new_solution = np.argsort(new_quantile_index_new_solution)
#     z_information_update = [0] * 1000
#     for index in subset_indice_new_solution[:random_size-k+1]:
#         z_information_update[index] = 1

#     z_difference = sum(abs(z_information[i]-z_information_update[i]) for i in range(random_size))
#     print(z_difference)

#     z_information = z_information_update

#     m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated_1[i]*(1-z_information[i]) for i in range(random_size))
#     m.update()
#     m.addConstrs(s[i] <=  Big_M_updated_2[i]*z_information[i] for i in range(random_size))
#     m.update()

#     m.optimize()

#     s_solution_store_am = m.getAttr('x', s)
#     gamma_solution_store_am = m.getAttr('x', gamma)
#     x_solution_store_am = m.getAttr('x', x)
#     lambda_solution_store_am = m.getAttr('x', lambda_1)
    


#     iteration = iteration + 1


# In[ ]:





# In[ ]:


# Model_checking(x_solution_store_am)


# In[ ]:


# Big-M Method
#Big_M = [None] * random_size

# Big_M = [30]*random_size

def Model3_orginial_improve_lower_bound():
    print ("Begin to solve Big M")
    
    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 10*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
#     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
    s = m.addVars(random_size,lb=0,name="s")
    gamma = m.addVars(1,lb=0,name="gamma")
    lambda_1 = m.addVars(1,lb=0,name="lambda_1")
#     z_var_continuous = m.addVars(stengthen_difference,name="z")
   
    m.update()
    

    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)

    m.update()
    
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound)
    m.update()
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= value_alsox_sharp)
    m.update()
    
    m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated_further_1[i]*(1-z[i]) for i in stengthen_difference)
    m.update()
    m.addConstrs(s[i] <=  Big_M_updated_further_2[i]*z[i] for i in stengthen_difference)
    m.update()
    m.addConstrs(s[i]==0 for i in cut_z_indice)
    m.update()
#     m.addConstr(lambda_1[0]>=sum(x[i] for i in range(x_random_size)))
    m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))

    
#     m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
    m.update()

#     m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_appro[i])  for i in stengthen_difference)
    m.update()

    m.update()
    m.addConstrs(y[i]+gamma[0] <= s[i] for i in stengthen_difference)
    
    m.update()
    m.addConstrs(y[i]+gamma[0] <= 0 for i in cut_z_indice)
    
    m.update()
    m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
    m.update()
    
    m.addConstrs(z[i] == 0 for i in cut_z_indice)
    m.update()
    

    m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
    m.update()
    
#     m.addConstrs(z[stengthen_difference_update[cut_third_three[i][0]]]+z[stengthen_difference_update[cut_third_three[i][1]]]+z[stengthen_difference_update[cut_third_three[i][2]]] <= 2 for i in range(len(cut_third)))
   
#     m.update()
    
    m.addConstr(z.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    m.update()
    
#     m.addConstrs(z[stengthen_difference_update_ascending[satisfy_indice[i][0]]] == 1 for i in range(len(satisfy_indice)))
#     m.update()
    
    m.addConstrs(z[i]>=z_var_appro_solution_new[i] for i in range(random_size))
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    m.update()
   
#     for i in range(x_random_size):

#         x[i].start = solution_alsox_sharp[i]
        
#     for i in range(random_size):

#         z[i].start = z_var_appro_solution[i]


    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_VaR = m.getAttr('x', z)
    x_VaR = m.getAttr('x', x)
    gamma_Var = m.getAttr('x', gamma)
    y_Var = m.getAttr('x', y)
    s_Var = m.getAttr('x', s)
    
    
    
    #return pp,ii,kk
    
    return m.objVal,z_VaR,x_VaR,gamma_Var,y_Var,s_Var
    


# In[ ]:


big_m_improve_lower = Model3_orginial_improve_lower_bound()


# In[ ]:


Model_checking(big_m_improve_lower[2])


# In[ ]:


initial_x = big_m_improve_lower[2]
initial_z = big_m_improve_lower[1]
initial_gamma = big_m_improve_lower[3]
initial_y = big_m_improve_lower[4]
initial_s = big_m_improve_lower[5]


# In[ ]:





# In[ ]:


start=time.time() 
m = Model()
m.setParam('Seed', 2)
m.setParam('TimeLimit', 20*60)
# Create variables
m.update()
# s1 = m.addVars(random_size,name="s1")
x = m.addVars(x_random_size,lb=0,ub=1,name="x")
z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
# z = m.addVars(random_size,lb=0,ub=1,name="z")
#     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
s = m.addVars(random_size,lb=0,name="s")
gamma = m.addVars(1,lb=0,name="gamma")
lambda_1 = m.addVars(1,lb=0,name="lambda_1")
#     z_var_continuous = m.addVars(stengthen_difference,name="z")

m.update()


m.setObjective(gamma[0], GRB.MINIMIZE)

m.update()

m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound)
m.update()
m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= big_m_improve_lower[0])
m.update()

m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated_further_1[i]*(1-z[i]) for i in stengthen_difference)
m.update()
m.addConstrs(s[i] <=  Big_M_updated_further_2[i]*z[i] for i in stengthen_difference)
m.update()
m.addConstrs(s[i]==0 for i in cut_z_indice)
m.update()
#     m.addConstr(lambda_1[0]>=sum(x[i] for i in range(x_random_size)))
m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))


#     m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
m.update()

#     m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_appro[i])  for i in stengthen_difference)
m.update()

m.update()
m.addConstrs(y[i]+gamma[0] <= s[i] for i in stengthen_difference)

m.update()
m.addConstrs(y[i]+gamma[0] <= 0 for i in cut_z_indice)

m.update()
m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
m.update()

m.addConstrs(z[i] == 0 for i in cut_z_indice)
m.update()


m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
m.update()

#     m.addConstrs(z[stengthen_difference_update[cut_third_three[i][0]]]+z[stengthen_difference_update[cut_third_three[i][1]]]+z[stengthen_difference_update[cut_third_three[i][2]]] <= 2 for i in range(len(cut_third)))

#     m.update()

m.addConstr(z.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
m.update()

#     m.addConstrs(z[stengthen_difference_update_ascending[satisfy_indice[i][0]]] == 1 for i in range(len(satisfy_indice)))
#     m.update()

m.addConstrs(z[i]>=z_var_appro_solution_new[i] for i in range(random_size))
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
m.update()

# for i in range(x_random_size):

#     x[i].start = solution_alsox_sharp[i]

gamma[0].start = gamma_alsox_sharp
#     for i in range(random_size):

#         z[i].start = z_var_appro_solution[i]


m.optimize()

# Store solutions
#pp = m.getAttr('x', z)
z_VaR_relaxed = m.getAttr('x', z)
x_VaR_relaxed = m.getAttr('x', x)
y_Var_relaxed = m.getAttr('x', y)
s_Var_relaxed = m.getAttr('x', s)
gamma_Var_relaxed_lower = m.getAttr('ObjBound')

modeltime_gamma_lower = time.time() - start


# In[ ]:


gamma_Var_relaxed_lower 


# In[ ]:


start=time.time() 
m = Model()
m.setParam('Seed', 2)
m.setParam('TimeLimit', 20*60)
# Create variables
m.update()
# s1 = m.addVars(random_size,name="s1")
x = m.addVars(x_random_size,lb=0,ub=1,name="x")

# Create a dictionary to store the variables
# z = {}

# # Add the variables to the model
# for i in range(0, random_size):
#     if i <= int(0.95*random_size):
#         z[i] = m.addVar(vtype=GRB.BINARY, name=f"z_{i}")
#     else:
#         z[i] = m.addVar(lb=0,ub=1, name=f"z_{i}")
        
# z = m.addVars(random_size,lb=0,ub=1,name="z")

z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
# z = m.addVars(random_size,lb=0,ub=1,name="z")
#     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
s = m.addVars(random_size,lb=0,name="s")
gamma = m.addVars(1,lb=0,name="gamma")
lambda_1 = m.addVars(1,lb=0,name="lambda_1")
#     z_var_continuous = m.addVars(stengthen_difference,name="z")

m.update()


m.setObjective(gamma[0], GRB.MAXIMIZE)

m.update()

m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound)
m.update()
m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= big_m_improve_lower[0])
m.update()

m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated_further_1[i]*(1-z[i]) for i in stengthen_difference)
m.update()
m.addConstrs(s[i] <=  Big_M_updated_further_2[i]*z[i] for i in stengthen_difference)
m.update()
m.addConstrs(s[i]==0 for i in cut_z_indice)
m.update()
#     m.addConstr(lambda_1[0]>=sum(x[i] for i in range(x_random_size)))
m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))


#     m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
m.update()

#     m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_appro[i])  for i in stengthen_difference)
m.update()

m.update()
m.addConstrs(y[i]+gamma[0] <= s[i] for i in stengthen_difference)

m.update()
m.addConstrs(y[i]+gamma[0] <= 0 for i in cut_z_indice)

m.update()
m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
m.update()

m.addConstrs(z[i] == 0 for i in cut_z_indice)
m.update()


m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
m.update()

#     m.addConstrs(z[stengthen_difference_update[cut_third_three[i][0]]]+z[stengthen_difference_update[cut_third_three[i][1]]]+z[stengthen_difference_update[cut_third_three[i][2]]] <= 2 for i in range(len(cut_third)))

#     m.update()

m.addConstr(sum(z[i] for i in range(random_size)) >= random_size - math.ceil(random_size*(epsilon))+1)
m.update()

#     m.addConstrs(z[stengthen_difference_update_ascending[satisfy_indice[i][0]]] == 1 for i in range(len(satisfy_indice)))
#     m.update()

m.addConstrs(z[i]>=z_var_appro_solution_new[i] for i in range(random_size))
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
m.update()

for i in range(x_random_size):

    x[i].start = solution_alsox_sharp[i]

#     for i in range(random_size):

#         z[i].start = z_var_appro_solution[i]
gamma[0].start = gamma_alsox_sharp

m.optimize()

# Store solutions
#pp = m.getAttr('x', z)
z_VaR_relaxed = m.getAttr('x', z)
x_VaR_relaxed = m.getAttr('x', x)
gamma_Var_relaxed_upper = m.getAttr('ObjBound')
y_Var_relaxed = m.getAttr('x', y)
s_Var_relaxed = m.getAttr('x', s)

modeltime_gamma_upper = time.time() - start


# In[ ]:


gamma_Var_relaxed_upper


# In[ ]:





# In[ ]:


# Big-M Method
#Big_M = [None] * random_size

# Big_M = [30]*random_size

def Model3_orginial_improve_lower_bound_h():
    print ("Begin to solve Big M")
    
    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 8*60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
#     z_appro_1 = m.addVars(random_size,lb=0,ub=1,name="z_appro_1")
    
#     z_appro_1 = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro_1")
#     z_appro_2 = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro_2")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
    s = m.addVars(random_size,lb=0,name="s")
    gamma = m.addVars(1,lb=0,name="gamma")
    lambda_1 = m.addVars(1,lb=0,name="lambda_1")
    t = m.addVars(1,lb=0,name="t")
#     z_var_continuous = m.addVars(stengthen_difference,name="z")
   
    m.update()
    

    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)

    m.update()
    
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound)
  
    m.update()
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= value_alsox_sharp)
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= big_m_improve_lower[0])
#     m.update()
    
    m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated_further_1[i]*(1-z[i]) for i in stengthen_difference)
    m.update()
    m.addConstrs(s[i] <=  Big_M_updated_further_2[i]*z[i] for i in stengthen_difference)
    m.update()
    m.addConstrs(s[i] == 0 for i in cut_z_indice)
    m.update()
#     m.addConstr(lambda_1[0]>=sum(x[i] for i in range(x_random_size)))
    m.update()
    m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))
    m.update()  

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) >= b[i]   for i in cut_z_indice)
    m.update()    
#     m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated_1[i]*(1-z[i])  for i in stengthen_difference)
#     m.update()

#     m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
#     m.update()

    m.update()
    m.addConstrs(y[i]+gamma[0] <= s[i] for i in stengthen_difference)
    m.update()
    
#     m.addConstrs(y[i]+gamma[0] >= 0 for i in stengthen_difference)
    
#     m.update()
    m.addConstrs(y[i]+gamma[0] == 0 for i in cut_z_indice)
    
    m.update()
    
    m.addConstr(gamma[0] <=gamma_Var_relaxed_upper)
    
    m.update()
    
    m.addConstr(gamma[0] >=gamma_Var_relaxed_lower)
    m.update()
    
    
#     m.addConstr(gamma[0]<=800)
#     m.addConstrs(gamma[0]>= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  - Big_M_updated[i]*(1-z_appro_1[i]) for i in range(random_size))
#     m.update()
# # #     m.addConstrs(gamma[0]<= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated[i]*(1-z_appro_2[i]) for i in range(random_size))
# # #     m.update()
#     m.addConstr(sum(z_appro_1[i] for i in range(random_size))==k)
#     m.update()
#     m.addConstr(sum(z_appro_2[i] for i in range(random_size))==random_size-k+1)
#     m.update()
#     m.addConstrs(z_appro_1[i]+z_appro_2[i]>=1 for i in range(random_size))
#     m.update()
#     m.addConstrs(z_appro_2[i]<=z[i] for i in range(random_size))
#     m.update()
    
    
    
#     m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] <=0)
#     m.update()

    m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
#     m.addConstr(theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
   
    m.update()
    
    m.addConstrs(z[i] == 0 for i in cut_z_indice)
    m.update()

    m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
    m.update()
    
#     m.addConstrs(z[stengthen_difference_update[cut_third_three[i][0]]]+z[stengthen_difference_update[cut_third_three[i][1]]]+z[stengthen_difference_update[cut_third_three[i][2]]] <= 2 for i in range(len(cut_third)))
   
    m.update()
    
#     m.addConstrs(z[i] == 0 for i in zero_positions_1) 
    
#     m.addConstrs(z[violation_indice_new_difference_update[violation_indice[i][0]]]+z[violation_indice_new_difference_update[violation_indice[i][1]]] <= 1 for i in range(len(violation_indice)))
#     m.update()
    
#     m.addConstrs(z[stengthen_difference_update_ascending[satisfy_indice[i][0]]] == 1 for i in range(len(satisfy_indice)))
#     m.update()
    
    m.addConstr(z.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    
#     m.addConstrs(z[i]>=z_var_appro_solution[i] for i in range(random_size))
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    m.update()
   
#     initial_x = big_m_improve_lower[2]
#     for i in range(x_random_size):
# #         x[i].start = solution_alsox_sharp[i]

#         x[i].start = initial_x[i]
    
    for i in range(random_size):
        y[i].start = initial_y[i]
        s[i].start = initial_s[i]
#         z[i].start = initial_z[i]
        
    
    gamma[0].start = initial_gamma[0]
        
#     initial_z = big_m_improve_lower[1]
#     for i in range(random_size):

#         z[i].start = initial_z[i]
#     m.setParam('MIPFocus', 2)
#     m.setParam('Cuts', 2)
#     m.setParam('Method', 3)
#     m.setParam('Cutoff', big_m_improve_lower[0])
    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_VaR = m.getAttr('x', z)
    x_VaR = m.getAttr('x', x)
    gamma_Var_optimal = m.getAttr('x', gamma)
    
    #return pp,ii,kk
    
    return m.objVal,z_VaR,x_VaR,gamma_Var_optimal
    


# In[ ]:


start=time.time() 
final_big_m = Model3_orginial_improve_lower_bound_h()
modeltime_final_big_m = time.time() - start
# initial var solution 
# add cuts, inital alsox#


# In[ ]:





# In[ ]:


final_big_m[3]


# In[ ]:


modeltime_big_m_stengthen_update


# In[ ]:


modeltime_gamma_stengthen = modeltime_gamma_lower + modeltime_gamma_upper


# In[ ]:


model_preprocess_time = modeltime_cvar+modeltime_quantile+modeltime_alsox_sharp

model_our_method = model_preprocess_time + modeltime_second_dual_bound + modeltime_var_lower + modeltime_fixing + modeltime_gamma_stengthen+ modeltime_big_m_stengthen_update + modeltime_final_big_m


# In[ ]:


print(f'The value of model_preprocess_time is exactly {model_preprocess_time:.5f}.')
print(f'The value of modeltime_second_dual_bound is exactly {modeltime_second_dual_bound:.5f}.')
print(f'The value of modeltime_var_lower is exactly {modeltime_var_lower:.5f}.')
print(f'The value of modeltime big strengthen is exactly {modeltime_big_m_stengthen_update:.5f}.')
print(f'The value of modeltime gamma strengthen is exactly {modeltime_gamma_stengthen:.5f}.')
print(f'The value of modeltime_fixing is exactly {modeltime_fixing:.5f}.')
print(f'The value of modeltime_final_big_m is exactly {modeltime_final_big_m:.5f}.')
print(f'The value of model_our_method is exactly {model_our_method:.5f}.')
print(f'The number of cut z_i is exactly {len(cut_z_indice):.0f}.')
print(f'The number of cut z_i+z_j is exactly {len(cut_third):.0f}.')


# In[ ]:





# In[ ]:





# In[ ]:


# def find_maximum_elements(arr1, arr2):
#     # Make sure the arrays have the same length
#     if len(arr1) != len(arr2):
#         raise ValueError("Arrays must have the same length")

#     # Create a new array to store the maximum values
#     max_elements = []

#     # Iterate through the arrays and find the maximum for each element
#     for i in range(len(arr1)):
#         max_value = max(arr1[i], arr2[i])
#         max_elements.append(max_value)

#     return max_elements


# In[ ]:


# big_m_update_str = find_maximum_elements(Big_M_updated_further_2,Big_M_updated_further_1)


# In[21]:


def big_m_find_strengthen_only_1(xi_1,b_1):
    
    m = Model()
    
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
    z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    
    y = m.addVar(lb=0,name="y")
    
    m.setObjective(sum(xi_1[j]*x[j] for j in range(x_random_size)) - b_1, GRB.MAXIMIZE)
    
    # Add objective 
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()
#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) >= b[i]   for i in cut_z_indice)
#     m.update() 

#     m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)
    
    m.update()
    
    m.addConstrs( theta/epsilon*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in range(random_size))
    m.update()
    m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon)))
    m.update()
    
    m.params.OutputFlag=0
    
    m.optimize()
    
    return m.objVal
    
    
    
    
    
    
    
    
    


# In[22]:


def big_m_find_strengthen_only_2(xi_1,b_1):
    
    m = Model()
    
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    
    slack_t = m.addVar(lb=0,name="slack_t")
    
    z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    
    y = m.addVar(lb=0,name="y")
    
    m.setObjective( b_1-sum(xi_1[j]*x[j] for j in range(x_random_size)), GRB.MAXIMIZE)
    
    # Add objective 
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
#     m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
    m.update()
#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) >= b[i]   for i in cut_z_indice)
    m.update() 

#     m.addConstr( theta/epsilon*slack_t + sum(xi_2[j]*x[j] for j in range(x_random_size))  <= b_2)
    
    m.update()
    
    m.addConstrs( theta/epsilon*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_coefficient_VaR[i]*(1-z_appro[i])  for i in range(random_size))
    m.update()
    m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon)))
    m.update()
    
    m.params.OutputFlag=0
    
    m.optimize()
    
    return m.objVal
    
    
    
    
    
    
    
    
    


# In[23]:


start=time.time()
Big_M_updated_further_only_1 = [0]*random_size
for j in range(random_size):
    print('current j is', j)
 
    Big_M_updated_further_only_1[j]=max(big_m_find_strengthen_only_1(xi[j],b[j]),0)
    
modeltime_big_m_stengthen_further_only_1  = time.time() - start


# In[24]:


start=time.time()
Big_M_updated_further_only_2 = [0]*random_size
for j in range(random_size):
    print('current j is', j)
 
    Big_M_updated_further_only_2[j]=max(big_m_find_strengthen_only_2(xi[j],b[j]),0)
    
modeltime_big_m_stengthen_further_only_2  = time.time() - start


# In[25]:


model_Big_M_updated_further_only = modeltime_big_m_stengthen_further_only_1 + modeltime_big_m_stengthen_further_only_2


# In[26]:


# Big_M_updated_further_only_11 = [math.ceil(num * 100) / 100 for num in Big_M_updated_further_only_1]
# Big_M_updated_further_only_11 = np.floor(Big_M_updated_further_only_1)


# In[27]:


# Big-M Method
#Big_M = [None] * random_size

# Big_M = [30]*random_size

def Model3_big_m_strengthen():
    print ("Begin to solve Big M")
    
    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 4*60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
#     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
    s = m.addVars(random_size,lb=0,name="s")
    gamma = m.addVars(1,lb=0,name="gamma")
    lambda_1 = m.addVars(1,lb=0,name="lambda_1")
#     z_var_continuous = m.addVars(stengthen_difference,name="z")
   
    m.update()
    

    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)

    m.update()
    
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound_second_dual)
#     m.update()
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= big_m_improve_lower[0])
#     m.update()
    
    m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated_further_only_1[i]*(1-z[i]) for i in range(random_size))
    m.update()
    m.addConstrs(s[i] <= Big_M_updated_further_only_2[i]*z[i] for i in range(random_size))
    m.update()
#     m.addConstrs(s[i] == 0 for i in cut_z_indice)
#     m.update()
    m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))

    
#     m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
    m.update()

#     m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_appro[i])  for i in stengthen_difference)
    m.update()

    m.update()
    m.addConstrs(y[i]+gamma[0] <= s[i] for i in range(random_size))
    

    
    m.update()
    m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
    m.update()
    
#     m.addConstrs(z[i] == 0 for i in cut_z_indice)
    m.update()

#     m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
    m.update()
    
#     m.addConstr(z.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    
#     m.addConstrs(z[i]>=z_var_appro_solution[i] for i in range(random_size))
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    m.update()
    
#     initial_x = big_m_improve_lower[2]
#     for i in range(x_random_size):

#         x[i].start = solution_alsox_sharp[i]
        
#     initial_z = big_m_improve_lower[1]
#     for i in range(random_size):

#         z[i].start = initial_z[i]


    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_VaR = m.getAttr('x', z)
    x_VaR = m.getAttr('x', x)
    
    #return pp,ii,kk
    
    return m.objVal,z_VaR,x_VaR
    


# In[28]:


start=time.time() 
Model3_big_m_strengthen_value = Model3_big_m_strengthen()
modeltime_big_m_strengthen_only  = time.time() - start



# In[ ]:





# In[29]:


print(f'The value of model_Big_M_updated_further_only is exactly {model_Big_M_updated_further_only:.5f}.')
print(f'The value of modeltime_big_m_strengthen_only is exactly {modeltime_big_m_strengthen_only:.5f}.')


# In[30]:


# Big-M Method
#Big_M = [None] * random_size

# Big_M = [30]*random_size

def Model3_Vanilla():
    print ("Begin to solve Big M")
    
    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 4*60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
#     z_appro = m.addVars(random_size,vtype=GRB.BINARY,name="z_appro")
#     z_appro = m.addVars(random_size,lb=0,ub=1,name="z_appro")
    y = m.addVars(random_size,lb=-GRB.INFINITY,ub=0,name="y")
    s = m.addVars(random_size,lb=0,name="s")
    gamma = m.addVars(1,lb=0,name="gamma")
    lambda_1 = m.addVars(1,lb=0,name="lambda_1")
#     z_var_continuous = m.addVars(stengthen_difference,name="z")
   
    m.update()
    

    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)

    m.update()
    
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound_second_dual)
#     m.update()
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= big_m_improve_lower[0])
#     m.update()
    
    m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_updated[i]*(1-z[i]) for i in range(random_size))
    m.update()
    m.addConstrs(s[i] <=  Big_M_updated[i]*z[i] for i in range(random_size))
    m.update()
#     m.addConstrs(s[i] == 0 for i in cut_z_indice)
#     m.update()
    m.addConstrs(lambda_1[0]>=x[i] for i in range(x_random_size))

    
#     m.addConstrs(z_appro[i] <= z[i] for i in range(random_size))

#     m.addConstrs(  sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z[i])  for i in stengthen_difference)
    m.update()

#     m.addConstrs( theta/epsilon*lambda_1[0] + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + Big_M_updated[i]*(1-z_appro[i])  for i in stengthen_difference)
    m.update()

    m.update()
    m.addConstrs(y[i]+gamma[0] <= s[i] for i in range(random_size))
    

    
    m.update()
    m.addConstr(lambda_1[0]*theta - epsilon*gamma[0] - quicksum(y[i] for i in range(random_size))/float(random_size)<=0)
    m.update()
    
#     m.addConstrs(z[i] == 0 for i in cut_z_indice)
    m.update()

#     m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]] <= 1 for i in range(len(cut_third)))
    m.update()
    
    m.addConstr(z.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    
#     m.addConstrs(z[i]>=z_var_appro_solution[i] for i in range(random_size))
#     m.addConstr(z_appro.sum() >= random_size - math.ceil(random_size*(epsilon))+1)
    m.update()
   
#     initial_x = big_m_improve_lower[2]
#     for i in range(x_random_size):

#         x[i].start = solution_alsox_sharp[i]
        
#     initial_z = big_m_improve_lower[1]
#     for i in range(random_size):

#         z[i].start = initial_z[i]


    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_VaR = m.getAttr('x', z)
    x_VaR = m.getAttr('x', x)
    
    #return pp,ii,kk
    
    return m.objVal,z_VaR,x_VaR
    


# In[31]:


start=time.time() 
big_m_Vanilla = Model3_Vanilla()
modeltime_big_m_Vanilla  = time.time() - start



# In[ ]:


print(f'The value of modeltime_big_m_Vanilla is exactly {modeltime_big_m_Vanilla:.5f}.') 


# In[ ]:




