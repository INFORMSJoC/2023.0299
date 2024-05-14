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
theta = 0.20


# In[5]:


filename = '1-7-5-500-5.txt'
f = open("1-7-5-500-5.txt", "r")
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


# In[19]:


Big_M = []
for i in range(random_size):
    
#     Big_M.append(abs(math.ceil(sum(x for x in xi[i][:] if x > 0) -b[i]))) 
    Big_M.append(max(abs(math.ceil(sum(x for x in xi[i][:] if x > 0) -b[i])),b[i])) 
    
Big_M_coefficient = Big_M 


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
    Big_M_strength_only.append(eta_value[np.argsort(eta_value)[random_size-math.floor(random_size*epsilon)+1]])
modeltime_big_m_stengthen_only = time.time() - start
Big_M_coefficient_regular = Big_M_strength_only


# In[21]:


Big_M_updated_further_only_2 = b


# In[22]:


model_Big_M_updated_further_only = modeltime_big_m_stengthen_only


# In[23]:


print(f'The running time of modeltime_big_m_strengthen_only is exactly {model_Big_M_updated_further_only:.5f}.')


# In[24]:


# Big_M_updated_further_only_11 = [math.ceil(num * 100) / 100 for num in Big_M_updated_further_only_1]
# Big_M_updated_further_only_11 = np.floor(Big_M_updated_further_only_1)


# In[25]:


# Big-M Method
#Big_M = [None] * random_size

# Big_M = [30]*random_size

def Model3_big_m_strengthen():
    print ("Begin to solve Big M")
    
    m = Model()
    m.setParam('Seed', 4)
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
    
    m.addConstrs(s[i] <= b[i] - sum(xi[i][j]*x[j]  for j in range(x_random_size))  + Big_M_coefficient_regular[i]*(1-z[i]) for i in range(random_size))
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
    


# In[26]:


# start=time.time() 
# Model3_big_m_strengthen_value = Model3_big_m_strengthen()
# modeltime_big_m_strengthen_only  = time.time() - start



# In[ ]:





# In[27]:


print(f'The value of model_Big_M_updated_further_only is exactly {model_Big_M_updated_further_only:.5f}.')


# In[ ]:




