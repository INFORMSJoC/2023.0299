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
epsilon = 0.10
theta = 0.10


# In[5]:


filename = '1-7-5-1000-5.txt'
f = open("1-7-5-1000-5.txt", "r")
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
    
    


# In[11]:


np.size(c)


# In[12]:


np.size(b)


# In[13]:


np.size(xi)


# In[14]:


np.size(xi[0])


# In[15]:


random_size = np.size(b)
x_random_size = np.size(c)


# In[16]:


k=math.floor(random_size*epsilon)


# In[17]:


c = -c


# In[18]:


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
    s = m.addVars(random_size,lb=-GRB.INFINITY,name="s")
#     s = m.addVars(random_size,lb=-GRB.INFINITY,name="s")
  #  z = m.addVars(random_size,name="z")
    beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")
    slack_t = m.addVar(lb=0,name="slack_t")
    m.update()
        
    
    # Add objective    
    m.setObjective(sum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    # Add constraints
    m.addConstrs(  s[i] >= beta for i in range(random_size))
    m.update()
    m.addConstr( slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
#     m.addConstrs(slack_t >= x[i] for i in range(x_random_size))
    m.update()

   
    m.addConstrs( theta*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + s[i]  for i in range(random_size))
    m.update()
    m.addConstr( s.sum()/random_size - (1-epsilon)*beta <=0)
    m.update()
    # Solve the problem
    m.update()
        
    m.optimize()
    

#     Store solutions
#     ppp = m.getAttr('x', beta)
    CVAR_s = m.getAttr('x', s)
    CVAR_x = m.getAttr('x', x)
    
#     for i in range(random_size):
#      #   g_x_si[i] = g_x + si[i] 
#         c_x = 0
#         for j in range(x_random_size):
#             c_x += a[i,j]*kkk[j]
#         c_x_xi[i] = c_x  
    
#     for i in range(random_size):
#         if c_x_xi[i]-4 <= 0:
#             z[i] = 1
#         else:
#             z[i] = 0
# #     return pp,ii,kk
#     print('Obj: %g' % m.objVal)
    aaaa=m.objVal
    for v in m.getVars():
        print('%s %g' % (v.varName, v.x))
    return m.objVal,CVAR_s,CVAR_x
    
    


# In[19]:


start=time.time()
CVaR=Model4()
CVAR=CVaR[0]
CVAR_s=CVaR[1]
CVAR_x=CVaR[2]
modeltime_cvar= time.time() - start


# In[20]:


# quantitle bounds

def Model_q():

    print ("Begin to solve model quantitle relaxtion")

    m = Model()
    
    m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1, name="x")
    slack_t = m.addVar(lb=0,name="slack_t")
    m.update()
    
    # Set functions
    
    #obj = 0
    g_x_xi = [None] * random_size
    f_x = 0
    g_x = 0
    
    m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
#     m.addConstrs(slack_t >= x[i] for i in range(x_random_size))
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
    
    con1=m.addConstr(theta*slack_t + g_x_xi[0] <= b[0] )
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
        con1=m.addConstr(theta*slack_t + g_x_xi[i]  <= b[i] )
        m.update()
        m.params.OutputFlag=0    
        m.optimize()
        kkk = m.getAttr('x', x)
        obj = 0
        for j in range(x_random_size):
            obj += c[j]*kkk[j]
        result[i]=obj
    return result
        



# In[21]:


start=time.time()
quantile_result=Model_q()
sorted_array=np.sort(quantile_result,-1)
v_q = sorted_array[random_size-k-1]
print(v_q)
modeltime_quantile = time.time() - start


# In[22]:


v_q


# In[23]:


CVAR


# In[24]:


violation = math.floor(random_size*epsilon)


# In[25]:


CVAR


# In[26]:


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
    s = m.addVars(random_size,lb=-GRB.INFINITY,name="s")
  #  z = m.addVars(random_size,name="z")
    beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")
    slack_t = m.addVar(lb=0,name="slack_t")
    m.update()
    

    # Add objective    
    m.setObjective(s.sum()/random_size - (1-epsilon)*beta, GRB.MINIMIZE)
    m.update()
    # Add constraints
    m.addConstrs(  s[i] >= beta for i in range(random_size))
    m.update()
    m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
#     m.addConstrs(slack_t >= x[i] for i in range(x_random_size))
    m.update()

   
    m.addConstrs( theta*slack_t + sum(xi[i][j]*x[j] for j in range(x_random_size)) <= b[i] + s[i]  for i in range(random_size))
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
    
    
    support=0
    for j in current_s:
        if current_s[j]>1e-6:
            support +=1
    if support > violation:
        v_lower_upper = current_bound
#         t1=con1.SARHSUp
    else:
        v_upper_bound = current_bound
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

    
    
    
#     aaaa=m.objVal
    
#     for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
    return v_upper_bound,current_x
    
    


# In[27]:


start=time.time()    
alsox_sharp = Model_alsox_sharp()
modeltime_alsox_sharp = time.time() - start   
value_alsox_sharp = alsox_sharp[0]
solution_alsox_sharp = alsox_sharp[1]


# In[28]:


v_q


# In[29]:


(CVAR-v_q)/abs(v_q)


# In[30]:


(value_alsox_sharp-v_q)/abs(v_q)


# In[31]:


value_alsox_sharp


# In[32]:


CVAR


# In[33]:


value_approximation_bound = value_alsox_sharp
# value_approximation_bound = min(CVAR,alsox_sharp_result[0])
# value_approximation_bound = CVAR


# In[34]:


sorted_array_indice = [i for i, value in enumerate(quantile_result) if value > value_approximation_bound]


# In[35]:


cut_z_indice = sorted_array_indice


# In[36]:


stengthen_difference = [element for element in range(random_size) if element not in cut_z_indice]


# In[37]:


len(cut_z_indice)


# In[ ]:





# In[38]:


# #  lower bound improvment
# start=time.time()   
# delta_1 = 1e-1
# current_bound = v_q
# # v_upper_bound = CVAR
# # current_bound = -20000

# delta_t = 10


# m = Model()
# #m.setParam('TimeLimit', 60*60)
# # Create variables
# m.update()
# # s1 = m.addVars(random_size,name="s1")
# x = m.addVars(x_random_size,lb=0,ub=1,name="x")

# y = m.addVar(lb=-GRB.INFINITY,name="y")

# slack_t = m.addVars(random_size,lb=0,name="slack_t")

# # slack_t = m.addVar(lb=0,name="slack_t")

# mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")

# w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")

# # z = m.addVars(random_size,vtype=GRB.BINARY,name="z")

# z = m.addVars(random_size,lb=0,ub=1,name="z")

# m.update()

# m.setObjective(y, GRB.MINIMIZE)
# # Add objective 

# # m.addConstr( slack_t >= sum(x[i] for i in range(x_random_size)))
# # m.addConstrs( slack_t >= x[i] for i in range(x_random_size))
# m.update()
# # m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
# m.update()
# m.addConstrs(slack_t[i]*slack_t[i] >= sum(mu[j,i]*mu[j,i] for j in range(x_random_size))  for i in range(random_size))
# # m.addConstrs( theta*slack_t + sum(xi[i][j]*mu[j,i] for j in range(x_random_size)) <=  b[i]*z[i]   for i in range(random_size))
# m.update()
# m.addConstrs(theta*slack_t[i]+ sum( xi[i][j]*mu[j,i] for j in range(x_random_size)) <= b[i]*(z[i]) for i in range(random_size)) 
# m.update()

# # Solve the problem
# # m.update()

# m.addConstrs( x[i] == mu[i,j] + w[i,j] for i in range(x_random_size) for j in range(random_size) )
# m.update()
# # m.addConstrs( mu[j,i]<=1-z[j]  for j in range(x_random_size) for i in range(random_size) )
# m.addConstrs( w[j,i] <= 1 - z[i]  for j in range(x_random_size) for i in range(random_size) )
# m.update()
# m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in range(random_size))
# m.update()
# m.addConstr(z.sum() >= random_size - math.ceil(random_size*epsilon))
# m.update()

# constr1 = m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + current_bound * (1-z[j])  for j in  range(random_size))
# m.update()
# constr2 = m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + current_bound * (z[j])  for j in  range(random_size) )
# m.update()

# # m.params.OutputFlag=0
# m.params.BarHomogeneous=1 
# m.optimize()

# # current_s = m.getAttr('x', s)
# current_x = m.getAttr('x', x)

# dual_bound = m.objVal
# print('bound is ', dual_bound)

# ite = 0 

# while delta_t >= delta_1:
    
#     m.remove(m.getConstrs()[-(2*random_size)-1:-1])
# #     m.remove(m.getConstrs()[-2])
#     delta_t = abs(dual_bound - current_bound)
    
#     current_bound = dual_bound
    
# #     if m.status == 2:
# #         current_bound = dual_bound
# #     else:
# #         v_lower_upper = current_bound

# #     current_bound = (v_lower_upper + v_upper_bound)/2.0
    
# #     print('current bound is', current_bound)
#     constr1 = m.addConstrs( y >= sum(c[i] * mu[i,j] for i in range(x_random_size)) + current_bound * (1-z[j])  for j in  range(random_size))
#     m.update()
#     constr2 = m.addConstrs( y >= sum(c[i] * w[i,j] for i in range(x_random_size)) + current_bound * (z[j])  for j in  range(random_size) )
#     m.update()
# #     constr1.RHS = current_bound
# #     constr2.RHS = current_bound
# #     for i in range(random_size):
# #         s[i].start = current_s[i]
#     for i in range(x_random_size):
#         x[i].start = current_x[i]  
#     m.params.OutputFlag=0
#     print('current bound is', current_bound)
#     m.optimize()
#     dual_bound = m.objVal
    
#     ite = ite + 1 
    
# current_bound = dual_bound
# print('with iterations', ite)
# #     Store solutions
# #     ppp = m.getAttr('x', beta)
# # first_dual_bound_s = m.getAttr('x', s)
# # first_dual_bound_x = m.getAttr('x', x)

# # dual_bound = m.objVal

# modeltime_second_dual_bound = time.time() - start   


# In[39]:


# lower_bound_second_dual = current_bound


# In[40]:


# (value_alsox_sharp-lower_bound_second_dual)/abs(lower_bound_second_dual)


# In[ ]:





# In[41]:


def big_m_infty_norm_2(aa,bb,right,alpha):
#     difference = a - b
#     difference_positive_indice = [i for i, value in enumerate(difference) if value > 0]
    
#     # find alpha
#     alpha_upper_bound = sympy.symbols('alpha_bound')
#     eqn1 = sum(pow(a[i] -alpha_upper_bound* b[i],2) for i in difference_positive_indice) - pow((alpha_upper_bound-1)*theta,2)
#     alpha_value=sympy.solve([eqn1], [alpha_upper_bound])


    difference = aa - bb * alpha
    difference[difference<=0]=0
    difference_positive_indice = [i for i, value in enumerate(difference) if value > 0]
    difference_ascending_index = np.argsort(difference)
    S_3 = [element for element in difference_ascending_index if element in difference_positive_indice]
#         print(S_3)
    iteration = 0
    sum_square = 0
    if len(S_3)!=0:
        while iteration <= len(S_3)-1 and sum_square + pow(difference[S_3[iteration]],2) <= pow((alpha-1)*theta,2):
            sum_square = sum_square + pow(difference[S_3[iteration]],2)
            difference[S_3[iteration]] = 0  
            iteration = iteration + 1

        value = alpha*right + sum(difference[S_3[i]] for i in range(len(S_3))) - math.sqrt(iteration) * math.sqrt(pow((alpha-1)*theta,2)-sum_square)
    else:
        value= sum(aa)
    value_return_1 = value

#     value_return_1 = 0
#     value_return = min(value_return_1,value_return_2)
#     print('current value 1 ',value_return_1 )
#     print('current value 2 ',value_return_2 )
#     print(value_return)
    return value_return_1
    



# In[42]:


def big_m_infty_norm_2_less_than_one(aa,bb,right,alpha):
    alpha_new = alpha
    difference_new = aa - bb * alpha_new
    difference_non_negative_indice = [i for i, value in enumerate(difference_new) if value >= 0]
    difference_descending_index = np.argsort(difference_new)[::-1]
    current_index = len(difference_non_negative_indice)
    current_value = (1-alpha_new)*theta*math.sqrt(current_index) + sum(difference_new[i] for i in difference_non_negative_indice)
#         print('current value ',current_value )
    if current_index <= x_random_size-1:
        current_value_update = (1-alpha_new)*theta*math.sqrt(current_index+1) + sum(difference_new[i] for i in difference_descending_index[:current_index+1])
        # print('current value update ',current_value_update )
        while current_value<=current_value_update and current_index<=x_random_size-1:
            current_index = current_index + 1
            current_value = current_value_update
            current_value_update = (1-alpha_new)*theta*math.sqrt(current_index+1) + sum(difference_new[i] for i in difference_descending_index[:current_index+1])
        #     print('current value ',current_value )
        #     print('current value update ',current_value_update )
    else:
        current_value = sum(aa)

    current_value = alpha_new*right + current_value
#         value_new.append(current_value)    
    value_return_2 = current_value
#     print('current value 2',value_return_2 )
    return value_return_2


# In[43]:


def goldenopt(x_l, x_u, tol, max_iter):
    iter_ctr = 0
    while (x_u-x_l)>tol and iter_ctr<max_iter:
        x_r1 = x_u - 0.618*(x_u - x_l)
        x_r2 = x_l + 0.618*(x_u - x_l)
#         print(iter_ctr, x_l, x_r1, x_r2, x_u)
        if x_r1<=1:
            value_x_r1 = big_m_infty_norm_2_less_than_one(aa,bb,right,x_r1)
        else:
            value_x_r1 = big_m_infty_norm_2(aa,bb,right,x_r1)    
        if x_r2<=1:
            value_x_r2 = big_m_infty_norm_2_less_than_one(aa,bb,right,x_r2)
        else:
            value_x_r2 = big_m_infty_norm_2(aa,bb,right,x_r2)    
        if value_x_r1 <= value_x_r2:
            x_u = x_r2
            x_r2 = x_r1
            x_r1 = x_u - 0.618*(x_u - x_l)
        else:
            x_l = x_r1
            x_r1 = x_r2
            x_r2  = x_l + 0.618*(x_u - x_l)
        iter_ctr = iter_ctr + 1
    return (x_r1 + x_r2)/2.0, iter_ctr


# In[44]:


# strengthen
start=time.time()
Big_M_updated = [0]*random_size
for j in stengthen_difference:
    eta_value=[]
    print('current j is', j)
    for i in stengthen_difference:
        if theta*math.sqrt(x_random_size)+sum(xi[i]) <= b[i]:
            eta_value.append(max(sum(xi[j]) - b[j],0))
        else:
            if j==i:
                eta_value.append(max(sum(xi[j]) - b[j],0))
            else:
                aa = xi[j]
                bb = xi[i]
                right = b[i]
#                 eta_value.append(big_m_infty_norm_2(aa,bb,right))
                x_l_init = 0
                x_u_init = 2.0
                tolerance = 5e-2
                max_iter = 50
                x_f, iter_ctr = goldenopt( x_l_init, x_u_init, tolerance, max_iter)
                if x_f<=1:
                    value_final = big_m_infty_norm_2_less_than_one(aa,bb,right,x_f)
                else:
                    value_final = big_m_infty_norm_2(aa,bb,right,x_f)   
                eta_value.append(value_final)
#             eta_value.append(big_m_infty_norm_2(xi[i],xi[j]))
        
    Big_M_updated[j]=max(eta_value[np.argsort(eta_value)[math.floor(random_size*epsilon)+1]]-b[j],0)
modeltime_big_m_stengthen= time.time() - start



# In[45]:


# start=time.time()
# Big_M_updated = [0]*random_size
# for j in stengthen_difference:
#     eta_value=[]
#     print('current j is', j)
#     for i in stengthen_difference:
#         if sum(xi[i]) <= b[i]:
#             eta_value.append(max(sum(xi[j]) - b[j],0))
# #             print("Special:",i)
#         else:  
#             division_sort = xi[i][np.argsort(np.divide(xi[j],xi[i]))][::-1]

#             division_index = np.argsort(np.divide(xi[j],xi[i]))[x_random_size-np.where(division_sort.cumsum() > b[i])[0][0]:x_random_size]

#             eta_value.append(max(sum(xi[j][division_index]) - b[j],0))
# #             print(i)
#     Big_M_updated[j]=(1+theta)*eta_value[np.argsort(eta_value)[len(stengthen_difference)-math.floor(random_size*epsilon)+len(cut_z_indice)+1]]
# modeltime_big_m_stengthen= time.time() - start



# In[46]:


# subset_indice = np.argsort(quantile_result)[random_size-2*k:random_size]
# # subset_indice = np.argsort(scenario_value_alsox_sharp)[random_size-2*k:random_size]
# stengthen_difference_orginial = [element for element in subset_indice if element not in cut_z_indice]
# stengthen_difference_orginial=stengthen_difference_orginial[::-1]


# In[47]:


# start=time.time()
# cut_third_one = []
# # value_approximation_bound = value_alsox
# print ("Begin to solve second dual bound")

# m = Model()
# m.setParam('Seed', 2)
# m.setParam('TimeLimit', 60*60)
# # Create variables
# m.update()
# # s1 = m.addVars(random_size,name="s1")
# x = m.addVars(x_random_size,lb=0,ub=1, name="x")
# mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")
# w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")
# z = m.addVars(random_size,lb=0,ub=1,name="z")
# slack_t = m.addVars(random_size,lb=0,name="slack_t")
# # slack_t_mu = m.addVars(random_size,lb=0,name="slack_t_mu")
# # slack_t_w = m.addVars(random_size,lb=0,name="slack_t_w")

# m.update()



# # Add objective    
# m.setObjective(0, GRB.MINIMIZE)
# m.update()
# #     m.addConstr(z[3]==1)
# #     m.addConstr(z[4]==1)
# #     m.addConstr(z[5]==1)
# constr1 = m.addConstrs( sum(c[i] * mu[i,j] for i in range(x_random_size)) <= value_approximation_bound * (z[j])  for j in stengthen_difference)
# m.update()
# constr2 = m.addConstrs( sum(c[i] * w[i,j] for i in range(x_random_size)) <= value_approximation_bound * (1-z[j])  for j in stengthen_difference )
# m.update()

# # Add constraints
# m.addConstrs( x[i] == mu[i,j]+w[i,j] for i in range(x_random_size) for j in stengthen_difference )
# m.update()
# # m.addConstrs( mu[j,i]<=1-z[j]  for j in range(x_random_size) for i in range(random_size) )
# m.addConstrs( w[j,i]<=1-z[i]  for j in range(x_random_size) for i in stengthen_difference )
# m.update()
# m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in stengthen_difference)
# m.update()
# # m.addConstrs( sum( xi[i][j]*mu[j,i] for j in range(x_random_size)) <= b[i]*(z[i]) for i in range(random_size)) 
# # m.update()
# # m.addConstrs(sum(xi[i][j]*x[j]  for j in range(x_random_size)) >= b[i] for i in cut_z_indice)

# m.addConstrs(theta*slack_t[i]+ sum( xi[i][j]*mu[j,i] for j in range(x_random_size)) <= b[i]*(z[i]) for i in stengthen_difference) 
# m.update()
# # constr3 = m.addConstrs( sum( xi[i][j]*w[j,i] for j in range(x_random_size)) <= b[i]*z[i] for i in range(random_size))
# m.update()
# #     m.addConstrs(0 <= x[i] <= 1 for i in range(x_random_size))
# m.addConstr(z.sum() >= math.ceil(random_size*(1-epsilon)))
# m.addConstrs(slack_t[i] >= mu[j,i] for j in range(x_random_size)  for i in stengthen_difference)

# # m.addConstrs(slack_t[0] - slack_t_mu[j] - slack_t_w[j] == 0 for j in range(random_size))
# # 

# m.optimize()

# z_solution_store = m.getAttr('x', z)
# mu_solution_store = m.getAttr('x', mu)
# w_solution_store = m.getAttr('x', w)
# x_solution_store = m.getAttr('x', x)

# # subset_indice = np.argsort(quantile_result)[random_size-int(k/2):random_size]
# # subset_indice = range(random_size)
# iteration = 0
# # for i in range(0,10):
# # #     m.addConstr(z[i-1]==1)
# # #     m.getConstrs()[-1]
# #     m.setParam('TimeLimit', 80)  
# #     print(stengthen_difference_update[i],stengthen_difference_update[i])
# #     m.addConstr(z[stengthen_difference_update[i]]+z[stengthen_difference_update[i]]==2)
# #     m.update()
# # #     m.params.OutputFlag=0  

# #     x.start = x_solution_store 
# #     mu.start = mu_solution_store
# #     w.start = w_solution_store
# #     z.start = z_solution_store
# # #     m.setParam('TimeLimit', 100)  
# #     m.params.BarHomogeneous=1 
# # #         m.params.NumericFocus=3
# #     m.optimize()
# #     print(m.status)
# #     print(stengthen_difference_update[i],stengthen_difference_update[i])
# # #     print(m.status)
# # #         if m.status == 4:
# #     if m.status == 12 or m.status == GRB.INFEASIBLE:
# #         print("INFEASIBLE")
# # #         cut_third.append([i,i])
# # #         m.addConstr(sum(xi[stengthen_difference_update[i]][j]*x[j]  for j in range(x_random_size)) >= b[stengthen_difference_update[i]])
# #         cut_z_indice.append(stengthen_difference_update[i])
# #     m.remove(m.getConstrs()[-1])
# #     iteration += 1
# #     print(iteration)
# for i in range(0,math.ceil(k*0.1)):       
# # for i in range(0,k-len(cut_z_indice)-20):
# #     m.addConstr(z[i-1]==1)
# #     m.getConstrs()[-1]
# #     m.setParam('TimeLimit', 10)  
# #     for j in range(1,2):
#     print(stengthen_difference_orginial[i])
#     m.addConstr(z[stengthen_difference_orginial[i]]==1)
#     m.update()
#     m.params.OutputFlag=0  

#     x.start = x_solution_store 
#     mu.start = mu_solution_store
#     w.start = w_solution_store
#     z.start = z_solution_store
#     m.setParam('TimeLimit', 80)  
#     m.params.BarHomogeneous=1 
# #         m.params.NumericFocus=3
#     m.optimize()
#     print(m.status)
#     print(stengthen_difference_orginial[i])
# #     print(m.status)
# #         if m.status == 4:
#     if m.status == 12 or m.status == GRB.INFEASIBLE:
#         print("INFEASIBLE")
#         cut_third_one.append([i,i+j])

#     m.remove(m.getConstrs()[-1])
#     iteration += 1
#     print(iteration)
        
# # for i in range(k-len(cut_z_indice)-40,k-len(cut_z_indice)-10):
# # #     m.addConstr(z[i-1]==1)
# # #     m.getConstrs()[-1]
# # #     m.setParam('TimeLimit', 10)  
# #     for j in range(1,2):
# #         print(stengthen_difference_update[i],stengthen_difference_update[i+j])
# #         m.addConstr(z[stengthen_difference_update[i]]+z[stengthen_difference_update[i+j]]==2)
# #         m.update()
# #     #     m.params.OutputFlag=0  

# #         x.start = x_solution_store 
# #         mu.start = mu_solution_store
# #         w.start = w_solution_store
# #         z.start = z_solution_store
# #         m.setParam('TimeLimit', 40)  
# #         m.params.BarHomogeneous=1 
# # #         m.params.NumericFocus=3
# #         m.optimize()
# #         print(m.status)
# #         print(stengthen_difference_update[i],stengthen_difference_update[i+j])
# #     #     print(m.status)
# # #         if m.status == 4:
# #         if m.status == 12 or m.status == GRB.INFEASIBLE:
# #             print("INFEASIBLE")
# #             cut_third.append([i,i+j])
# #         m.remove(m.getConstrs()[-1])
# #         iteration += 1
# #         print(iteration)
# modeltime_fixing_one = time.time() - start




# In[ ]:





# In[ ]:





# In[48]:


# for i in range(len(cut_third_one)):
#     cut_z_indice.append(stengthen_difference_orginial[cut_third_one[i][0]])
    


# In[49]:


# subset_indice = np.argsort(quantile_result)[random_size-2*k:random_size]
# stengthen_difference_update = [element for element in subset_indice if element not in cut_z_indice]
# stengthen_difference_update=stengthen_difference_update[::-1]
# # value_approximation_bound = value_alsox


# In[50]:


# start=time.time()
# cut_third = []
# # value_approximation_bound = value_alsox
# print ("Begin to solve second dual bound")

# m = Model()
# m.setParam('Seed', 2)
# m.setParam('TimeLimit', 60*60)
# # Create variables
# m.update()
# # s1 = m.addVars(random_size,name="s1")
# x = m.addVars(x_random_size,lb=0,ub=1, name="x")
# mu = m.addVars(x_random_size,random_size,lb=0,ub=1, name="mu")
# w = m.addVars(x_random_size,random_size,lb=0,ub=1, name="w")
# z = m.addVars(random_size,lb=0,ub=1,name="z")
# slack_t = m.addVars(random_size,lb=0,name="slack_t")
# # slack_t_mu = m.addVars(random_size,lb=0,name="slack_t_mu")
# # slack_t_w = m.addVars(random_size,lb=0,name="slack_t_w")

# m.update()



# # Add objective    
# m.setObjective(0, GRB.MINIMIZE)
# m.update()
# #     m.addConstr(z[3]==1)
# #     m.addConstr(z[4]==1)
# #     m.addConstr(z[5]==1)
# constr1 = m.addConstrs( sum(c[i] * mu[i,j] for i in range(x_random_size)) <= value_approximation_bound * (z[j])  for j in stengthen_difference)
# m.update()
# constr2 = m.addConstrs( sum(c[i] * w[i,j] for i in range(x_random_size)) <= value_approximation_bound * (1-z[j])  for j in stengthen_difference )
# m.update()

# # Add constraints
# m.addConstrs( x[i] == mu[i,j]+w[i,j] for i in range(x_random_size) for j in stengthen_difference )
# m.update()
# # m.addConstrs( mu[j,i]<=1-z[j]  for j in range(x_random_size) for i in range(random_size) )
# m.addConstrs( w[j,i]<=1-z[i]  for j in range(x_random_size) for i in stengthen_difference )
# m.update()
# m.addConstrs( mu[j,i] <= z[i] for j in range(x_random_size) for i in stengthen_difference)
# m.update()
# # m.addConstrs( sum( xi[i][j]*mu[j,i] for j in range(x_random_size)) <= b[i]*(z[i]) for i in range(random_size)) 
# # m.update()
# # m.addConstrs(sum(xi[i][j]*x[j]  for j in range(x_random_size)) >= b[i] for i in cut_z_indice)

# m.addConstrs(theta*slack_t[i]+ sum( xi[i][j]*mu[j,i] for j in range(x_random_size)) <= b[i]*(z[i]) for i in stengthen_difference) 
# m.update()
# # constr3 = m.addConstrs( sum( xi[i][j]*w[j,i] for j in range(x_random_size)) <= b[i]*z[i] for i in range(random_size))
# m.update()
# #     m.addConstrs(0 <= x[i] <= 1 for i in range(x_random_size))
# m.addConstr(z.sum() >= math.ceil(random_size*(1-epsilon)))
# m.addConstrs(slack_t[i] >= mu[j,i] for j in range(x_random_size)  for i in stengthen_difference)

# # m.addConstrs(slack_t[0] - slack_t_mu[j] - slack_t_w[j] == 0 for j in range(random_size))
# # 

# m.optimize()

# z_solution_store = m.getAttr('x', z)
# mu_solution_store = m.getAttr('x', mu)
# w_solution_store = m.getAttr('x', w)
# x_solution_store = m.getAttr('x', x)

# # subset_indice = np.argsort(quantile_result)[random_size-int(k/2):random_size]
# # subset_indice = range(random_size)
# iteration = 0
# # for i in range(0,10):
# # #     m.addConstr(z[i-1]==1)
# # #     m.getConstrs()[-1]
# #     m.setParam('TimeLimit', 80)  
# #     print(stengthen_difference_update[i],stengthen_difference_update[i])
# #     m.addConstr(z[stengthen_difference_update[i]]+z[stengthen_difference_update[i]]==2)
# #     m.update()
# # #     m.params.OutputFlag=0  

# #     x.start = x_solution_store 
# #     mu.start = mu_solution_store
# #     w.start = w_solution_store
# #     z.start = z_solution_store
# # #     m.setParam('TimeLimit', 100)  
# #     m.params.BarHomogeneous=1 
# # #         m.params.NumericFocus=3
# #     m.optimize()
# #     print(m.status)
# #     print(stengthen_difference_update[i],stengthen_difference_update[i])
# # #     print(m.status)
# # #         if m.status == 4:
# #     if m.status == 12 or m.status == GRB.INFEASIBLE:
# #         print("INFEASIBLE")
# # #         cut_third.append([i,i])
# # #         m.addConstr(sum(xi[stengthen_difference_update[i]][j]*x[j]  for j in range(x_random_size)) >= b[stengthen_difference_update[i]])
# #         cut_z_indice.append(stengthen_difference_update[i])
# #     m.remove(m.getConstrs()[-1])
# #     iteration += 1
# #     print(iteration)
# for i in range(0,math.ceil(k*0.2)):       
# # for i in range(0,k-len(cut_z_indice)-20):
# #     m.addConstr(z[i-1]==1)
# #     m.getConstrs()[-1]
# #     m.setParam('TimeLimit', 10)  
#     for j in range(1,2):
#         print(stengthen_difference_update[i],stengthen_difference_update[i+j])
#         m.addConstr(z[stengthen_difference_update[i]]+z[stengthen_difference_update[i+j]]==2)
#         m.update()
#         m.params.OutputFlag=0  

#         x.start = x_solution_store 
#         mu.start = mu_solution_store
#         w.start = w_solution_store
#         z.start = z_solution_store
#         m.setParam('TimeLimit', 80)  
#         m.params.BarHomogeneous=1 
# #         m.params.NumericFocus=3
#         m.optimize()
#         print(m.status)
#         print(stengthen_difference_update[i],stengthen_difference_update[i+j])
#     #     print(m.status)
# #         if m.status == 4:
#         if m.status == 12 or m.status == GRB.INFEASIBLE:
#             print("INFEASIBLE")
#             cut_third.append([i,i+j])

#         m.remove(m.getConstrs()[-1])
#         iteration += 1
#         print(iteration)
        
# # for i in range(k-len(cut_z_indice)-40,k-len(cut_z_indice)-10):
# # #     m.addConstr(z[i-1]==1)
# # #     m.getConstrs()[-1]
# # #     m.setParam('TimeLimit', 10)  
# #     for j in range(1,2):
# #         print(stengthen_difference_update[i],stengthen_difference_update[i+j])
# #         m.addConstr(z[stengthen_difference_update[i]]+z[stengthen_difference_update[i+j]]==2)
# #         m.update()
# #     #     m.params.OutputFlag=0  

# #         x.start = x_solution_store 
# #         mu.start = mu_solution_store
# #         w.start = w_solution_store
# #         z.start = z_solution_store
# #         m.setParam('TimeLimit', 40)  
# #         m.params.BarHomogeneous=1 
# # #         m.params.NumericFocus=3
# #         m.optimize()
# #         print(m.status)
# #         print(stengthen_difference_update[i],stengthen_difference_update[i+j])
# #     #     print(m.status)
# # #         if m.status == 4:
# #         if m.status == 12 or m.status == GRB.INFEASIBLE:
# #             print("INFEASIBLE")
# #             cut_third.append([i,i+j])
# #         m.remove(m.getConstrs()[-1])
# #         iteration += 1
# #         print(iteration)
# modeltime_fixing_two = time.time() - start




# In[51]:


# cut_third


# In[52]:


# len(cut_third)


# In[53]:


# modeltime_fixing = modeltime_fixing_one + modeltime_fixing_two


# In[54]:


# modeltime_fixing


# In[55]:


len(cut_z_indice)


# In[56]:


# number_of_cuts = len(cut_z_indice) + len(cut_third)


# In[57]:


# stengthen_difference = [element for element in range(random_size) if element not in cut_z_indice]
# 


# In[58]:


# Big-M Method
#Big_M = [None] * random_size
# Big_M = []
# for i in range(random_size):
   
#     Big_M.append(math.ceil(sum(x for x in a[i,:] if x > 0) -bb[i])) 
# Big_M = [30]*random_size

def Model3_var():
    print ("Begin to solve model 3")
    
    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 4*60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(stengthen_difference,vtype=GRB.BINARY,name="z")
#     y = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="y")
#     s = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="s")
#     gamma = m.addVars(1,lb=0,name="gamma")
    slack_t = m.addVar(lb=0,name="slack_t")
   
    m.update()

    # Set functions
    
    #obj = 0
#     g_x_xi = [None] * len(z)
#     f_x = 0
#     g_x = 0
    
    # Set functions
#     for i in range(x_random_size):
      
#         f_x += c[i]*x[i]
#         m.update()
    
#     for i in range(random_size):
#         g_x_si[i] = g_x +       SSS   
#     for i in range(random_size):
#      #   g_x_si[i] = g_x + si[i] 
#         g_x = 0
#         for j in range(x_random_size):
#             g_x += xi[i][j]*x[j]
#             m.update()
#         g_x_xi[i] = g_x  
    
    # Add objective    
    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= v_q)
    
    m.update()
    
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= value_alsox_sharp)
    m.update()
    
    m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
    
    m.update()
    # Add constraints
   
#     m.addConstrs(sum(xi[i][j]*x[j]  for j in range(x_random_size)) >= b[i] for i in cut_z_indice)
    
    m.update()
    
    m.addConstrs( slack_t*theta + sum(xi[i][j]*x[j]  for j in range(x_random_size)) <= b[i]   + Big_M_updated[i]*z[i] for i in stengthen_difference)
    m.update()


    m.addConstr(z.sum() <= math.floor(random_size*epsilon)-len(cut_z_indice))
    m.update()
#     m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]]>= 1 for i in range(len(cut_third)))

    # Solve the problem
    m.update()
    
    for i in range(x_random_size):

        x[i].start = solution_alsox_sharp[i]
#     m.params.BarHomogeneous=1     
    m.optimize()
    

    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_VaR = m.getAttr('x', z)
    x_VaR = m.getAttr('x', x)
    
    #return pp,ii,kk
    
    return m.objVal,z_VaR,x_VaR
    


# In[59]:


start=time.time()
Model3_var_value = Model3_var()
z_VaR = Model3_var_value[1]
x_VaR = Model3_var_value[2]
# value_z = Model3_cut()[0]
# kk=Model3_cut()[1]
# obj_cut = 0
# for i in range(x_random_size):
#     obj_cut += c[i]*kk[i]
modeltime_Model3_var = time.time() - start
# print(f'The value of optimal value t is exactly {obj_cut_stengthen_updated:.5f}.')
print(f'The value of running time is approximately {modeltime_Model3_var:.3f} s.')
# print(f'The value of total running time is approximately {modeltime_cut_stengthen_updated+modeltime_big_m_stengthen+modeltime_fixing:.3f} s.')


# In[60]:


modeltime_preprocessing = modeltime_quantile + modeltime_cvar


# In[61]:


print(f'The value of modeltime preprocessing is approximately {modeltime_preprocessing:.3f} s.')


# In[62]:


print(f'The value of modeltime modeltime_big_m_stengthen is approximately {modeltime_big_m_stengthen:.3f} s.')


# In[63]:


len(cut_z_indice)


# In[64]:


# print(f'The value of modeltime fixing is approximately {modeltime_fixing:.3f} s.')


# In[65]:


# print(f'The number of cuts is approximately {number_of_cuts:.3f} s.')


# In[66]:


# strengthen
start=time.time()
Big_M_strength_only = [0]*random_size
for j in range(random_size):
    eta_value=[]
    print('current j is', j)
    for i in range(random_size):
        if theta*math.sqrt(x_random_size)+sum(xi[i]) <= b[i]:
            eta_value.append(max(sum(xi[j]) - b[j],0))
        else:
            if j==i:
                eta_value.append(max(sum(xi[j]) - b[j],0))
            else:
                aa = xi[j]
                bb = xi[i]
                right = b[i]
#                 eta_value.append(big_m_infty_norm_2(aa,bb,right))
                x_l_init = 0
                x_u_init = 2.0
                tolerance = 5e-2
                max_iter = 50
                x_f, iter_ctr = goldenopt( x_l_init, x_u_init, tolerance, max_iter)
                if x_f<=1:
                    value_final = big_m_infty_norm_2_less_than_one(aa,bb,right,x_f)
                else:
                    value_final = big_m_infty_norm_2(aa,bb,right,x_f)   
                eta_value.append(value_final)
#             eta_value.append(big_m_infty_norm_2(xi[i],xi[j]))
        
    Big_M_strength_only[j]=max(eta_value[np.argsort(eta_value)[math.floor(random_size*epsilon)+1]]-b[j],0)
modeltime_big_m_stengthen_only= time.time() - start



# In[67]:


# start=time.time()
# Big_M_strength_only = []
# for j in range(random_size):
#     eta_value=[]
#     for i in range(random_size):
#         if sum(xi[i]) <= b[i]:
#             eta_value.append(max(sum(xi[j]) - b[j],0))
# #             print("Special:",i)
#         else:  
#             division_sort = xi[i][np.argsort(np.divide(xi[j],xi[i]))][::-1]

#             division_index = np.argsort(np.divide(xi[j],xi[i]))[x_random_size-np.where(division_sort.cumsum() > b[i])[0][0]:x_random_size]

#             eta_value.append(max(sum(xi[j][division_index]) - b[j],0))
# #             print(i)
#     Big_M_strength_only.append((1+theta)*eta_value[np.argsort(eta_value)[random_size-math.floor(random_size*epsilon)+1]])
# modeltime_big_m_stengthen_only = time.time() - start


# In[68]:


# Big-M Method
#Big_M = [None] * random_size
# Big_M = []
# for i in range(random_size):
   
#     Big_M.append(math.ceil(sum(x for x in a[i,:] if x > 0) -bb[i])) 
# Big_M = [30]*random_size

def Model3_strength_only():
    #Model
    print ("Begin to solve model 3")

    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 4*60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
    slack_t = m.addVar(lb=0,name="slack_t")
#     y = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="y")
#     s = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="s")
#     gamma = m.addVars(1,lb=0,name="gamma")
#     lambda_1 = m.addVars(1,lb=0,name="lambda_1")
   
    m.update()

    # Set functions
    
    #obj = 0
#     g_x_xi = [None] * len(z)
#     f_x = 0
#     g_x = 0
    
    # Set functions
#     for i in range(x_random_size):
      
#         f_x += c[i]*x[i]
#         m.update()
    
#     for i in range(random_size):
#         g_x_si[i] = g_x +       SSS   
#     for i in range(random_size):
#      #   g_x_si[i] = g_x + si[i] 
#         g_x = 0
#         for j in range(x_random_size):
#             g_x += xi[i][j]*x[j]
#             m.update()
#         g_x_xi[i] = g_x  
    
    # Add objective    
    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    
#     m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= lower_bound_second_dual)
    
#     m.update()
    # Add constraints
    m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
#     m.addConstrs(slack_t >= x[i] for i in range(x_random_size))
    m.update()
   
#     m.addConstrs(sum(xi[i][j]*x[j]  for j in range(x_random_size)) >= b[i] for i in cut_z_indice)
#     m.update()
    
    m.addConstrs(slack_t*theta+sum(xi[i][j]*x[j]  for j in range(x_random_size)) <= b[i]   + Big_M_strength_only[i]*z[i] for i in range(random_size))
    m.update()


    m.addConstr(z.sum() <= math.floor(random_size*epsilon))
    m.update()
#     m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]]>= 1 for i in range(len(cut_third)))

    # Solve the problem
    m.update()
#     m.params.BarHomogeneous=1     
    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_strength_only = m.getAttr('x', z)
    x_strength_only = m.getAttr('x', x)
    
    #return pp,ii,kk
    
    return m.objVal,z_strength_only,x_strength_only


# In[69]:


start=time.time()
Model3_strength_only_value = Model3_strength_only()
# z_VaR = Model3_var_value[1]
# x_VaR = Model3_var_value[2]
# value_z = Model3_cut()[0]
# kk=Model3_cut()[1]
# obj_cut = 0
# for i in range(x_random_size):
#     obj_cut += c[i]*kk[i]
modeltime_Model3_strength_only = time.time() - start
# print(f'The value of optimal value t is exactly {obj_cut_stengthen_updated:.5f}.')
print(f'The value of running time is approximately {modeltime_Model3_strength_only:.3f} s.')
# print(f'The value of total running time is approximately {modeltime_cut_stengthen_updated+modeltime_big_m_stengthen+modeltime_fixing:.3f} s.')


# In[70]:


print(f'The value of Big-M stengthening only is approximately {modeltime_big_m_stengthen_only:.3f} s.')


# In[71]:


Big_M = []
for i in range(random_size):
   
    Big_M.append(math.ceil(sum(x for x in xi[i][:] if x > 0) -b[i])) 
# Big_M = [30]*random_size


# In[72]:


# Big-M Method
#Big_M = [None] * random_size
# Big_M = []
# for i in range(random_size):
   
#     Big_M.append(math.ceil(sum(x for x in a[i,:] if x > 0) -bb[i])) 
# Big_M = [30]*random_size

def Model3_cut_only():
    #Model
    print ("Begin to solve model 3")

    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 4*60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(stengthen_difference,vtype=GRB.BINARY,name="z")
    slack_t = m.addVar(lb=0,name="slack_t")
#     y = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="y")
#     s = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="s")
#     gamma = m.addVars(1,lb=0,name="gamma")
#     lambda_1 = m.addVars(1,lb=0,name="lambda_1")
   
    m.update()

    # Set functions
    
    #obj = 0
#     g_x_xi = [None] * len(z)
#     f_x = 0
#     g_x = 0
    
    # Set functions
#     for i in range(x_random_size):
      
#         f_x += c[i]*x[i]
#         m.update()
    
#     for i in range(random_size):
#         g_x_si[i] = g_x +       SSS   
#     for i in range(random_size):
#      #   g_x_si[i] = g_x + si[i] 
#         g_x = 0
#         for j in range(x_random_size):
#             g_x += xi[i][j]*x[j]
#             m.update()
#         g_x_xi[i] = g_x  
    
    # Add objective    
    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
#     m.addConstrs(slack_t >= x[i] for i in range(x_random_size))
    m.update()
    
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) >= v_q)
    
    m.update()
        
    m.addConstr(sum(c[i]*x[i]  for i in range(x_random_size)) <= value_alsox_sharp)
    m.update()
    # Add constraints
   
    m.addConstrs(sum(xi[i][j]*x[j]  for j in range(x_random_size)) >= b[i] for i in cut_z_indice)
    m.update()
    
    m.addConstrs(slack_t*theta + sum(xi[i][j]*x[j]  for j in range(x_random_size)) <= b[i]   + Big_M[i]*z[i] for i in stengthen_difference)
    m.update()


    m.addConstr(z.sum() <= math.floor(random_size*epsilon)-len(cut_z_indice))
    m.update()
#     m.addConstrs(z[stengthen_difference_update[cut_third[i][0]]]+z[stengthen_difference_update[cut_third[i][1]]]>= 1 for i in range(len(cut_third)))

    # Solve the problem
    m.update()
    for i in range(x_random_size):

        x[i].start = solution_alsox_sharp[i]
#     m.params.BarHomogeneous=1     
    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_cut_only = m.getAttr('x', z)
    x_cut_only = m.getAttr('x', x)
    
    #return pp,ii,kk
    
    return m.objVal,z_cut_only,x_cut_only


# In[73]:


# start=time.time()
# Model3_cut_only_value = Model3_cut_only()
# # z_VaR = Model3_var_value[1]
# # x_VaR = Model3_var_value[2]
# # value_z = Model3_cut()[0]
# # kk=Model3_cut()[1]
# # obj_cut = 0
# # for i in range(x_random_size):
# #     obj_cut += c[i]*kk[i]
# modeltime_Model3_cut_only = time.time() - start
# # print(f'The value of optimal value t is exactly {obj_cut_stengthen_updated:.5f}.')
# print(f'The value of running time is approximately {modeltime_Model3_cut_only:.3f} s.')
# # print(f'The value of total running time is approximately {modeltime_cut_stengthen_updated+modeltime_big_m_stengthen+modeltime_fixing:.3f} s.')


# In[74]:


# print(f'The value of modeltime fixing is approximately {modeltime_fixing:.3f} s.')


# In[ ]:





# In[ ]:





# In[75]:


# Big-M Method
#Big_M = [None] * random_size
# Big_M = []
# for i in range(random_size):
   
#     Big_M.append(math.ceil(sum(x for x in a[i,:] if x > 0) -bb[i])) 
# Big_M = [30]*random_size

def Model3():
    #Model
    print ("Begin to solve model 3")

    m = Model()
    m.setParam('Seed', 2)
    m.setParam('TimeLimit', 4*60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
    slack_t = m.addVar(lb=0,name="slack_t")
#     y = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="y")
#     s = m.addVars(stengthen_difference,lb=-GRB.INFINITY,name="s")
#     gamma = m.addVars(1,lb=0,name="gamma")
#     lambda_1 = m.addVars(1,lb=0,name="lambda_1")
   
    m.update()

    # Set functions
    
    #obj = 0
#     g_x_xi = [None] * len(z)
#     f_x = 0
#     g_x = 0
    
    # Set functions
#     for i in range(x_random_size):
      
#         f_x += c[i]*x[i]
#         m.update()
    
#     for i in range(random_size):
#         g_x_si[i] = g_x +       SSS   
#     for i in range(random_size):
#      #   g_x_si[i] = g_x + si[i] 
#         g_x = 0
#         for j in range(x_random_size):
#             g_x += xi[i][j]*x[j]
#             m.update()
#         g_x_xi[i] = g_x  
    
    # Add objective    
    m.setObjective(sum(c[i]*x[i]  for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    # Add constraints
    m.addConstr(slack_t*slack_t >= sum(x[i]*x[i] for i in range(x_random_size)))
#     m.addConstrs(slack_t >= x[i] for i in range(x_random_size))
    m.update()
   
    
    m.addConstrs(slack_t*theta+ sum(xi[i][j]*x[j]  for j in range(x_random_size)) <= b[i]   + Big_M[i]*z[i] for i in range(random_size))
    m.update()


    m.addConstr(z.sum() <= math.floor(random_size*epsilon))
    m.update()

    # Solve the problem
    m.update()
#     m.params.BarHomogeneous=1     
    m.optimize()
    
    # Store solutions
    #pp = m.getAttr('x', z)
    z_org= m.getAttr('x', z)
    x_org = m.getAttr('x', x)
    
    #return pp,ii,kk
    
    return m.objVal,z_org,x_org


# In[76]:


start=time.time()
Model3_big=Model3()
obj_big_m = Model3_big[0]
value_z = Model3_big[1]
kk=Model3_big[2]
# obj = 0
# for i in range(x_random_size):
#     obj += c[i]*kk[i]
modeltime_big_m = time.time() - start
print(f'The value of optimal value t is exactly {obj_big_m:.5f}.')
print(f'The value of running time is approximately {modeltime_big_m:.3f} s.')


# In[ ]:





# In[ ]:




