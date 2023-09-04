import numpy as np
from scipy.optimize import fsolve

from sympy import symbols, Eq, solve


import random
import matplotlib.pyplot as plt


def generate_r0(m,wild_type_r0, evolved_r0, mode = 'Jackpot',valley = 0.7):
	if mode == 'Jackpot':
		all_r0 = [wild_type_r0 for i in range(m-1)]
		all_r0.append(evolved_r0)
	elif mode == 'Additive':
		all_r0 = []
		avg_fitness = (wild_type_r0+evolved_r0)/(2*(m-1))
		for i in range(m-1):
			all_r0.append(wild_type_r0+(i*avg_fitness))
		all_r0.append(evolved_r0)
	else:
		valley_r0 = valley*wild_type_r0
		all_r0 = [wild_type_r0]
		for i in range(m-2):
			all_r0.append(valley_r0)
		all_r0.append(evolved_r0)
	return all_r0
def equations(variables, mu, r0_list):
	## variables will be list of m elements, one for each s_i
	## we generate the pgf on the fly 
	all_pgf = []
	for idx, var in enumerate(variables):
		if idx != len(variables)-1:
			pgf =  np.exp(-(1-mu)*r0_list[idx]*(1-var))*np.exp(-mu*r0_list[idx]*(1-variables[idx+1]))-var
			all_pgf.append(pgf)
		else:
			pgf = np.exp(-r0_list[idx]*(1-var))-var
			all_pgf.append(pgf)
	return all_pgf

def calculate_emergence_prob(m,wild_type_r0, evolved_r0,mu, mode = 'Jackpot',valley_depth = 0.7):
	all_r0 = generate_r0(m,wild_type_r0, evolved_r0, mode, valley_depth)
	initial_guess = [0.1 for i in range(m)]
	result = fsolve(equations, initial_guess, args=(mu, all_r0))
	print(all_r0, result)
	return (1 - result[0])

def simulate(m,mu, wild_type_r0, evolved_r0,mode, nindv,ngen, nsim):
	all_sim_results = {}
	r0_list = generate_r0(m,wild_type_r0, evolved_r0, mode = 'Jackpot',valley = 0.7)
	if nsim == 1:
		for i in range(len(r0_list)):
			all_sim_results[i] = [0 for _ in range(ngen)]
		all_sim_results[0][0] = 0 
	else: ## collect the number of type m in each generation by simulation
		all_sim_results['type_m'] = [] 
	all_sim_results['T_to_E'] = []
	all_sim_results['Emerged'] = []
	all_sim_results['Extinct'] = []
	for k in range(nsim):
		sim_dict = {}
		time_to_emerge = -1
		time_to_extinction = -1
		emerge_flag = -1
		## initialize sim dictionary
		for i in range(len(r0_list)):
			sim_dict[i] = [0 for _ in range(ngen)]
		sim_dict[0][0] = nindv
		total_pop = nindv
		for j in range(ngen-1): ## for each generation
			if total_pop >= 10000:
				break
			else:
				for i in range(len(r0_list)): ## for each type i
					if sim_dict[i][j] != 0: ## if this generation has individuals with type j 
						total_pop = sim_dict[i][j]
						if total_pop >= 10000:
							break 
						else:
							for m in range(sim_dict[i][j]): # for each individual with type i at generation j 
								if i != len(r0_list)-1:
									same_i_inf = np.random.poisson((1-mu)*r0_list[i])
									next_i_inf = np.random.poisson(mu*r0_list[i])
									# print(i,j,same_i_inf, next_i_inf)
									sim_dict[i][j+1] += same_i_inf
									sim_dict[i+1][j+1] += next_i_inf
									# print(sim_dict)
									if i+1 == len(r0_list) - 1 and sim_dict[i+1][j+1] >0:
										time_to_emerge = j+1 
										emerge_flag = 1
								else:
									same_i_inf = np.random.poisson(r0_list[i])
									sim_dict[i][j+1] += same_i_inf

				cnt = 0 
				for i in range(len(r0_list)):
					cnt += sim_dict[i][j+1]
				if cnt == 0: ## if every type has gone extinct, break the simulation
					time_to_extinction = j+1
					break
		if nsim ==1:
			for i in range(len(r0_list)): ## for each r0
				all_sim_results[i] = np.array(all_sim_results[i]) + np.array(sim_dict[i])
		else:
			all_sim_results['type_m'].append(sim_dict[len(r0_list)-1])
		if time_to_emerge != -1:
			all_sim_results['T_to_E'].append(time_to_emerge) ## time to emergence in every simulation 
		all_sim_results['Emerged'].append(emerge_flag) ## if emerged in this simulation 
		if time_to_extinction != -1:
			all_sim_results['Extinct'].append(time_to_extinction)
	if nsim == 1:
		for i in range(len(r0_list)): ## for each generation
			all_sim_results[i] = all_sim_results[i]/nsim ## averaging over all simulations
	return all_sim_results 

def mu_r0_combo(mu_list, wild_r0_list):
	all_combos = [(mu, r0) for mu in mu_list for r0 in wild_r0_list]
	return all_combos 

def simulate_all_combos(m,mu_list, wild_r0_list, evolved_r0,mode, nindv,ngen, nsim):
	combos = mu_r0_combo(mu_list, wild_r0_list)
	combo_result = {}
	for pair in combos:
		mu, wild_type_r0 = float(pair[0]), float(pair[1])
		sim_reuslt = simulate(m,mu, wild_type_r0, evolved_r0,mode, nindv,ngen, nsim)
		combo_result[pair] = sim_reuslt
	return combo_result

# Call fsolve
if __name__ == '__main__':
	m,mu, wild_type_r0, evolved_r0,mode, nindv,ngen, nsim = 5,0.01,3,5,'Jackpot',10,100,5
	all_r = simulate(m,mu, wild_type_r0, evolved_r0,mode, nindv,ngen, nsim)
	print(all_r)


