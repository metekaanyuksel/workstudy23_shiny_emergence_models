import pandas as pd 
import numpy as np
from scipy.integrate import odeint
import time 
from functools import partial
from multiprocessing import Pool

def model(time, state, parms):
	## state and parms are dictionaries 
    S, R, IG00, IG01, IG10, IG11 = tuple(state) ## Assuming this is inputted in the right order 
    betaG00, betaG01, betaG10, betaG11, gamma,delta,rho,b,d,alpha, mu00_01, mu11_10, mu10_11, mu01_00,  mu11_01,mu01_11, mu10_00, mu00_10, sigma = tuple(list(parms.values()))

    dS = b*(S+IG00+IG01+IG10+IG11+R)*(1-alpha*(S+IG00+IG01+IG10+IG11+R)) - (rho*(S+IG00+IG01+IG10+IG11+R)+d)*S + delta*R - S*(betaG00*IG00+betaG01*IG01+betaG10*IG10+betaG11*IG11)
    dR = gamma*(IG00+IG01+IG10+IG11) - (rho*(S+IG00+IG01+IG10+IG11+R)+d)*R - delta*R
    dIG00 = betaG00*S*IG00 - mu00_01*IG00 - mu00_10*IG00 + mu10_00*IG10 + mu01_00*IG01 - gamma*IG00 - (rho*(S+IG00+IG01+IG10+IG11+R)+d)*IG00 - sigma*IG00*IG11 + sigma*IG01*IG10
    dIG01 = betaG01*S*IG01 - mu01_00*IG01 - mu01_11*IG01 + mu11_01*IG11 + mu00_01*IG00 - gamma*IG01 - (rho*(S+IG00+IG01+IG10+IG11+R)+d)*IG01 - sigma*IG01*IG10 + sigma*IG00*IG11
    dIG10 = betaG10*S*IG10 - mu10_00*IG10 - mu10_11*IG10 + mu00_10*IG00 + mu11_10*IG11 - gamma*IG10 - (rho*(S+IG00+IG01+IG10+IG11+R)+d)*IG10 - sigma*IG10*IG01 + sigma*IG11*IG00
    dIG11 = betaG11*S*IG11 - mu11_01*IG11 - mu11_10*IG11 + mu01_11*IG01 + mu10_11*IG10 - gamma*IG11 - (rho*(S+IG00+IG01+IG10+IG11+R)+d)*IG11 - sigma*IG11*IG00 + sigma*IG10*IG01

    return [dS, dR, dIG00, dIG01, dIG10, dIG11]


def fixed_params(N, I0_fraction = 1, time_max=365*25):
    I0 = 1
    state = {'S': N - I0, 'R': 0, 'IG00': I0, 'IG01': 0, 'IG10': 0, 'IG11': 0}
    # time_max = 10 #365*10
    times = np.arange(0, time_max + 1, 1)
    return state, times

### define parameters under following assumptions: symmetric mutation, recomb. w/ no within host differences in fitness
def define_params(params_dict):

    gamma = params_dict['turnover']*params_dict['gamma']
    delta = params_dict['turnover']*params_dict['delta']

    alpha = (params_dict['b'] - params_dict['turnover']) / (params_dict['b'] * params_dict['Neq'])  # density-dependent birth rate
    rho = (params_dict['turnover'] - params_dict['d']) / params_dict['Neq']  # density-dependent death rate
    beta = params_dict['R0'] * \
            (rho * (params_dict['b'] + gamma) \
            + alpha * params_dict['b'] * (gamma \
            + params_dict['d'])) / (params_dict['b'] - params_dict['d'])  # transmission rate for 00 genotype

    h = 1  # scaling to adjust cost of carrying 1s at both loci in the reservoir
    betaG00 = beta
    betaG01 = beta * (1 - params_dict['c2'])  # transmission rate for 01 genotype
    betaG10 = beta * (1 - params_dict['c1'])  # transmission rate for 10 genotype
    betaG11 = beta * (1 - params_dict['c1']) * (1 - h * params_dict['c2'])  # transmission rate for 11 genotype
    mu00_01 = params_dict['mu2']
    mu11_10 = params_dict['nu2']
    mu10_11 = params_dict['mu2']
    mu01_00 = params_dict['nu2']
    mu11_01 = params_dict['nu1']
    mu01_11 = params_dict['mu1']
    mu10_00 = params_dict['nu1']
    mu00_10 = params_dict['mu1']

    params_to_use = {
        'betaG00': betaG00,
        'betaG01': betaG01,
        'betaG10': betaG10,
        'betaG11': betaG11,
        'gamma': gamma,
        'delta': delta,
        'rho': rho,
        'b': params_dict['b'],
        'd': params_dict['d'],
        'alpha': alpha,
        'mu00_01': mu00_01,
        'mu11_10': mu11_10,
        'mu10_11': mu10_11,
        'mu01_00': mu01_00,
        'mu11_01': mu11_01,
        'mu01_11': mu01_11,
        'mu10_00': mu10_00,
        'mu00_10': mu00_10,
        'sigma': params_dict['sigma']
    }

    return params_to_use

# params_to_use = define_params(turnover = 0.001)

# sigma_values = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1]
# turnover_values = np.arange(0.0001, 0.0002, 1) #np.arange(0.0001, 0.001, 0.00001)
# gamma_values = np.arange(0.001, 0.002, 0.001) #np.arange(0.001, 0.1, 0.001)

def create_meshgrid(p1,p1_name, p2, p2_name):
    '''
    input from slider range - p1 and p2 are tuple (min,max)
    '''
    p1_values = np.linspace(p1[0],p1[1],20)
    p2_values = np.linspace(p2[0],p2[1],20)
    p1, p2 = np.meshgrid(p1_values, p2_values)
    p1 = p1.reshape(-1,1)
    p2 = p2.reshape(-1,1)
    return [p1_name,p2_name],np.hstack((p1, p2))


def turnover_loop(param_name, values_all, N,time_max, i):
    default = {'mu1':0.0001, 'mu2':0.0001, 'nu1':0.0001, 'nu2':0.0001,
                'sigma': 0, 'c1':0.9, 'c2':0.9, 'gamma':1/30,
                'delta':1/90, 'b':1/365, 'd':1/(10*365),
                'turnover':0.00049, 'Neq' : 250, 'R0':10}
    if 'I0_fraction' not in param_name:
        default[param_name[0]] = values_all[i,0]
        default[param_name[1]] = values_all[i,1]
        params_to_use = define_params(default)
        state, times = fixed_params(N,time_max=time_max)
    else:
        I0_frac_idx = param_name.index('I0_fraction') ## either 0 or 1 
        p_idx = 1 - I0_frac_idx ## the other param's index
        I0_frac = values_all[i, I0_frac_idx]
        default[param_name[p_idx]] = values_all[i,p_idx]
        params_to_use = define_params(default)
        state, times = fixed_params(N, I0_frac,time_max)

    def ode_func(state, time):
        return model(time, state, params_to_use)

    out = odeint(ode_func, list(state.values()), times)
    # print(out.shape)
    # print(times.shape)
    out = np.column_stack((times, out)) ### odeint does not output time directly like ode in R
    out_df = pd.DataFrame(out, columns=["time", "S", "R", "IG00", "IG01", "IG10", "IG11"])
    out_df["I"] = out_df["IG00"] + out_df["IG10"] + out_df["IG01"] + out_df["IG11"]
    out_df["N"] = out_df["S"] + out_df["I"] + out_df["R"]
    out_df["sigmaIG11_in"] = default['sigma'] * out_df["IG01"] * out_df["IG10"]
    out_df["sigmaIG11_out"] = -default['sigma'] * out_df["IG11"] * out_df["IG00"]
    out_df["sigmaIG11"] = out_df["sigmaIG11_in"] + out_df["sigmaIG11_out"]
    out_df["x11"] = out_df["IG11"] / out_df["I"]
    out_df["x10"] = out_df["IG10"] / out_df["I"]
    out_df["x01"] = out_df["IG01"] / out_df["I"]
    out_df["x00"] = out_df["IG00"] / out_df["I"]
    out_df["p1"] = (out_df["IG10"] + out_df["IG11"]) / out_df["I"]
    out_df["p2"] = (out_df["IG01"] + out_df["IG11"]) / out_df["I"]
    out_df["D"] = out_df["x11"] - out_df["p1"] * out_df["p2"]
    out_df["sigma"] = default['sigma']
    out_df["turnover"] = default['turnover']
    out_df["gamma"] = default['gamma']

    out_new = out_df[out_df["time"] == time_max].reset_index(drop=True)

    return out_new

def steady_state(N, p1,p1_name,p2,p2_name,time_max=10):
    turnover_param_name, turnover_values_to_loop_over = create_meshgrid(p1, p1_name, p2, p2_name)
    for i in range(turnover_values_to_loop_over.shape[0]):
        pool = Pool(processes = 1) ## specify the # of cores to use 
    
        # Apply the turnover_loop function to each row of turnover_values_to_loop_over in parallel
        data = pool.starmap(turnover_loop,
                        [(turnover_param_name,
                         turnover_values_to_loop_over, N, time_max,
                        i) for i in range(turnover_values_to_loop_over.shape[0])])
        
        # Concatenate the resulting dataframes into a single dataframe
        dat = pd.concat(data, ignore_index=True)
        
        # Close the pool of workers
        pool.close()
        # print(dat[['I','S','R']])
        pool.join()
        return dat,turnover_values_to_loop_over 

def turnover_loop_by_time(param_name, values_all,N, time_max, i):
    default = {'mu1':0.0001, 'mu2':0.0001, 'nu1':0.0001, 'nu2':0.0001,
                'sigma': 0, 'c1':0.4, 'c2':0.6, 'gamma':1/30,
                'delta':1/90, 'b':1/365, 'd':1/(10*365),
                'turnover':0.001, 'Neq' : 250, 'R0':3}
    if param_name != 'I0_fraction':
        # print(values_all[i, 0])
        default[param_name] = values_all[i,0]
        # print(default[param_name])
        params_to_use = define_params(default)
        state, times = fixed_params(N, time_max=time_max)
    else:
        params_to_use = define_params(default)
        state, times = fixed_params(N,values_all[i,0], time_max)
    print(params_to_use, default)
    def ode_func(state, time):
        return model(time, state, params_to_use)

    out = odeint(ode_func, list(state.values()), times)
    out = np.column_stack((times, out)) ### odeint does not output time directly like ode in R
    # print(out)
    out_df = pd.DataFrame(out, columns=["time", "S", "R", "IG00", "IG01", "IG10", "IG11"])
    out_df["I"] = out_df["IG00"] + out_df["IG10"] + out_df["IG01"] + out_df["IG11"]
    out_df["N"] = out_df["S"] + out_df["I"] + out_df["R"]
    out_df["sigmaIG11_in"] = default['sigma'] * out_df["IG01"] * out_df["IG10"]
    out_df["sigmaIG11_out"] = -default['sigma'] * out_df["IG11"] * out_df["IG00"]
    out_df["sigmaIG11"] = out_df["sigmaIG11_in"] + out_df["sigmaIG11_out"]
    out_df["x11"] = out_df["IG11"] / out_df["I"]
    out_df["x10"] = out_df["IG10"] / out_df["I"]
    out_df["x01"] = out_df["IG01"] / out_df["I"]
    out_df["x00"] = out_df["IG00"] / out_df["I"]
    out_df["p1"] = (out_df["IG10"] + out_df["IG11"]) / out_df["I"]
    out_df["p2"] = (out_df["IG01"] + out_df["IG11"]) / out_df["I"]
    out_df["D"] = out_df["x11"] - out_df["p1"] * out_df["p2"]
    out_df["sigma"] = default['sigma']
    out_df["turnover"] = default['turnover']
    out_df["gamma"] = default['gamma']
    out_df[param_name] = values_all[i,0]
    return out_df

def simu_by_time(N, p1,p1_name,time_max=10):
    # state, times = fixed_params(N, I0,time_max)
    # turnover_param_name, turnover_values_to_loop_over = create_meshgrid(p1, p1_name, p2, p2_name)
    turnover_param_name = p1_name
    turnover_values_to_loop_over = np.array([p1])
    turnover_values_to_loop_over = turnover_values_to_loop_over.reshape(-1,1)
    for i in range(turnover_values_to_loop_over.shape[0]):
        pool = Pool(processes = 1) ## specify the # of cores to use 
        print('start {}:'.format(i), time.time(), turnover_values_to_loop_over)
        # Apply the turnover_loop function to each row of turnover_values_to_loop_over in parallel
        data = pool.starmap(turnover_loop_by_time,
                        [(turnover_param_name,
                         turnover_values_to_loop_over, N, time_max,
                        i) for i in range(turnover_values_to_loop_over.shape[0])])
        print('end:', time.time())
        # Concatenate the resulting dataframes into a single dataframe
        dat = pd.concat(data, ignore_index=True)
        
        # Close the pool of workers
        pool.close()
        pool.join()
    return dat,turnover_values_to_loop_over 



# Create a pool of workers
# if __name__  == '__main__':
#     sigma_values = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1]
#     turnover_values = np.arange(0.0001, 0.0002, 0.0001) #np.arange(0.0001, 0.001, 0.00001)
#     gamma_values = np.arange(0.001, 0.004, 0.001) #np.arange(0.001, 0.1, 0.001)
#     # Create a meshgrid of the three arrays
#     sigma, turnover, gamma = np.meshgrid(sigma_values, turnover_values, gamma_values)

#     # Reshape the arrays into a single column
#     sigma = sigma.reshape(-1, 1)
#     turnover = turnover.reshape(-1, 1)
#     gamma = gamma.reshape(-1, 1)

#     # Combine the arrays column-wise
#     turnover_values_to_loop_over = np.hstack((sigma, turnover, gamma))
#     print(turnover_values_to_loop_over)
    # start_time = time.time()
    # pool = Pool(processes = 1) ## specify the # of cores to use 
    
    # # Apply the turnover_loop function to each row of turnover_values_to_loop_over in parallel
    # data = pool.map(turnover_loop, range(turnover_values_to_loop_over.shape[0]))
    
    # # Concatenate the resulting dataframes into a single dataframe
    # dat = pd.concat(data, ignore_index=True)
    
    # # Close the pool of workers
    # pool.close()
    # pool.join()
    
    # end_time = time.time()
    
    # elapsed_time = end_time - start_time
    
    # print('time',elapsed_time)
