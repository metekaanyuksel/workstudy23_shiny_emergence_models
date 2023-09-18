from shiny import App, render, ui,reactive
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import deterministic_recombination as dr
import antia as an 
import matplotlib.colors as colors
from matplotlib.gridspec import SubplotSpec


app_ui = ui.page_fluid(
    ui.panel_title("Models of evolutionary emergence"),
    ui.p(
                """
                This web application allows users to play with parameters of two models of evolutionary emergence.
                The first addresses how mutation at two loci, selection, and recombination between genotypes
                (00, 01, 10, and 11) shape the evolutionary dynamics of a pathogen on route to emergence, and how
                the number of infected individuals (incidence) aqnd 
                genotype frequencies depend on the life history of the resevoir host. The second model
                is due to Antia et al. (2003), in which branching processes are used to study how the probability that
                a pathogen adapts to a novel host depends on its initial degree of maladaptation, fitness landscape,
                mutation rate, and initial size of introduction. These models differ in key ways: the first addresses
                processes that are unfolding within the resevoir, and the second processes that are unfolding
                after the pathogen has sucessfully spilled over into a novel host (e.g., humans). In the first tab,
                users can play with parameters to solve the underlying system of equations to steady-state, or to plot
                variables over time. In the second tab, users can input parameters to visualize how the probability
                of emergence from a single introduction of a maladapted strain depends on those parameters,
                and do full simulations of the branching process from any introduction size."""
            ),
    ui.navset_tab_card(

    ui.nav('Spillcombinavoir!',
      ui.p(
                """
                """
            ),
      ui.navset_tab_card(
        ui.nav("Steady-state solutions",
        ui.p(
                """
                """
            ),
        ui.layout_sidebar(
        ui.panel_sidebar(
            ui.row(
            ui.column(3, ui.input_numeric("n", "Initial Population Size", value=100)),
            ui.column(3, ui.input_numeric("n_eq", "Equilibrium Population Size", value=250)),
            # ui.column(3, ui.input_numeric("i_0", "Initial Infected Individuals",value=1)),
            ),
            ui.markdown("""Select the parameters for plotting"""),
            ui.input_selectize("x1", "Parameter 1", ['c1','c2','log(mu1)','log(mu2)','log(nu1)', 'log(nu2)', 'sigma','turnover','Gamma','Delta'], selected = 'c1'),
            ui.output_ui("ui_select_p1"),
            ui.input_selectize("x2", "Parameter 2", ['c1','c2','log(mu1)','log(mu2)','log(nu1)', 'log(nu2)','sigma','turnover','Damma','Delta'], selected ='c2'),
            ui.output_ui("ui_select_p2"),
            ui.input_selectize("x5", "Dependent variable", ['I','IG11','IG10','IG00','IG01','p1','p2','x10','x11','x01','x00','D'],multiple = True),
            ui.input_action_button("run", "Run simulation", class_="btn-primary w-100")
          ),
      ui.panel_main(
        ui.output_plot("plot" ,width="800px", height="1000px"))
        # ui.panel_conditional(
        #   "input.run >= 1", 
        #   ui.output_plot("plot" ,width="800px", height="1000px"))
      )),
        ui.nav("Solutions through time",
        ui.p(
                """
                """
            ),
        ui.layout_sidebar(
        ui.panel_sidebar(
          ui.row(
            ui.column(3, ui.input_numeric("n_f", "Initial Population Size", value=100)),
            ui.column(3, ui.input_numeric("n_eq_f", "Equilibrium Population Size", value=250)),
            # ui.column(3, ui.input_numeric("i_0_f", "Initial Infected Individuals",value=1)),
            ui.column(3, ui.input_numeric("t_max_f", "Maximum Time (in days)", value=10)),
            ),
            ui.markdown("""Select the parameters for plotting"""),
            ui.input_selectize("x1_f", "Parameter 1", ['c1','c2','log(mu1)','log(mu2)','log(nu1)', 'log(nu2)','sigma','turnover','Gamma','Delta'], selected = 'c1'),
            # ui.output_ui("ui_select_p1_f"),
            ui.input_text("x1_f_val", "Parameter values (separate by ',')"),
            ui.input_selectize("x5_f", "Dependent variable", ['I','IG11','IG10','IG00','IG01','p1','p2','x10','x11','x01','x00','D'],multiple = True),
            ui.input_action_button(
                        "run_f", "Run simulation", class_="btn-primary w-100")
          ), 
        ui.panel_main(
          ui.output_plot("full_plot" ,width="800px", height="1000px")
          # ui.panel_conditional(
          # "input.run_f >= 1", 
          # ui.output_plot("full_plot" ,width="800px", height="1000px")
          # )
        ),
      ))
      )),
    ui.nav("Adaptation to Novel Host - Antia et al. (2003)",
      ui.p(
                """
                """
            ),
      ui.navset_tab_card(
      ui.nav('Emergence probability as function of initial degree of maladapation (R0)',
      ui.p(
                """
                """
            ),
      ui.layout_sidebar(
      ui.panel_sidebar(
            ui.row(
            ui.column(5, ui.input_text("mu", "Mutation rate (separate by ',')")),
            ui.column(5, ui.input_numeric("m", "Number of Types", value=5)),
            ),
            ui.markdown("""Select the R0 mode"""),
            ui.input_radio_buttons("r0_mode", "R0 Mode", ["Jackpot", "Additive", "Fitness Valley"], inline=True),
            ui.panel_conditional("input.r0_mode === Fitness Valley",
                                 ui.input_numeric("v", "Valley Depth", value=0.7)),
            ui.row(
            ui.column(5, ui.input_slider("r0_wild", "Initial R0 Range",value=(0.3, 1.1), min=0.0, max=2.0),
            ui.column(5, ui.input_numeric("r0_m", "Final R0", value=2)),
            ),
            ui.input_action_button("run_a", "Run simulation", class_="btn-primary w-100")
            )),
        ui.panel_main(
          ui.output_plot("plot_a" ,width="500px", height="500px"))
        # ui.panel_conditional(
        #   "input.run_a >= 1", 
        #   ui.output_plot("plot_a" ,width="500px", height="500px"))
      )
      ),

      ui.nav('Emergence probability as a function of R0 and mutation rate (mu)',
        ui.p(
                """
                """
            ),
         ui.layout_sidebar(
          ui.panel_sidebar(
            ui.row(
            ui.column(5, ui.input_slider("mu_1", "Mutation rate",value = (0.001,0.1),  min=0, max=0.1)),
            ui.column(5, ui.input_numeric("m_1", "Number of Types", value=5)),
            ),
            ui.markdown("""Select the R0 mode"""),
            ui.input_radio_buttons("r0_mode_1", "R0 Mode", ["Jackpot", "Additive", "Fitness Valley"], inline=True),
            ui.panel_conditional("input.r0_mode_1 === Fitness Valley",
                                 ui.input_numeric("v_1", "Valley Depth", value=0.7)),
            ui.row(
            ui.column(5, ui.input_slider("r0_wild_1", "Initial R0 Range",value=(0.3, 1.1), min=0.0, max=2.0),
            ui.column(5, ui.input_numeric("r0_m_1", "Final R0", value=2)),
            ),
            ui.input_action_button("run_a_1", "Run simulation", class_="btn-primary w-100")
            )),
        ui.panel_main(
          ui.output_plot("plot_a_heatmap" ,width="500px", height="500px")
        # ui.panel_conditional(
        #   "input.run_a_1 >= 1", 
        #   ui.output_plot("plot_a_heatmap" ,width="500px", height="500px"))
      )
      )
        ),
      ui.nav('Full simulation: branching process',
        ui.p(
                """
                """
            ),
         ui.layout_sidebar(
          ui.panel_sidebar(
            ui.row(
            ui.column(5, ui.input_text("mu_bp", "Mutation rate (separate by ',')")), ### string split by ','
            ui.column(5, ui.input_numeric("m_bp", "Number of Types", value=2)),
            ),
            ui.row(
            ui.column(5, ui.input_numeric("nindv_bp", "Number of Individuals",value = 10)),
            ui.column(5, ui.input_selectize("nsim_bp", "Number of Simulations",[1, 100,1000], selected=1)),
            ui.column(5, ui.input_numeric("ngen_bp", "Number of Generations", value=100)),
            ),
            ui.markdown("""Select the R0 mode"""),
            ui.input_radio_buttons("r0_mode_bp", "R0 Mode", ["Jackpot", "Additive", "Fitness Valley"], inline=True),
            ui.panel_conditional("input.r0_mode === Fitness Valley",
                                 ui.input_numeric("v_bp", "Valley Depth", value=0.7)),
            ui.row(
            ui.column(5, ui.input_text("r0_wild_bp", "Initial R0 (separate by ',')"), ### string split by ','
            ui.column(5, ui.input_numeric("r0_m_bp", "Final R0", value=2)),
            ),
            ui.input_action_button("run_a_bp", "Run simulation", class_="btn-primary w-100")
            )),
        ui.panel_main(
          ui.output_plot("plot_a_bp" ,width="900px", height="1500px")
        # ui.panel_conditional(
        #   "input.run_a_bp >= 1", 
        #   ui.output_plot("plot_a_bp" ,width="900px", height="1500px"))
        )
      ))))))



def server(input, output, session):

    @reactive.event(input.run)
    def steady_state_simu():
        N = input.n()
        # I0 = input.i_0()
        p1, p1_name= input.x3(),input.x1().lower()
        p2, p2_name = input.x4(),input.x2().lower()
        if p1_name == 'log(mu1)' or p1_name == 'log(nu1)'or p1_name =='log(mu2)' or p1_name =='log(nu2)':
          p1_name = p1_name[4:7]
          p1 = (10**(p1[0]), 10**(p1[1]))
        if p2_name == 'log(mu1)' or p2_name == 'log(nu1)'or p2_name =='log(mu2)' or p2_name =='log(nu2)':
          p2_name = p2_name[4:7]
          p2 = (10**(p2[0]), 10**(p2[1]))
        steady_state_data, ind_var = dr.steady_state(N, p1,p1_name,p2,p2_name,365*10)
        ind_var_df = pd.DataFrame(ind_var, columns=[p1_name, p2_name])
        return steady_state_data, ind_var_df
    @output
    @render.plot()    
    def plot() -> object:
        steady_state_data, df = steady_state_simu()
        n_inputs = len(input.x5())
        fig,ax = plt.subplots(n_inputs, 1)
        temp_x1, temp_x2 = input.x1().lower(), input.x2().lower()
        if input.x1() == 'log(mu1)' or input.x1()  == 'log(nu1)'or input.x1()  =='log(mu2)' or input.x1()  =='log(nu2)':
          temp_x1 = input.x1()[4:7]
        if input.x2() == 'log(mu1)' or input.x2()  == 'log(nu1)'or input.x2()  =='log(mu2)' or input.x2()  =='log(nu2)':
          temp_x2  = input.x2()[4:7] 
        for i, val in enumerate(input.x5()):
          dep_var = steady_state_data[val]
          complete_df = pd.concat([df, dep_var], axis=1)
          pivot_df = complete_df.pivot_table(index=temp_x1, columns=temp_x2, values=val)

          # Extract x and y values from the pivot table
          x_values = pivot_df.index
          y_values = pivot_df.columns

          # Extract z values as a 2D array from the pivot table
          z_values = pivot_df.values
          if n_inputs >1:
            temp = ax[i].pcolormesh(y_values, x_values, z_values, cmap='viridis')
            fig.colorbar(temp, label=val, ax = ax[i])
            ax[i].set_ylabel(temp_x1)
            ax[i].set_xlabel(temp_x2)
            ax[i].set_title('Plot of {} values by {} and {}'.format(val, temp_x1, temp_x2))
          else:
            temp = ax.pcolormesh(y_values, x_values, z_values, cmap='viridis')
            fig.colorbar(temp, label=val, ax = ax)
            ax.set_ylabel(temp_x1)
            ax.set_xlabel(temp_x2)
            ax.set_title('Plot of {} values by {} and {}'.format(val, temp_x1, temp_x2))

        return fig

    @reactive.event(input.run_f)
    def steady_state_full():
        N = input.n_f()
        # I0 = input.i_0_f()
        t_max = input.t_max_f()
        all_steady_state_data = {}
        p1, p1_name = [float(i) for i in input.x1_f_val().split(',')], input.x1_f().lower()
        for p in p1:
          if input.x1_f().lower() == 'log(mu1)' or input.x1_f().lower() == 'log(nu1)'or input.x1_f().lower() =='log(mu2)' or input.x1_f().lower() =='log(nu2)':
            p1_name = input.x1_f().lower()[4:7]
            p = 10**(p)
          print(p)
          steady_state_data, ind_var = dr.simu_by_time(N, p, p1_name,t_max)
          # ind_var_df = pd.DataFrame(ind_var, columns=[p1_name])
          all_steady_state_data[p] = steady_state_data
        # print(all_steady_state_data)
        return all_steady_state_data

    @output
    @render.plot()    
    def full_plot() -> object:
        all_steady_state_data = steady_state_full()
        n_inputs = len(input.x5_f())
        fig,ax = plt.subplots(n_inputs, 1)
        temp_x1_f = input.x1_f().lower()
        if input.x1_f() == 'log(mu1)' or input.x1_f()  == 'log(nu1)'or input.x1_f()  =='log(mu2)' or input.x1_f()  =='log(nu2)':
          temp_x1_f  = input.x1_f()[4:7]
        for p in all_steady_state_data:
          for i, val in enumerate(input.x5_f()):
            dep_var = all_steady_state_data[p][val]
            complete_df = pd.concat([all_steady_state_data[p]['time'],
                                    all_steady_state_data[p][temp_x1_f],
                                    dep_var], axis=1)
            # pivot_df = complete_df.pivot_table(index='time', columns=temp_x1_f, values=val)

            # # Extract x and y values from the pivot table
            # x_values = pivot_df.index
            # y_values = pivot_df.columns
            # # Extract z values as a 2D array from the pivot table
            # z_values = pivot_df.values
            # print(x_values[-1],y_values[-1], z_values[-1,-1])

            if n_inputs >1:
              temp = ax[i].plot(all_steady_state_data[p]['time'], dep_var, label = f'{temp_x1_f} = {p}')
              # fig.colorbar(temp, label=val, ax = ax[i])
              ax[i].set_xlabel('Time')
              ax[i].set_ylabel(val)
              ax[i].legend(loc = 'lower right')
              ax[i].set_title('Plot of {} values by {} and {}'.format(val, 'time',temp_x1_f))
              if val == 'D':
                ax[i].set_yscale('log')
            else:
              temp = ax.plot(all_steady_state_data[p]['time'], dep_var, label = f'{temp_x1_f} = {p}')
              # fig.colorbar(temp, label=val, ax = ax)
              ax.set_xlabel('Time')
              ax.set_ylabel(val)
              ax.legend(loc = 'lower right')
              ax.set_title('Plot of {} values by {} and {}'.format(val, 'time',temp_x1_f))
              if val == 'D':
                ax.set_yscale('log')
        return fig

    @reactive.event(input.run_a)
    def calc_all_prob():
      mu_dict = {}
      for i in [float(i) for i in input.mu().split(',')]:
        mu_dict[i] = [] 
        r0_list = np.linspace(input.r0_wild()[0],input.r0_wild()[1],50)
        for j in r0_list:
          em_prob = an.calculate_emergence_prob(input.m(),j, input.r0_m(),float(i), input.r0_mode(),input.v())
          mu_dict[i].append(em_prob)
      return mu_dict

    @output
    @render.plot()
    def plot_a():
      mu_dict = calc_all_prob()
      fig,ax = plt.subplots()
      for key, value in mu_dict.items():
        value = np.array(value)
        value[value <= 10e-9] = 0
        ax.plot(np.linspace(input.r0_wild()[0],input.r0_wild()[1],50), value, label=f'$\mu$={key}')
      ax.set_xlabel('R0')
      ax.set_ylabel('Emergence Probability')
      ax.legend(bbox_to_anchor=(1, 0.5))
      ax.set_yscale('log')
      return fig

    @reactive.event(input.run_a_1)
    def calc_all_prob_1():
      mu_dict = {}
      for i in np.linspace(input.mu_1()[0],input.mu_1()[1],50):
        mu_dict[i] = [] 
        r0_list = np.linspace(input.r0_wild_1()[0],input.r0_wild_1()[1],50)
        for j in r0_list:
          em_prob = an.calculate_emergence_prob(input.m_1(),j, input.r0_m_1(),float(i), input.r0_mode_1(),input.v_1())
          mu_dict[i].append(em_prob)
      return r0_list, mu_dict
    @output
    @render.plot()
    def plot_a_heatmap():
      r0_list, mu_dict = calc_all_prob_1()
      labels = list(mu_dict.keys())
      values = list(mu_dict.values())
      heatmap_data = np.array(values)
      fig,ax = plt.subplots()
      temp = ax.pcolormesh(r0_list, labels, heatmap_data, cmap='viridis',shading='auto',norm=colors.LogNorm(vmin=10e-7, vmax=1.))
      fig.colorbar(temp, ax = ax)
      ax.set_xlabel('R0')
      ax.set_ylabel('Mutation Rate')
      return fig 

    @reactive.event(input.run_a_bp)
    def simulate_bp():
      mu_list, wild_r0_list = [float(i) for i in input.mu_bp().split(',')] , [float(i) for i in input.r0_wild_bp().split(',')]
      m_type, nindv, nsim, ngen, r0_mode, r0_final = input.m_bp(),input.nindv_bp(),int(input.nsim_bp()),input.ngen_bp(),input.r0_mode_bp(),input.r0_m_bp()
      all_combos = an.mu_r0_combo(mu_list, wild_r0_list)
      all_combos_result = an.simulate_all_combos(m_type,mu_list, wild_r0_list, r0_final,r0_mode, nindv,ngen, nsim)
      return all_combos_result, (len(mu_list), len(wild_r0_list))
    @output
    @render.plot()
    def plot_a_bp():
      all_results, combo = simulate_bp()
      nsim = int(input.nsim_bp())
      rng = np.random.default_rng()
      mu_list, wild_r0_list = [float(i) for i in input.mu_bp().split(',')] , [float(i) for i in input.r0_wild_bp().split(',')]

      def create_subtitle(fig: plt.Figure, grid: SubplotSpec, title: str):
        "Sign sets of subplots with title"
        row = fig.add_subplot(grid)
        # the '\n' is important
        row.set_title(f'{title}\n', fontweight='semibold')
        # hide subplot
        row.set_frame_on(False)
        row.axis('off')

      rows = len(mu_list)*3
      cols = len(wild_r0_list)
      fig, axs = plt.subplots(rows, cols, figsize=(9, 9))
      for idx, (outer_key, inner_dict) in enumerate(all_results.items()):
        for inner_key, values in inner_dict.items():
          if isinstance(inner_key, int):
            print(values)
            print(len(values))
            mask=~np.isnan(values)
            maskedValues=values[mask]
            axs.flat[idx].plot(maskedValues, label=f'Type {inner_key+1}')
            axs.flat[idx].set_title(f'$\mu$={outer_key[0]}, $R_0$={outer_key[1]}')
            axs.flat[idx].legend(loc = 'upper right')
          elif inner_key == 'type_m':
            for v in values:
              if values.index(v) in rng.integers(nsim, size=5):
                al = 1
              else:
                al = 0.2 
              axs.flat[idx].plot(v, alpha = al)
            axs.flat[idx].set_title(f'$\mu$={outer_key[0]}, $R_0$={outer_key[1]}')
          elif inner_key == 'T_to_E':
            axs.flat[(len(mu_list)*len(wild_r0_list))+(idx)].hist(values, bins = 20)
            axs.flat[(len(mu_list)*len(wild_r0_list))+(idx)].set_title(f'$\mu$={outer_key[0]}, $R_0$={outer_key[1]}')
          elif inner_key == 'Extinct':
            axs.flat[(len(mu_list)*len(wild_r0_list))*2+(idx)].hist(values, bins = 20)
            axs.flat[(len(mu_list)*len(wild_r0_list))*2+(idx)].set_title(f'$\mu$={outer_key[0]}, $R_0$={outer_key[1]}')
      # fig.tight_layout()
      grid = plt.GridSpec(rows, cols)
      if nsim == 1:
        create_subtitle(fig, grid[0, ::], 'Average Number of Individuals of Each Type')
      else:
        create_subtitle(fig, grid[0, ::], f'Number of Individuals of Type {input.m_bp()}')
      create_subtitle(fig, grid[len(mu_list), ::], 'Time to Final Mutation')
      create_subtitle(fig, grid[len(mu_list)*2, ::], 'Time to Extinction')
      fig.tight_layout()
      return fig

      


    @output
    @render.ui
    def ui_select_p1():
      x = input.x1()
      if x == 'Gamma' or x == 'Delta':
        v_min, v_max = 0.0, 50.0
        v_range = (25,27)
      elif x == 'c1' or x == 'c2':
        v_min, v_max = 0.0, 1.0
        v_range = (0.5,0.7)
      elif x == 'turnover':
        v_min, v_max = 0.0001, 0.01
        v_range = (0.005,0.008)
      elif x == 'sigma':
        v_min, v_max = 0, 0.1
        v_range = (0.03, 0.05)
      elif x == 'log(mu1)' or x == 'log(nu1)'or x =='log(mu2)' or x=='log(nu2)': #'log(mu1)','log(mu2)','log(nu1)', 'log(nu2)'
        v_min, v_max = -8.0, -1.0
        v_range = (-5,-3)
      return ui.input_slider("x3", "Range for {}".format(x), value=v_range, min=v_min, max=v_max)
    @output
    @render.ui
    def ui_select_p2():
      x = input.x2()
      '''
      Gamma, Delta (capital): 0-50 with a defaults of 25. 
      R0: 2.5-10 with a default of 5. 
      Costs: 0-1 with defaults of 0.5. 
      Turnover: 0.0001-0.01 with a default of 0.01. 
      sigma: 0-1 with a default of 0.1. 
      Mutation rates: 10^-8-10^-2 with default of 10^-5.
      '''
      if x == 'Gamma' or x == 'Delta':
        v_min, v_max = 0.0, 50.0
        v_range = (25,27)
      elif x == 'c1' or x == 'c2':
        v_min, v_max = 0.0, 1.0
        v_range = (0.5,0.7)
      elif x == 'turnover':
        v_min, v_max = 0.0001, 0.01
        v_range = (0.005,0.008)
      elif x == 'sigma':
        v_min, v_max = 0, 0.1
        v_range = (0.03, 0.05)
      elif x == 'log(mu1)' or x == 'log(nu1)'or x =='log(mu2)' or x=='log(nu2)':
        v_min, v_max = -8.0, -1.0
        v_range = (-5,-3)
      return ui.input_slider("x4", "Range for {}".format(x), value=v_range, min=v_min, max=v_max)
    @output
    @render.ui
    def ui_select_p1_f():
      x = input.x1_f()
      if x == 'Gamma' or x == 'Delta':
        v_min, v_max = 0.0, 50.0
        v_range = (25,27)
      elif x == 'c1' or x == 'c2':
        v_min, v_max = 0.0, 1.0
        v_range = (0.5,0.7)
      elif x == 'turnover':
        v_min, v_max = 0.0001, 0.01
        v_range = (0.005,0.008)
      elif x == 'sigma':
        v_min, v_max = 0, 0.1
        v_range = (0.03, 0.05)
      elif x == 'log(mu1)' or x == 'log(nu1)'or x =='log(mu2)' or x=='log(nu2)':
        v_min, v_max = -8.0, -1.0
        v_range = (-5,-3)
      return ui.input_slider("x4_f", "Range for {}".format(x), value=v_range, min=v_min, max=v_max)


app = App(app_ui, server)
