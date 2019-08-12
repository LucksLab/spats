import matplotlib.pyplot as plt

from nb import colors, cotrans_matrix, cotrans_matrix_data, spats_run_data

def plot_sl_counts(path = None):
    run_data = spats_run_data(path)
    row = run_data.single_profile
    plt.xlim([0, run_data.n + 1])
    plt.plot(row.x_axis, row.treated, color = colors.red, label = 'f+')
    plt.plot(row.x_axis, row.untreated, color = colors.blue, label = 'f-')
    plt.title("Total Treated/Untreated Counts")
    plt.xlabel("Site")
    plt.ylabel("# of Stops")
    plt.legend()
    return plt

def plot_reactivity(reactivity_type):
    run_data = spats_run_data()
    row = run_data.single_profile
    reactivity = getattr(row, reactivity_type)
    plt.xlim([0, run_data.n + 1])
    plt.ylim([0, max(reactivity)])
    plt.plot(row.x_axis, reactivity, 1, label=reactivity_type)
    #plt.bar(row.x_axis, reactivity, 1, label=reactivity_type)
    plt.xlabel("Nucleotide Position (nt)")
    plt.ylabel(reactivity_type)
    plt.legend()
    return plt

def plot_sl_rho():
    return plot_reactivity('rho')

def plot_sl_beta():
    return plot_reactivity('beta')

def plot_sl_theta():
    return plot_reactivity('theta')

def plot_sl_muts():
    run_data = spats_run_data()
    row = run_data.single_profile
    plt.xlim([0, run_data.n + 1])
    #plt.plot(row.x_axis, row.treated_counts, color = "red", label = 's+')
    #plt.plot(row.x_axis, row.untreated_counts, color = "blue", label = 's-')
    plt.plot(row.x_axis, row.treated_muts, color = "orange", label = 'mut+')
    plt.plot(row.x_axis, row.untreated_muts, color = "purple", label = 'mut-')
    plt.title("Total Treated/Untreated Mutations")
    plt.xlabel("Site")
    plt.ylabel("# of Mutations")
    plt.legend()
    return plt

def plot_sl_muts_reactivity():
    return plot_reactivity('r_mut')

def plot_cotrans_counts():
    run_data = spats_run_data()
    plt.xlim([run_data.min_length, run_data.n + 1]) # Change x-axis here
    plt.plot(run_data.all_sites, run_data.total_treated_counts, color = colors.red, label = '(+)')
    plt.plot(run_data.all_sites, run_data.total_untreated_counts, color = colors.blue, label = '(-)')
    plt.title("Total Treated/Untreated Counts")
    plt.xlabel("RNA Length")
    plt.ylabel("# of Stops")
    plt.legend()
    return plt

def plot_cotrans_c():
    run_data = spats_run_data()
    plt.xlim([run_data.min_length, run_data.n + 1]) #Change x-axis here
    plt.plot(run_data.all_sites, run_data.c_values, color = colors.black, label = "c")
    plt.plot(run_data.all_sites, [0.4 for i in run_data.all_sites], color = colors.red, label = "Recommended Cutoff")
    ax = plt.gca(); ax.yaxis.grid(True)
    plt.title("c Values")
    plt.xlabel("RNA Length")
    plt.ylabel("c")
    plt.legend()
    return plt

def plot_cotrans_treated():
    matrix_data = cotrans_matrix_data('treated_counts', max_val = 5000)
    plt.matshow(matrix_data,cmap='jet')
    ax = plt.gca()
    ax.grid(color='grey',linestyle='-',linewidth='0.5')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel("Nucleotide (nt)")
    plt.ylabel("RNA Length (nt)")
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Counts')
    plt.title('Treated Channel Counts')
    return plt

def plot_cotrans_untreated():
    matrix_data = cotrans_matrix_data('untreated_counts', max_val = 5000)
    plt.matshow(matrix_data,cmap='jet')
    ax = plt.gca()
    ax.grid(color='grey',linestyle='-',linewidth='0.5')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel("Nucleotide (nt)")
    plt.ylabel("RNA Length (nt)")
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Counts')
    plt.title('Untreated Channel Counts')
    return plt

def plot_cotrans_treated_muts():
    matrix_data = cotrans_matrix_data('treated_muts', max_val = 5000)
    plt.matshow(matrix_data,cmap='jet')
    ax = plt.gca()
    ax.grid(color='grey',linestyle='-',linewidth='0.5')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel("Nucleotide (nt)")
    plt.ylabel("RNA Length (nt)")
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Counts')
    plt.title('Treated Channel Mutations')
    return plt

def plot_cotrans_untreated_muts():
    matrix_data = cotrans_matrix_data('untreated_muts', max_val = 5000)
    plt.matshow(matrix_data,cmap='jet')
    ax = plt.gca()
    ax.grid(color='grey',linestyle='-',linewidth='0.5')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel("Nucleotide (nt)")
    plt.ylabel("RNA Length (nt)")
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Counts')
    plt.title('Untreated Channel Mutations')
    return plt

def plot_matrix(reactivity_type, max_val = 0.025):
    matrix_data = cotrans_matrix_data(reactivity_type, max_val = max_val)
    plt.matshow(matrix_data)
    ax = plt.gca()
    ax.grid(color='grey',linestyle='-',linewidth='0.5')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel("Nucleotide (nt)")
    plt.ylabel("RNA Length (nt)")
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label(reactivity_type)
    return plt

def plot_cotrans_rho():
    return plot_matrix('rho', 4.0)

def plot_cotrans_beta():
    return plot_matrix('beta', 0.025)

def plot_cotrans_theta():
    return plot_matrix('theta', 0.025)

def plot_matrix_muts(reactivity_type, max_val = 0.025):
    matrix_data = cotrans_matrix_data(reactivity_type, max_val = max_val)
    plt.matshow(matrix_data)
    ax = plt.gca()
    ax.grid(color='grey',linestyle='-',linewidth='0.5')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel("Nucleotide (nt)")
    plt.ylabel("RNA Length (nt)")
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label(reactivity_type)
    return plt

def plot_cotrans_r_mut():
    return plot_matrix_muts('r', 0.025)

def plot_cotrans_mu_mut():
    return plot_matrix_muts('mu', 0.025)

def plot_cotrans_beta_mut():
    return plot_matrix_muts('beta', 0.025)

def plot_cotrans_treated_mut():
    return plot_matrix_muts('treated_mut', 0)

def plot_cotrans_untreated_mut():
    return plot_matrix_muts('untreated_mut', 0)

