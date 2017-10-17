
import base64
import datetime
import nbformat as nbf
import os


##########################################################
## methods for live-updating notebook
##   (should only be called from nbutil module)
##########################################################

# where = at_bottom, above, below (default)
_insert_code_cell_template = """
var code = IPython.notebook.insert_cell_{}('code');
code.set_text(atob('{}'));
{}
"""

def _load_ipython():
    import IPython.display as Idp
    from IPython.utils.py3compat import str_to_bytes, bytes_to_str

def create_code_cell(code, execute = True, where = 'below'):
    _load_ipython()
    encoded_code = bytes_to_str(base64.b64encode(str_to_bytes(code)))
    Idp.display(Idp.Javascript(_insert_code_cell_template.format(where, encoded_code, "code.execute();" if execute else "")))

def create_html_cell(html):
    _load_ipython()
    Idp.display(Idp.HTML(html))

def create_json_cell(json_dict):
    _load_ipython()
    Idp.display(Idp.JSON(json_dict))





##########################################################
## class for modifying on-disk notebook
##########################################################

_NBF_VERSION = 4

class Notebook(object):

    def __init__(self, path = None):
        self.path = path
        if path and os.path.exists(path):
            self.load(path)
        else:
            self._nb = nbf.v4.new_notebook(cells=[])

    def load(self, path):
        self._nb = nbf.read(open(path, 'r'), _NBF_VERSION)

    def save(self, path = None):
        nbf.write(self._nb, open(path or self.path, 'w'), _NBF_VERSION)

    def _stamp(self, dateOnly = False):
        return datetime.datetime.now().strftime('%Y/%m/%d' if dateOnly else '%Y/%m/%d %I:%M%p')

    def is_empty(self):
        return not(bool(self._nb.cells))

    def add_code_cell(self, code):
        self._nb.cells.append(nbf.v4.new_code_cell(code.strip()))
        return self

    def add_md_cell(self, md):
        self._nb.cells.append(nbf.v4.new_markdown_cell(md.strip()))
        return self

    def add_metadata(self, metadata):
        metadata['date'] = metadata.get('date') or self._stamp(True)
        filled_metadata = { k : metadata.get(k, "").strip() for k in _metadata_template_keys }
        self.add_md_cell(metadata_template.format(**filled_metadata))
        self._nb.metadata['spats_info'] = metadata
        return self

    def add_initializer(self):
        return self.add_code_cell(initializer_code_template)

    def add_preseq(self):
        return self.add_md_cell(preseq_md_template.format(self._stamp())).add_code_cell(preseq_code_template)

    def add_spats_run(self, cotrans):
        self.add_md_cell(spats_run_md_template.format(self._stamp()))
        if cotrans:
            self.add_code_cell(cotrans_counts_template)
            self.add_code_cell(cotrans_c_value_template)
            self.add_code_cell(cotrans_matrix_treated_template)
            self.add_code_cell(cotrans_matrix_untreated_template)
            self.add_md_cell(cotrans_reactivity_matrix_md_template)
            self.add_code_cell(cotrans_reactivity_matrix_template)
            self.add_md_cell(cotrans_reactivity_analysis_md_template)
            self.add_code_cell(cotrans_reactivity_analysis_template)
            self.add_code_cell(cotrans_reactivity_analysis_2_template)
        else:
            pass
        return self



##########################################################
## notebook templates
##########################################################


metadata_template = """
{name}
------

Date: {date}

Author: {author}

Folding Conditions:

* Buffer: {buffer}

* Temperature: {temp}

* Salt Concentration: {salt}

Probe: {probe}

Lot Numbers:

* Adapter: {adapter}

* Enzyme: {enzyme}

* Reagent: {reagent}
"""

_metadata_template_keys = [ 'name', 'date', 'author', 'buffer', 'temp', 'salt', 'probe', 'adapter', 'enzyme', 'reagent' ]

initializer_code_template = """
from spats_shape_seq.nb import *
import matplotlib.pyplot as plt
%matplotlib inline
%load_ext rpy2.ipython
"""

preseq_md_template = """
Pre-Sequencing
--------------
{}
"""

preseq_code_template = """
pre_data = preseq_data()
default_ylim = pre_data.max_val
default_xlim = max(pre_data.x_axis)
plt.ylim([0, default_ylim]) #Change y-axis here
plt.xlim([0, default_xlim]) #Change x-axis here
plt.style.use('fsa')
plt.plot(pre_data.x_axis, pre_data.treated, color = colors.green, label = '(+)')
plt.plot(pre_data.x_axis, pre_data.untreated, color = colors.black, label = '(-)')
plt.plot(pre_data.x_axis, pre_data.base, color = colors.cyan, label = 'Standard')
plt.xlabel('Elution Time (pixel)')
plt.ylabel('Intensity (AU)')
plt.gcf().set_size_inches(36, 8)
plt.legend()
plt.show()
print "Treated 1st moment: {}".format(first_moment(pre_data.treated))
print "Untreated 1st moment: {}".format(first_moment(pre_data.untreated))
"""

spats_run_md_template = """
SPATS Run
--------------
{}
"""

cotrans_counts_template = """
run_data = spats_run_data()
plt.style.use('fsa')
plt.xlim([run_data.min_length, run_data.n + 1]) # Change x-axis here
plt.plot(run_data.all_sites, run_data.total_treated_counts, color = colors.red, label = '(+)')
plt.plot(run_data.all_sites, run_data.total_untreated_counts, color = colors.blue, label = '(-)')
plt.title("Total Treated/Untreated Counts")
plt.xlabel("RNA Length")
plt.ylabel("# of Stops")
plt.legend()
plt.gcf().set_size_inches(36, 8)
plt.show()
"""

cotrans_c_value_template = """
run_data = spats_run_data()
plt.style.use('fsa')
plt.xlim([run_data.min_length, run_data.n + 1]) #Change x-axis here
plt.plot(run_data.all_sites, run_data.c_values, color = colors.black, label = "c")
plt.plot(run_data.all_sites, [0.4 for i in run_data.all_sites], color = colors.red, label = "Recommended Cutoff")
ax = plt.gca(); ax.yaxis.grid(True)
plt.title("c Values")
plt.xlabel("RNA Length")
plt.ylabel("c")
plt.legend()
plt.gcf().set_size_inches(36, 8)
plt.show()
"""

cotrans_matrix_treated_template = """
matrix_data = cotrans_matrix_data('treated_counts', max_val = 5000)
plt.style.use('fsa') #Todo - change to be a custom matrix style to remove large fonts
plt.matshow(matrix_data,cmap='jet')
ax = plt.gca()
ax.grid(color='grey',linestyle='-',linewidth='0.5')
ax.xaxis.set_ticks_position('bottom')
plt.xlabel("Nucleotide (nt)")
plt.ylabel("RNA Length (nt)")
cbar = plt.colorbar(fraction=0.046, pad=0.04)
cbar.set_label('Counts')
plt.title('Treated Channel Counts')
plt.gcf().set_size_inches(12, 12)
plt.show()
"""

cotrans_matrix_untreated_template = """
matrix_data = cotrans_matrix_data('untreated_counts', max_val = 5000)
plt.style.use('fsa') #Todo - change to be a custom matrix style to remove large fonts
plt.matshow(matrix_data,cmap='jet')
ax = plt.gca()
ax.grid(color='grey',linestyle='-',linewidth='0.5')
ax.xaxis.set_ticks_position('bottom')
plt.xlabel("Nucleotide (nt)")
plt.ylabel("RNA Length (nt)")
cbar = plt.colorbar(fraction=0.046, pad=0.04)
cbar.set_label('Counts')
plt.title('Untreated Channel Counts')
plt.gcf().set_size_inches(12, 12)
plt.show()
"""

cotrans_reactivity_matrix_md_template = """
Reactivity Matrix
------

Comments:

 - XXX
"""

cotrans_reactivity_matrix_template = """
reactivity_type = 'beta' #Can be 'rho', 'beta', 'theta'
matrix_data = cotrans_matrix_data(reactivity_type, max_val = 0.025) #Suggested max_val: beta = 0.025, rho = 4
#Todo - interface that allows to get subsets of the matrix from cotrans_matrix_data (to create subset plots)
plt.style.use('fsa') #Todo - change to be a custom matrix style to remove large fonts
plt.matshow(matrix_data)
ax = plt.gca()
ax.grid(color='grey',linestyle='-',linewidth='0.5')
ax.xaxis.set_ticks_position('bottom')
plt.xlabel("Nucleotide (nt)")
plt.ylabel("RNA Length (nt)")
cbar = plt.colorbar(fraction=0.046, pad=0.04)
cbar.set_label(reactivity_type)
plt.gcf().set_size_inches(12, 12)
plt.show()
"""

cotrans_matrix_template_html = """
cotrans_matrix(data_type = 'rho', max_val = 4.0, flags = False)
"""

cotrans_reactivity_analysis_md_template = """
Reactivity Analysis
------

Comments:

 - XXX
"""

cotrans_reactivity_analysis_template = """
run_data = spats_run_data()
reactivity_type = 'beta' #Can be 'rho', 'beta', 'theta'
end_lengths = [99,115]
row_data = [run_data.row(end_length) for end_length in end_lengths]
plt.style.use('fsa')

#Reactivity bar plots
xmax = 0
ymax = 0
bar_width = 1./float(len(end_lengths))
for i in range(len(row_data)):
    row = row_data[i]
    xmax = max(row.x_axis)
    ymax = max(ymax,max(row.beta))
    plt.bar([x + i*bar_width for x in row.x_axis], 
            row.beta, 
            bar_width,
            label=end_lengths[i]) #TODO - interface that allows y-data to be reatrieved by using the reactivity_type variable
#Todo - set xticks better with minor at each position
plt.xlim([0, xmax])
plt.ylim(0, ymax)
plt.title(reactivity_type+", length = {}".format(end_lengths))
plt.xlabel("Nucleotide Position (nt)")
plt.ylabel(reactivity_type)
plt.legend()
plt.gcf().set_size_inches(36, 8)
plt.show()

#Individual +/- fragment distribution plots
plt.style.use('fsa')
f, sub_axes = plt.subplots(len(row_data), sharex=True, sharey=True)
sub_axes[0].set_title("Treated/Untreated Fragment Distribution, length = {}".format(end_lengths))
for i in range(len(row_data)):
    row = row_data[i]
    sub_axes[i].plot(row.x_axis, normalize(row.treated), color = colors.red, label = 'f+')
    sub_axes[i].plot(row.x_axis, normalize(row.untreated), color = colors.blue, label = 'f-')
    sub_axes[i].set_ylabel("%")
plt.xlim([0, run_data.n + 1])
plt.legend()
f.set_size_inches(36, 2*len(row_data))
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) #mashing plots together
plt.xlabel("Nucleotide Position(nt)")
plt.show()
"""

cotrans_reactivity_analysis_2_template = """
run_data = spats_run_data()
reactivity_type = 'beta' #Can be 'rho', 'beta', 'theta'
sites = [22,30,45]
col_data = [run_data.column(site) for site in sites]
plt.style.use('fsa')
ymax = 0
for i in range(len(col_data)):
    col = col_data[i]
    xmax = max(col.x_axis)
    ymax = max(ymax,max(col.beta))
    plt.plot(col.x_axis, col.beta,label=sites[i]) #TODO - interface that allows y-data to be reatrieved by using the reactivity_type variable
plt.xlim([run_data.min_length, run_data.n + 1])
plt.title("Reactivity Traces (Beta) {} ".format(sites)) #Todo - change beta label based on reactivity_type
plt.xlabel("RNA Length (nt)")
plt.ylabel(reactivity_type)
plt.legend()
plt.gcf().set_size_inches(36, 8)
plt.show()
"""
