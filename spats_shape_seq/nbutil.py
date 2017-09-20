
import base64
import datetime
import IPython.display as Idp
import nbformat as nbf
import os

from IPython.utils.py3compat import str_to_bytes, bytes_to_str


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
def create_code_cell(code, execute = True, where = 'below'):
    encoded_code = bytes_to_str(base64.b64encode(str_to_bytes(code)))
    Idp.display(Idp.Javascript(_insert_code_cell_template.format(where, encoded_code, "code.execute();" if execute else "")))

def create_html_cell(html):
    Idp.display(Idp.HTML(html))

def create_json_cell(json_dict):
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
        nb = self.add_md_cell(spats_run_md_template.format(self._stamp()))
        if cotrans:
            nb.add_code_cell(cotrans_counts_template)
            nb.add_code_cell(cotrans_c_value_template)
        else:
            pass
        return nb



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
plt.plot(pre_data.x_axis, pre_data.treated, color = colors.green)
plt.plot(pre_data.x_axis, pre_data.untreated, color = colors.black)
plt.plot(pre_data.x_axis, pre_data.base, color = colors.cyan)
plt.ylim([0, 30000])
plt.gcf().set_size_inches(36, 8)
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
plt.plot(run_data.all_sites, run_data.total_treated_counts, color = colors.red)
plt.plot(run_data.all_sites, run_data.total_untreated_counts, color = colors.blue)
plt.title("Total Treated/Untreated Counts")
plt.legend(["f+", "f-"])
plt.xlabel("Length")
plt.ylabel("# of Stops")
plt.xlim([run_data.min_length, run_data.n + 1])
plt.gcf().set_size_inches(36, 8)
plt.show()
"""

cotrans_c_value_template = """
run_data = spats_run_data()
plt.plot(run_data.all_sites, run_data.c_values, color = colors.black)
plt.title("c Values")
plt.xlabel("Length")
plt.ylabel("c")
plt.xlim([run_data.min_length, run_data.n + 1])
plt.gcf().set_size_inches(36, 8)
plt.show()
"""
