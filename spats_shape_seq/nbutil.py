
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
        return datetime.datetime.now().strftime('%Y/%m/%d' if dateOnly else '%Y/%m/%d %H:%M')

    def is_empty(self):
        return not(bool(self._nb.cells))

    def add_code_cell(self, code):
        self._nb.cells.append(nbf.v4.new_code_cell(code.strip()))
        return self

    def add_md_cell(self, md):
        self._nb.cells.append(nbf.v4.new_markdown_cell(md.strip()))
        return self

    def add_metadata(self, metadata):
        metadata = { k : metadata.get(k, "") for k in _metadata_template_keys }
        metadata['date'] = metadata.get('date', self._stamp(True))
        return self.add_md_cell(metadata_template.format(**metadata))

    def add_initializer(self):
        return self.add_code_cell(initializer_code_template)

    def add_preseq(self):
        return self.add_md_cell(preseq_md_template.format(self._stamp())).add_code_cell(preseq_code_template)



##########################################################
## templates
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
%load_ext rpy2.ipython
from spats_shape_seq.nb import *
import matplotlib.pyplot as plt
%matplotlib inline
"""

preseq_md_template = """
Pre-Sequencing
--------------
{}
"""

preseq_code_template = """
pre_data = nb_load_pre_data()
pre_colors = [ [0,1,0], [0,0,0], [0,1,1] ]
x_values = range(len(pre_data[0]))
for i in range(len(pre_data)):
    plt.plot(x_values, pre_data[i], color = pre_colors[i % len(pre_colors)])
plt.ylim([0, 30000])
plt.gcf().set_size_inches(36, 8)
plt.show()
print "Treated (green) 1st moment: {}".format(first_moment(pre_data[0]))
print "Untreated (black) 1st moment: {}".format(first_moment(pre_data[1]))
"""
