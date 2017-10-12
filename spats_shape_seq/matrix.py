_template_css = """
.matrix_table {
overflow: hidden;
background-color: #ccc;
}

.matrix_table td {
width: 2px; height: 3px;
}

.matrix_table td:not(.no_site):hover {
background-color: green;
}

.matrix_table tr {
background-color: #fff;
}

.matrix_table tr:hover {
  background-color: #ffa;
  border-style: solid;
  border-color: #fff;
  border-width: 1px;
}

.matrix_table td {
  position: relative;
}

.matrix_table td:not(.no_site):hover::after {
  content: "";
  position: absolute;
  background-color: #ffa;
  border-style: solid;
  border-color: #fff;
  border-width: 1px;
  left: 0;
  top: -5000px;
  height: 10000px;
  width: 100%;
  z-index: 1;
}

.no_site {
  background-color: #ccc;
  z-index: 2;
}
"""

_base_template = """
<div>
 <style>{css}</style>
 <div class="matrix_div">
  <div class="table_div">
   <table class="matrix_table">
{table}
   </table>
  </div>
 </div>
</div>
"""

_row_template = '<tr>{}</tr>\n'

def matrix_html(min_length, n, profiles):
    return _base_template.format(css = _template_css, table = _make_table(min_length, n))

def _make_table(min_length, n):
    res = ""
    for end in range(min_length, n + 1):
        res += _row_template.format(('<td/>' * (end + 1)) + '<td class="no_site" colspan="{}"/>'.format(n + 1 - end))
    return res
