
define([
    'jquery',
    'base/js/namespace',
    'base/js/events'
], function(
    $,
    Jupyter,
    events
) {
    "use strict";

    var update_input_visibility = function () {
        Jupyter.notebook.execute_all_cells();
        Jupyter.notebook.get_cells().forEach(function(cell) {
            if (cell.input) {
                var codeCell = cell.element.find("div.inner_cell")[0];
                console.log("Found element", cell, codeCell);
                codeCell.style.display = 'none';
                cell.element.find("div.input_prompt").click(function(){
                    codeCell.style.display = (codeCell.style.display == 'block' ? 'none' : 'block');
                });
            }
        })
    };

    var run_all_cells = function() {
        console.log("Auto-running all cells-below...");
        Jupyter.actions.call('jupyter-notebook:run-all-cells');
        Jupyter.actions.call('jupyter-notebook:save-notebook');
    }

    var load_ipython_extension = function() {
        if (Jupyter.notebook !== undefined  &&  Jupyter.notebook._fully_loaded)
            update_input_visibility();
        events.on("notebook_loaded.Notebook", update_input_visibility);
        events.on("kernel_ready.Kernel", run_all_cells);
    };

    return {
        load_ipython_extension : load_ipython_extension
    };
});
