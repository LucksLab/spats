
import ast
import ConfigParser
import csv
import datetime
import json
import multiprocessing
import os
import shutil
import subprocess
import sys
import time

import spats_shape_seq
import spats_shape_seq.nbutil as nbutil
from spats_shape_seq import Spats
from spats_shape_seq.parse import abif_parse, fastq_handle_filter
from spats_shape_seq.reads import ReadsData, ReadsAnalyzer


class SpatsTool(object):

    def __init__(self, path):
        self.path = path or os.getcwd()
        self.config = None
        self.cotrans = False
        self._skip_log = False
        self._no_config_required_commands = [ "doc", "help", "init", "viz" ]
        self._private_commands = [ "viz" ]
        self._temp_files = []
        self._r1 = None
        self._r2 = None
        self._parse_config()

    def _parse_config(self):
        config_path = os.path.join(self.path, "spats.config")
        if not os.path.exists(config_path):
            return
        parser = ConfigParser.SafeConfigParser()
        parser.read(config_path)
        config = {}
        for section in parser.sections():
            if section not in config:
                config[section] = {}
            for name, value in parser.items(section):
                config[section][name] = value
        if config and 'spats' in config:
            self.config = config['spats']
            if 'cotrans' in self.config:
                cotrans = self.config['cotrans']
                self.cotrans = bool(ast.literal_eval(cotrans))
        self.metadata = config.get('metadata', {})

    def _module_path(self):
        return os.path.dirname(spats_shape_seq.__file__)

    def _spats_path(self):
        return os.path.normpath(os.path.join(self._module_path(), ".."))

    def _native_tool(self, tool_name):
        bin_path = os.path.join(self._spats_path(), "native", "bin", tool_name)
        return bin_path if os.path.exists(bin_path) else None

    def _add_note(self, note):
        self._notes.append(note)
        print(":{}".format(note))

    def _load_r1_r2(self, ret_r1):
        def decomp(rx):
            base, ext = os.path.splitext(os.path.basename(rx))
            if ext.lower() == '.gz':
                self._add_note("decompressing {}".format(rx))
                out = os.path.join(self.path, base + ".tmp")
                subprocess.check_call('gzip -d -c "{}" > "{}"'.format(rx, out), cwd = self.path, shell = True)
                self._temp_files.append(out)
                self._sentinel("decompress {}".format('R1' if ret_r1 else 'R2'))
                return out
            return rx
        self._r1, self._r2 = decomp(self.config['r1']), decomp(self.config['r2'])
        return self._r1 if ret_r1 else self._r2

    @property
    def r1(self):
        return self._r1 or self._load_r1_r2(True)

    @property
    def r2(self):
        return self._r2 or self._load_r1_r2(False)


    def _run(self, args):
        if not args:
            print("Command required. Try 'spats_tool help'.")
            return

        command = args[0]
        self._command_args = args[1:]
        if not self.config and command not in self._no_config_required_commands:
            print("Missing spats.config")
            return

        self._notes = []

        self.start = time.time()
        hdlr = getattr(self, command, None)
        if not hdlr:
            print("Invalid command: {}".format(command))
            return
        try:
            hdlr()
            failure = False
        except Exception, e:
            print("** Command {} failed. ({})".format(command, e))
            failure = True

        if not failure and not self._skip_log:
            delta = self._sentinel("{} complete".format(command))
            self._log(command, delta)

        for f in self._temp_files:
            if os.path.exists(f):
                os.remove(f)

    def _sentinel(self, label):
        delta = time.time() - self.start
        self._add_note("{} @ {:.2f}s".format(label, delta))
        return delta

    def _log(self, command, delta):
        stamp = datetime.datetime.now().strftime('%Y/%m/%d %H:%M')
        with open(os.path.join(self.path, 'spats.log'), 'ab') as outfile:
            outfile.write("{} : {}, {:.2f}s\n".format(stamp, command, delta))
            for note in self._notes:
                outfile.write("   - {}\n".format(note))
            outfile.write("\n")

    def init(self):
        """Set up a spats_tool folder.
        """
        self._skip_log = True
        config_path = os.path.join(self.path, "spats.config")
        if not os.path.exists(config_path):
            open(config_path, 'wb').write(_spats_config_template)
            self._add_note("Created default spats.config, please edit before running tools.")
        else:
            self._add_note("** spats.config already exists, not overwriting!")

    def reads(self):
        """Perform reads analysis on the r1/r2 fragment pairs, for use with the visualization tool.
        """

        db_name = self._reads_file()
        if os.path.exists(db_name):
            self._add_note("** removing previous reads file")
            os.remove(db_name)

        native_tool = self._native_tool('reads')
        if native_tool:
            self._add_note("using native reads")
            subprocess.check_call([native_tool, self.config['target'], self.r1, self.r2, db_name], cwd = self.path)

        data = ReadsData(db_name)
        if not native_tool:
            self._add_note("using python reads")
            data.parse(self.config['target'], self.r1, self.r2)

        analyzer = ReadsAnalyzer(data, cotrans = self.cotrans)
        self._update_run_config(analyzer.run)
        analyzer.process_tags()
        self._add_note("tags processed to {}".format(os.path.basename(db_name)))

    def pre(self):
        """Process the pre-sequence data file.
        """
        if 'preseq' not in self.config:
            raise Exception("Missing 'preseq' file in spats.config")
        pre_files = [ f.strip() for f in self.config['preseq'].split(',') ]
        for filename in pre_files:
            key, _ = os.path.splitext(os.path.basename(filename))
            pre_name = self._pre_file(key)
            if os.path.exists(pre_name):
                self._add_note("** removing previous preseq file")
                os.remove(pre_name)
            preseq_data = abif_parse(filename, fields = [ 'DATA2', 'DATA3', 'DATA105' ])
            open(pre_name, 'wb').write(json.dumps(preseq_data))
            self._add_note("pre-sequencing data processed to {}".format(os.path.basename(pre_name)))
            self._notebook().add_preseq(key).save()

    def _run_plots(self):
        self._notebook().add_spats_run(self.cotrans, True).save()

    def run(self):
        """Process the SPATS data for the configured target(s) and r1/r2 fragment pairs.
        """
        run_name = self._run_file()
        if os.path.exists(run_name):
            self._add_note("** removing previous run file")
            os.remove(run_name)

        native_tool = self._native_tool('cotrans')
        if native_tool and not self.cotrans:
            self._add_note("skipping native tool due to non-cotrans run")
            native_tool = None

        spats = Spats(cotrans = self.cotrans)
        if self._update_run_config(spats.run) and native_tool:
            self._add_note("skipping native tool due to custom config")
            native_tool = None

        if native_tool:
            self._add_note("using native cotrans processor")
            subprocess.check_call([native_tool, self.config['target'], self.r1, self.r2, run_name], cwd = self.path)
        else:
            self._add_note("using python processor")
            spats.addTargets(self.config['target'])
            spats.process_pair_data(self.r1, self.r2)
            spats.store(run_name)
        self._add_note("wrote output to {}".format(os.path.basename(run_name)))
        self._notebook().add_spats_run(self.cotrans, spats.run.count_mutations).save()

    def _update_run_config(self, run):
        custom_config = False
        sentinel = '_-=*< sEnTiNeL >*-=_'
        for key, value in self.config.iteritems():
            if key in [ "r1", "r2", "preseq", "target", "cotrans" ]:
                continue
            if sentinel != getattr(run, key, sentinel):
                try:
                    val = ast.literal_eval(value)
                except:
                    val = value
                setattr(run, key, val)
                self._add_note("config set {} = {}".format(key, val))
                custom_config = True
            else:
                self._add_note("warning: unknown config {}".format(key))
        return custom_config

    def _spats_file(self, base_name):
        return os.path.join(self.path, '{}.spats'.format(base_name))

    def _pre_file(self, key):
        return self._spats_file('pre_{}'.format(key))

    def _run_file(self):
        return self._spats_file('run')

    def _reads_file(self):
        return self._spats_file('reads')

    def _notebook(self):
        nb = nbutil.Notebook(os.path.join(self.path, 'spats.ipynb'))
        if nb.is_empty():
            nb.add_metadata(self.metadata)
            nb.add_initializer()
        return nb

    def validate(self):
        """Validate the results of a previous 'process' run against a second (slower) algorithm.
        """

        run_name = self._run_file()
        if not os.path.exists(run_name):
            raise Exception("Run must be performed before validating")

        spats = Spats()
        spats.load(run_name)
        if spats.validate_results(self.r1, self.r2):
            self._add_note("Validation pass")
        else:
            self._add_note("Validation FAILURE")

    def _install_nbextensions(self):
        ext_out = subprocess.check_output(["jupyter", "nbextension", "list", "--user"])
        if "spats_shape_seq/main" not in ext_out:
            subprocess.check_call(["jupyter", "nbextension", "install", "--user", "--py", "spats_shape_seq"])
            subprocess.check_call(["jupyter", "nbextension", "enable",  "--user", "--py", "spats_shape_seq"])

    def _install_jupyter_browser_fix(self):
        jup_conf_path = os.path.expanduser('~/.jupyter/jupyter_notebook_config.py')
        jup_conf_line = "c.NotebookApp.browser = u'open %s'\n"
        if os.path.exists(jup_conf_path):
            jup_conf = open(jup_conf_path, 'rb').read()
            if 'c.NotebookApp.browser' in jup_conf:
                return
            jup_conf += "\n" + jup_conf_line
        else:
            jup_conf = jup_conf_line
        open(jup_conf_path, 'wb').write(jup_conf)

    def _install_matplotlib_styles(self):
        import matplotlib as mpl
        conf_dir = mpl.get_configdir()
        if not os.path.exists(conf_dir):
            os.mkdir(conf_dir)
        style_dir = os.path.join(conf_dir, 'stylelib')
        if not os.path.exists(style_dir):
            os.mkdir(style_dir)
        static_path = os.path.join(self._module_path(), 'static', 'styles')
        for style in os.listdir(static_path):
            target_path = os.path.join(style_dir, style)
            if not os.path.exists(target_path):
                shutil.copyfile(os.path.join(static_path, style), target_path)

    def nb(self):
        """Launch the Jupyter notebook.
        """
        self._skip_log = True

        self._install_nbextensions()
        self._install_matplotlib_styles()
        self._install_jupyter_browser_fix()

        try:
            process = subprocess.Popen(["jupyter", "notebook", "-y", "spats.ipynb"], cwd = self.path)
            process.wait()
        except KeyboardInterrupt:
            process.terminate()
            time.sleep(0.4)
            process.wait()

    def viz(self):
        """Launch the visualization tool UI.
        """

        self._skip_log = True
        if sys.platform != "darwin":
            raise Exception("Invalid platform for viz UI: {}".format(sys.platform))

        subprocess.check_call(["make", "vizprep"], cwd = self._spats_path())
        def viz_handler():
            from viz.ui import SpatsViz
            sv = SpatsViz()
            sv.stop_on_disconnect = True
            sv.start()
            sv.waitFor()
        def uiclient_handler():
            bin_path = os.path.join(self._spats_path(), "bin", "UIClient.app", "Contents", "MacOS", "UIClient")
            try:
                subprocess.call([bin_path], cwd = self._spats_path())
            except KeyboardInterrupt:
                pass
        viz_worker = multiprocessing.Process(target = viz_handler, args = [])
        viz_worker.start()
        uiclient_worker = multiprocessing.Process(target = uiclient_handler, args = [])
        time.sleep(0.1)
        uiclient_worker.start()
        try:
            while viz_worker.is_alive() and uiclient_worker.is_alive():
                viz_worker.join(0.1)
                uiclient_worker.join(0.1)
        except KeyboardInterrupt:
            pass
        if viz_worker.is_alive():
            viz_worker.join()
        if uiclient_worker.is_alive():
            uiclient_worker.terminate()

    def dump(self):
        """Dump data. Provide 'reads' or 'run' as an argument to dump the
        indicated type of data.
        """
        self._skip_log = True
        if not self._command_args:
            raise Exception("Dump requires a type (either 'reads' or 'run').")
        dump_type = self._command_args[0]
        handler = getattr(self, "_dump_" + dump_type, None)
        if not handler:
            raise Exception("Invalid dump type: {}".format(dump_type))
        handler()

    def _dump_prefixes(self):
        run_name = self._run_file()
        if not os.path.exists(run_name):
            raise Exception("Run must be run before attempting dump")

        spats = Spats()
        spats.load(run_name)
        countinfo = spats.counters.counts_dict()
        total = float(countinfo['total_pairs']) / 100.0
        for mask in spats.run.masks:
            prefixes = []
            for key in sorted([k for k in countinfo.keys() if k.startswith('prefix_{}_'.format(mask))], key = lambda k : countinfo[k], reverse = True):
                prefixes.append((key[12:], float(countinfo[key]) / total, countinfo[key]))
            output_path = os.path.join(self.path, 'prefixes_{}.csv'.format(mask))
            self._write_csv(output_path, [ "Tag", "Percentage", "Count" ], prefixes)

    def _dump_reads(self):
        reads_name = self._reads_file()
        if not os.path.exists(reads_name):
            raise Exception("Reads must be run before attempting dump")
        data = ReadsData(reads_name)
        db = data.pair_db
        counts = db.tag_counts(1)
        total = float(db.count()) / 100.0
        keys = sorted(counts.keys(), key = lambda x : counts[x], reverse = True)
        data = [ [ key, float(counts[key])/total, counts[key] ] for key in keys ]
        output_path = os.path.join(self.path, 'reads.csv')
        self._write_csv(output_path, [ "Tag", "Percentage", "Count" ], data)

    def _dump_run(self):
        run_name = self._run_file()
        if not os.path.exists(run_name):
            raise Exception("Run must be run before attempting dump")

        spats = Spats()
        spats.load(run_name)
        profiles = spats.compute_profiles()
        mutations = spats.run.count_mutations
        if mutations:
            headers = [ "L", "site", "nt", "f+", "f-", "mut+", "mut-", "beta", "mu", "r", "c" ]
        else:
            headers = [ "L", "site", "nt", "f+", "f-", "beta", "theta", "rho", "c" ]
        data = []

        if self.cotrans:
            tgt = spats.targets.targets[0]
            tseq = tgt.seq
            for key in profiles.cotrans_keys():
                end = int(key.split('_')[-1])
                prof = profiles.profilesForTargetAndEnd(tgt.name, end)
                for i in range(end + 1):
                    if mutations:
                        data.append([ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i], prof.treated_muts[i], prof.untreated_muts[i], prof.beta[i], prof.mu[i], prof.r_mut[i], prof.c ])
                    else:
                        data.append([ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i], prof.beta[i], prof.theta[i], prof.rho[i], prof.c ])
            output_path = os.path.join(self.path, '{}.csv'.format(tgt.name))
            self._write_csv(output_path, headers, data)
            empty_cell = ''
            if mutations:
                keys = [ 'treated', 'untreated', 'treated_mut', 'untreated_mut', 'beta', 'mu', 'r' ]
            else:
                keys = [ 'treated', 'untreated', 'beta', 'theta', 'rho' ]
            for key in keys:
                cotrans_keys = profiles.cotrans_keys()
                ncols = 0
                mat = []
                for pkey in cotrans_keys:
                    end = int(pkey.split('_')[-1])
                    prof = profiles.profilesForTargetAndEnd(tgt.name, end)
                    vals = getattr(prof, key)
                    if not ncols:
                        ncols = len(cotrans_keys) + len(vals)
                    if len(vals) < ncols:
                        vals += ([empty_cell] * (ncols - len(vals)))
                    mat.append(vals)
                self._write_csv('{}_{}_mat.csv'.format(tgt.name, key), None, mat)
        else:
            for tgt in spats.targets.targets:
                tseq = tgt.seq
                end = len(tgt.seq)
                prof = profiles.profilesForTarget(tgt)
                data = []
                for i in range(end + 1):
                    if mutations:
                        data.append([ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i], prof.treated_muts[i], prof.untreated_muts[i], prof.beta[i], prof.mu[i], prof.r_mut[i], prof.c ])
                    else:
                        data.append([ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i], prof.beta[i], prof.theta[i], prof.rho[i], prof.c ])
                output_path = os.path.join(self.path, '{}.csv'.format(tgt.name))
                self._write_csv(output_path, headers, data)


    def _write_csv(self, output_path, headers, data):
        with open(output_path, 'wb') as out_file:
            writer = csv.writer(out_file)
            if headers:
                writer.writerow(headers)
            for row in data:
                writer.writerow(row)
        self._add_note("Data dumped to {}".format(os.path.basename(output_path)))

    def handle_filter(self):
        """Generates an output of the demultiplexed positive and negative
           fastq files for R1 and R2.
        """
        self._skip_log = True
        files = fastq_handle_filter(self.r1, self.r2)
        self._add_note("Pairs filtered to: {}".format(", ".join(files)))

    def doc(self):
        """Show the spats documentation.
        """

        self._skip_log = True
        if spats_shape_seq._PUBLIC_RELEASE:
            import webbrowser
            webbrowser.open('http://spats.readthedocs.io/')
        else:
            subprocess.check_call(["make", "showdocs"], cwd = self._spats_path())

    def help(self):
        """Display available commands.
        """

        self._skip_log = True
        print("\nspats_tool v{}\nCommands:\n".format(spats_shape_seq._VERSION))
        for key in sorted(SpatsTool.__dict__.keys()):
            if key.startswith('_'):
                continue
            if key in self._private_commands and spats_shape_seq._PUBLIC_RELEASE:
                continue
            value = getattr(SpatsTool, key)
            if value.__doc__:
                print("  {}:  {}".format(key, value.__doc__.rstrip()))
        print("\nData commands require 'spats.config' to be present in the current working directory,")
        print("which requires a [spats] section header and should at least specify the 'target', 'r1',")
        print("and 'r2' configuration keys.\n")


_spats_config_template = """
[spats]

# set to True/False depending on whether this is a cotrans experiment
cotrans = True

# Required for pre-sequencing tool: the path to the ABIF (.fsa) file. Can also be a comma-separated list for multiple files.
#preseq = 

# Required for SPATS runs and reads analysis: the path to the targets FASTA file.
#target = 

# Required for SPATS runs and reads analysis: the paths to the R1/R2 data files.
#r1 = 
#r2 = 


# Known metadata. Recommended to provide all applicable fields.
[metadata]

# Experiment name
#name = 

# Author / Experimenter (initials)
#author = 

# Date: this will be filled in automatically (today's date/time) unless explicitly provided.
#date = 

# Folding conditions:
#buffer = 
#temperature = 
#salt = 

#probe = 

# Lot numbers:
#adapter = 
#enzyme = 
#reagent = 

"""

def run(args, path = None):
    SpatsTool(path)._run(args)


if __name__ == '__main__':
    run(sys.argv[1:])
