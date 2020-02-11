
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
from shutil import copyfile
import re
import ntpath

import spats_shape_seq
from spats_shape_seq import Spats
from spats_shape_seq.parse import abif_parse, fastq_handle_filter, FastFastqParser
from spats_shape_seq.reads import ReadsData, ReadsAnalyzer
from counters import Counters
from util import objdict_to_dict
from mask import PLUS_PLACEHOLDER, MINUS_PLACEHOLDER


class SpatsTool(object):

    def __init__(self, path):
        self.path = path or os.getcwd()
        self.config = None
        self.cotrans = False
        self._skip_log = False
        self._no_config_required_commands = [ "doc", "help", "init", "viz", "show", "extract_case", "add_case", "show_test_case" ]
        self._private_commands = [ "viz", "to_shapeware", "rerun" ]
        self._temp_files = []
        self._r1 = None
        self._r2 = None
        self._r1_plus = None
        self._r2_plus = None
        self._r1_minus = None
        self._r2_minus = None
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

    def _load_r1_r2(self, suffix = ''):
        def decomp(rx):
            base, ext = os.path.splitext(os.path.basename(rx))
            if ext.lower() == '.gz':
                self._add_note("decompressing {}".format(rx))
                out = os.path.join(self.path, base + ".tmp")
                subprocess.check_call('gzip -d -c "{}" > "{}"'.format(rx, out), cwd = self.path, shell = True)
                self._temp_files.append(out)
                self._sentinel("decompress {}".format(rx))
                return out
            return rx
        def singleOrList(rkey):
            # hack to keep r1_plus and r2_plus properties backwards compatible...
            res = [ decomp(r.strip()) for r in self.config[rkey].split(',') ]
            return res if len(res) > 1 else res[0]
        return singleOrList('r1' + suffix), singleOrList('r2' + suffix)

    @property
    def r1(self):
        if not self._r1:
            self._r1, self._r2 = self._load_r1_r2()
        return self._r1

    @property
    def r2(self):
        if not self._r2:
            self._r1, self._r2 = self._load_r1_r2()
        return self._r2

    @property
    def using_separate_channel_files(self):
        return bool(self.config.get('r1_plus'))

    @property
    def r1_plus(self):
        if not self._r1_plus:
            self._r1_plus, self._r2_plus = self._load_r1_r2('_plus')
        return self._r1_plus

    @property
    def r2_plus(self):
        if not self._r2_plus:
            self._r1_plus, self._r2_plus = self._load_r1_r2('_plus')
        return self._r2_plus

    @property
    def plus_channels(self):
        r1p, r2p = self.r1_plus, self.r2_plus
        if not isinstance(r1p, list):
            r1p = [ r1p ]
        if not isinstance(r2p, list):
            r2p = [ r2p ]
        if len(r1p) != len(r2p):
            raise Exception("r1/r2 plus channels do not correspond")
        return zip(r1p, r2p)

    @property
    def r1_minus(self):
        if not self._r1_minus:
            self._r1_minus, self._r2_minus = self._load_r1_r2('_minus')
        return self._r1_minus

    @property
    def r2_minus(self):
        if not self._r2_minus:
            self._r1_minus, self._r2_minus = self._load_r1_r2('_minus')
        return self._r2_minus


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
            #raise

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
            if self.using_separate_channel_files:
                pcs = self.plus_channels
                if len(pcs) > 1:
                    raise Exception("multiple positive r1 channel files not supported for reads tool.")
                data.parse(self.config['target'], [pcs[0][0], self.r1_minus], [pcs[0][1], self.r2_minus])
            else:
                data.parse(self.config['target'], [self.r1], [self.r2])

        analyzer = ReadsAnalyzer(data, cotrans = self.cotrans)
        self._update_run_config(analyzer.run)

        # xref https://trello.com/c/VMFyZjjg/286-handle-quality-parsing-in-reads-tool-if-configured
        # to do this, we'd need to parse quality to the db (nontrivial)
        # for now just force-disable this, as it's not required to do reads analysis
        analyzer.run.mutations_require_quality_score = None

        analyzer.process_tags()
        self._add_note("tags processed to {}".format(os.path.basename(db_name)))

    def rerun(self):
        """ Re-does reads analysis for pairs that have a tag into a new result set to allow comparing results across config options.
        """
        if not self._command_args  or  len(self._command_args) != 4:
            raise Exception("usage: spats_tool rerun tag spats_run result_set_name cmp_set_name\n\t- tag indicating the pairs to rerun\n\t- path to a spats database file\n\t- name for the new result set\n\t- name of result set to compare against.")
        tag = self._command_args[0]
        spatsdb = self._command_args[1]
        result_set_name = self._command_args[2]
        cmp_set_name = self._command_args[3]
        if not os.path.exists(spatsdb):
            raise Exception("spats_run file does not exist at path {}".format(spatsdb))

        native_tool = self._native_tool('reads')
        if native_tool:
            self._add_note("using native reads")
            subprocess.check_call([native_tool, self.config['target'], self.r1, self.r2, db_name], cwd = self.path)

        data = ReadsData(spatsdb)
        if data.pair_db.result_set_id_for_name(result_set_name):
            raise Exception("result_set_name '{}' already exists in spats db".format(result_set_name))
        cmp_set_id = data.pair_db.result_set_id_for_name(cmp_set_name)
        if not cmp_set_id:
            raise Exception("cmp_set_name '{}' does not exist in spats db".format(cmp_set_name))
        if not native_tool:
            self._add_note("using python reads")
            if self.using_separate_channel_files:
                raise Exception('rerun tool not supported with separate channel files')

        analyzer = ReadsAnalyzer(data, cotrans = self.cotrans)
        data.pair_db.load_run(analyzer.run)
        self._update_run_config(analyzer.run)

        analyzer.run._redo_tag = tag
        analyzer.run.result_set_name = result_set_name
        analyzer.run.cmp_set_id = cmp_set_id

        # xref https://trello.com/c/VMFyZjjg/286-handle-quality-parsing-in-reads-tool-if-configured
        # to do this, we'd need to parse quality to the db (nontrivial)
        # for now just force-disable this, as it's not required to do reads analysis
        analyzer.run.mutations_require_quality_score = None

        redo_tags = data.pair_db.tag_counts(cmp_set_id, [ tag ])
        self._add_note("rerunning {} pairs (incl. multiples) with tag '{}'...".format(redo_tags[tag], tag))
        analyzer.process_tags()
        self._add_note("tags processed to {}".format(os.path.basename(spatsdb)))

        num_done = data.pair_db.num_results(result_set_name)
        diff_failures = 0
        new_failures = 0
        fixed_failures = 0
        diff_sites = 0
        diff_ends = 0
        for res in data.pair_db.differing_results(result_set_name, cmp_set_name):
            assert(res[13] == res[14])
            multiplicity = min(int(res[13]), int(res[14]))
            if res[7] != res[12]:
                diff_failures += multiplicity
                if res[7] and not res[12]:
                    new_failures += multiplicity
                elif res[12] and not res[7]:
                    fixed_failures += multiplicity
            else:
                if res[5] != res[10]:
                    diff_sites += multiplicity
                if res[4] != res[9]:
                    diff_ends += multiplicity
        print("New result set '{}' for tagged pairs added to {}.".format(result_set_name, spatsdb))
        print("Quick comparison of {} unique tagged pairs with result_set '{}' yielded:  ".format(num_done, cmp_set_name))
        print("\t- {} pairs with different failures".format(diff_failures))
        print("\t    - {} failures were fixed".format(fixed_failures))
        print("\t    - {} failures were new".format(new_failures))
        print("\t- {} pairs with different sites".format(diff_sites))
        print("\t- {} pairs with different ends".format(diff_ends))


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
            nb = self._notebook()
            if nb:
                nb.add_preseq(key).save()

    def _run_plots(self):
        nb = self._notebook()
        if nb:
            nb.add_spats_run(self.cotrans, True).save()

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
            if self.using_separate_channel_files:
                self._add_note("using separate channel files.")
                self._add_note("processing minus/untreated pairs")
                spats.process_pair_data(self.r1_minus, self.r2_minus, force_mask = (spats.run.masks[1] if spats.run.masks[1] else MINUS_PLACEHOLDER))
                pcs = self.plus_channels
                if len(pcs) == 1:
                    self._add_note("processing plus/treated pairs")
                    spats.process_pair_data(pcs[0][0], pcs[0][1], force_mask = (spats.run.masks[0] if spats.run.masks[0] else PLUS_PLACEHOLDER))
                    spats.store(run_name)
                else:
                    if len(spats.run.masks) >= 2 and (not spats.run.masks[0] or not spats.run.masks[1]):
                        raise Exception("empty masks with split channels not supported for reads tool.")
                    spats.store(run_name)
                    for i, (r1_plus, r2_plus) in enumerate(pcs):
                        self._add_note("processing plus/treated pairs, set #{}".format(i + 1))
                        spats._processor.counters.reset()
                        spats.process_pair_data(r1_plus, r2_plus, force_mask = (spats.run.masks[0] if spats.run.masks[0] else PLUS_PLACEHOLDER))
                        spats.store(run_name + ".p{}".format(i + 1))
            else:
                spats.process_pair_data(self.r1, self.r2)
                spats.store(run_name)
        self._add_note("wrote output to {}".format(os.path.basename(run_name)))
        nb = self._notebook()
        if nb:
            nb.add_spats_run(self.cotrans, spats.run.count_mutations).save()

    def _update_run_config(self, run, dictionary = None):
        custom_config = False
        sentinel = '_-=*< sEnTiNeL >*-=_'
        for key, value in self.config.iteritems():
            if key in [ "r1", "r2", "r1_plus", "r2_plus", "r1_minus", "r2_minus", "preseq", "target", "cotrans" ]:
                continue
            if sentinel != getattr(run, key, sentinel):
                try:
                    val = ast.literal_eval(value)
                except:
                    val = value
                #print("overwriting with", key, val)
                setattr(run, key, val)
                self._add_note("config set {} = {}".format(key, val))
                if dictionary:
                    dictionary[key] = val
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
        try:
            import spats_shape_seq.nbutil as nbutil
        except:
            return None
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

        if not self._notebook():
            raise Exception('Notebook requires the jupyter and nbformat pacakges. Try "pip install nbformat jupyter".')

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
        """Dump data. Provide 'reads', 'run', 'prefixes', 'mut_counts' or 'indel_lens'
        as an argument to dump the indicated type of data.
        """
        self._skip_log = True
        if not self._command_args:
            raise Exception("Dump requires a type ('reads', 'run', 'prefixes', 'mut_counts' or 'indel_lens').")
        dump_type = self._command_args[0]
        handler = getattr(self, "_dump_" + dump_type, None)
        if not handler:
            raise Exception("Invalid dump type: {}".format(dump_type))
        if dump_type == 'reads':
            self._dump_reads()
        else:
            self._dump_wrapper(handler)

    def _dump_reads(self):
        reads_name = self._reads_file()
        if not os.path.exists(reads_name):
            raise Exception("Reads must be run before attempting dump")
        data = ReadsData(reads_name)
        db = data.pair_db
        counts = db.tag_counts(1)
        total = float(db.count()) / 100.0    # because of this, cannot dump "rerun" data
        keys = sorted(counts.keys(), key = lambda x : counts[x], reverse = True)
        data = [ [ key, float(counts[key])/total, counts[key] ] for key in keys ]
        output_path = os.path.join(self.path, 'reads.csv')
        self._write_csv(output_path, [ "Tag", "Percentage", "Count" ], data)

    def _dump_wrapper(self, handler):
        run_name = self._run_file()
        if not os.path.exists(run_name):
            raise Exception("Run must be run before attempting dump")
        partial = False
        base_file = ntpath.basename(run_name)
        for filep in os.listdir(self.path):
            rn = re.match("^{}\.p(\d+)$".format(base_file), filep)
            if rn:
                spats = Spats()
                spats.load(run_name)
                spats.merge(filep)
                handler(spats, "p{}_".format(rn.group(1)))
                partial = True
        if not partial:
            spats = Spats()
            spats.load(run_name)
            handler(spats)

    def _dump_prefixes(self, spats, fprefix = ""):
        countinfo = spats.counters.counts_dict()
        total = float(countinfo['total_pairs']) / 100.0
        for mask in spats.run.masks:
            prefixes = []
            keyprefix = "prefix_{}_".format(mask)
            for key in sorted([k for k in countinfo.keys() if k.startswith(keyprefix)], key = lambda k : countinfo[k], reverse = True):
                prefixes.append((key[len(keyprefix):], float(countinfo[key]) / total, countinfo[key]))
            output_path = os.path.join(self.path, '{}prefixes_{}.csv'.format(fprefix, mask))
            self._write_csv(output_path, [ "Tag", "Percentage", "Count" ], prefixes)
        total = float(countinfo['registered_pairs']) / 100.0
        for mask in spats.run.masks:
            prefixes = []
            keyprefix = "mapped_prefix_{}_".format(mask)
            for key in sorted([k for k in countinfo.keys() if k.startswith(keyprefix)], key = lambda k : countinfo[k], reverse = True):
                prefixes.append((key[len(keyprefix):], float(countinfo[key]) / total, countinfo[key]))
            output_path = os.path.join(self.path, '{}mapped_prefixes_{}.csv'.format(fprefix, mask))
            self._write_csv(output_path, [ "Tag", "Percentage", "Count" ], prefixes)

    def _dump_mut_counts(self, spats, prefix = ""):
        countinfo = spats.counters.counts_dict()
        mut_cnts = []
        for muts in sorted([int(k.split('_')[-1]) for k in countinfo.keys() if k.startswith('mut_count_')]):
            mut_cnts.append((muts, countinfo["mut_count_{}".format(muts)]))
        output_path = os.path.join(self.path, '{}mut_counts.csv'.format(prefix))
        self._write_csv(output_path, [ "Mutation Count", "Reads" ], mut_cnts)
        mut_cnts = []
        for muts in sorted([int(k.split('_')[-1]) for k in countinfo.keys() if k.startswith('mapped_mut_count_')]):
            mut_cnts.append((muts, countinfo["mapped_mut_count_{}".format(muts)]))
        output_path = os.path.join(self.path, '{}mapped_mut_counts.csv'.format(prefix))
        self._write_csv(output_path, [ "Mutation Count", "Reads" ], mut_cnts)

    def _dump_indel_lens(self, spats, prefix = ""):
        countinfo = spats.counters.counts_dict()
        ilen_cnt = []
        for lc in sorted([int(k.split('_')[-1]) for k in countinfo.keys() if k.startswith('mapped_indel_len_')]):
            ilen_cnt.append((lc, countinfo["mapped_indel_len_{}".format(lc)]))
        output_path = os.path.join(self.path, '{}mapped_indel_len_counts.csv'.format(prefix))
        self._write_csv(output_path, [ "Indel Length", "Reads" ], ilen_cnt)

    def _dump_run(self, spats, prefix = ""):
        profiles = spats.compute_profiles()
        mutations = spats.run.count_mutations
        indels = spats.run.handle_indels
        headers = [ "L", "site", "nt", "f+", "f-" ]
        if indels:
            headers += [ "ins+", "ins-", "del+", "del-" ]
        if mutations:
            headers += [ "mut+", "mut-", "beta", "mu", "r" ]
        else:
            headers += [ "beta", "theta", "rho" ]
        headers += [ "c thresh", "c" ]
        data = []
        sites_missing_reads = []

        if self.cotrans:
            tgt = spats.targets.targets[0]
            tseq = tgt.seq
            for key in profiles.cotrans_keys():
                end = int(key.split('_')[-1])
                prof = profiles.profilesForTargetAndEnd(tgt.name, end)
                for i in xrange(end + 1):
                    if 0 == prof.treated[i] and 0 == prof.untreated[i]:
                        sites_missing_reads.append( (tgt, end, i) )
                    datapt = [ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i] ]
                    if indels:
                        datapt += [ prof.treated_inserts[i], prof.untreated_inserts[i], prof.treated_deletes[i], prof.untreated_deletes[i] ]
                    if mutations:
                        datapt += [ prof.treated_muts[i], prof.untreated_muts[i], prof.beta[i], prof.mu[i], prof.r_mut[i] ]
                    else:
                        datapt += [ prof.beta[i], prof.theta[i], prof.rho[i] ]
                    datapt += [ prof.c_thresh, prof.c ]
                    data.append(datapt)
            output_path = os.path.join(self.path, '{}{}.csv'.format(prefix, tgt.name))
            self._write_csv(output_path, headers, data)
            empty_cell = ''
            keys = [ 'treated', 'untreated' ]
            if indels:
                keys += [ 'treated_inserts', 'untreated_inserts', 'treated_deletes', 'untreated_deletes' ]
            if mutations:
                keys += [ 'treated_mut', 'untreated_mut', 'beta', 'mu', 'r' ]
            else:
                keys += [ 'beta', 'theta', 'rho' ]
            cotrans_keys = profiles.cotrans_keys()
            for key in keys:
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
                self._write_csv('{}{}_{}_mat.csv'.format(prefix, tgt.name, key), None, mat)
        else:
            for tgt in spats.targets.targets:
                tseq = tgt.seq
                end = len(tgt.seq)
                prof = profiles.profilesForTarget(tgt)
                data = []
                for i in xrange(end + 1):
                    if 0 == prof.treated[i] and 0 == prof.untreated[i]:
                        sites_missing_reads.append( (tgt, end, i) )
                    datapt = [ end, i, tseq[i - 1] if i else '*', prof.treated[i], prof.untreated[i] ]
                    if indels:
                        datapt += [ prof.treated_inserts[i], prof.untreated_inserts[i], prof.treated_deletes[i], prof.untreated_deletes[i] ]
                    if mutations:
                        datapt += [ prof.treated_muts[i], prof.untreated_muts[i], prof.beta[i], prof.mu[i], prof.r_mut[i] ]
                    else:
                        datapt += [ prof.beta[i], prof.theta[i], prof.rho[i] ]
                    datapt += [ prof.c_thresh, prof.c ]
                    data.append(datapt)
                output_path = os.path.join(self.path, '{}{}.csv'.format(prefix, tgt.name))
                self._write_csv(output_path, headers, data)

        missing_targets = list(set([ s[0] for s in sites_missing_reads ]))
        for tgt in missing_targets:
            # only warn for sites "far from" the end (for now defined by 2 * minimum required match length)
            min_len = spats.run.minimum_target_match_length
            tgt_missing_sites = [ s for s in sites_missing_reads if s[0] == tgt and s[1] - s[2] > 2 * min_len ]
            num_missing = len(tgt_missing_sites)
            if 0 == num_missing:
                # we might have pruned them due to length
                continue
            if self.cotrans:
                tgt_missing_sites = [ "{}/{}".format(s[1], s[2]) for s in tgt_missing_sites ]
            else:
                tgt_missing_sites = [ str(s[2]) for s in tgt_missing_sites ]
            if len(tgt_missing_sites) > 20:
                tgt_missing_sites = tgt_missing_sites[:20]
                tgt_missing_sites.append("...")
            print(" ** Warning: target {} has 0 reads from both channels at {} sites: {} ".format(tgt.name, num_missing, ", ".join(tgt_missing_sites)))


    def _write_csv(self, output_path, headers, data):
        with open(output_path, 'wb') as out_file:
            writer = csv.writer(out_file)
            if headers:
                writer.writerow(headers)
            for row in data:
                writer.writerow(row)
        self._add_note("Data dumped to {}".format(os.path.basename(output_path)))

    def plot(self):
        """Plot data, provide an argument for the plot type. For
        single-length, provide 'counts', 'rho', 'beta', 'theta',
        'muts', or 'muts_reactivity'. For cotrans, provide 'counts',
        'c', 'treated', 'untreated', 'treated_muts', 'untreated_muts',
        'rho', 'beta', 'theta', 'r_mut', 'mu_mut', 'beta_mut',
        'treated_mut', 'untreated_mut'.
        """
        self._skip_log = True
        if not self._command_args:
            raise Exception("Plot requires a type.")
        plot_type = self._command_args[0]
        import plots
        handler = getattr(plots, "plot_" + ("cotrans_" if self.cotrans else "sl_") + plot_type, None)
        if not handler:
            raise Exception("Invalid plot type: {}".format(plot_type))
        plt = handler()
        plt.show()

    def handle_filter(self):
        """Generates an output of the demultiplexed positive and negative
           fastq files for R1 and R2.
        """
        self._skip_log = True
        counters = Counters()
        files = fastq_handle_filter(self.r1, self.r2, counters=counters)
        self._add_note("{} pairs filtered to:\n   {}".format(counters.RRRY + counters.YYYR, "\n   ".join(files)))
        self._add_note("{} RRRY pairs included.".format(counters.RRRY))
        self._add_note("{} YYYR pairs included.".format(counters.YYYR))
        self._add_note("{} pairs without matching handle were not included.".format(counters.no_mask))

    def _pp_channel_files(self, channel_files):
        treated = "RRRY"
        untreated = "YYYR"
        for cfp in channel_files:
            cf = os.path.basename(cfp)
            subprocess.check_call('gzip "{}"'.format(cfp), cwd = os.path.dirname(cfp), shell = True)
            if cf.startswith('RRRY'):
                treated = cf.split('_')[0]
            elif cf.startswith('YYYR'):
                untreated = cf.split('_')[0]
        return treated, untreated

    def to_shapeware(self):
        """Create a folder from a spats dataset suitable for running
           with the SHAPEware tool produced by Arrakis.
        """
        self._skip_log = True
        if not self._command_args:
            raise Exception("to_shapeware requires a path be specified for the ouptut folder.")
        output_folder = self._command_args[0]
        if os.path.exists(output_folder):
            raise Exception("to_shapeware output folder already exists at {}".format(output_folder))
        data_dir = os.path.join(output_folder, "raw_data")
        os.makedirs(data_dir)
        counters = Counters()
        channel_files = fastq_handle_filter(self.r1, self.r2, strip_mask=True, outpath=data_dir, counters=counters)
        treated, untreated = self._pp_channel_files(channel_files)
        ## Note: the denatured sample is used as a "baseline" in SHAPEware, 
        ## for example:  reactivity = (mutr_treated - mutr_untreated) / mutr_denatured
        ## Since we don't have denatured, we'll use the untreated sample as the baseline.
        denatured = untreated
        experiment_name = self.r1.split('_')[0]
        target_file = self.config['target']
        with open(target_file, 'r') as TF:
            target_name = TF.readline().strip()
            if target_name[0] == ">":
                target_name = target_name[1:]
        copyfile(target_file, os.path.join(output_folder, "ref_seqs.fa"))
        copyfile(target_file, os.path.join(output_folder, "ref_masks.fa"))
        input_sheet = os.path.join(output_folder, "input_sheet.csv")
        with open(input_sheet, 'w') as IS:
            IS.write("Experiment name,Replicate,Reference sequence,Ligand,Sample with SHAPE reagent,Sample without SHAPE reagent,Denatured sample\n")
            IS.write("{},1,{},None,{},{},{}\n".format(experiment_name, target_name, treated, untreated, denatured))
        self._add_note("Created SHAPEware folder at {}.".format(output_folder))
        self._add_note("   {} treated pairs included.".format(counters.RRRY))
        self._add_note("   {} untreated pairs included.".format(counters.YYYR))
        self._add_note("   {} pairs without matching handle were not included.".format(counters.no_mask))

    def compare(self):

        from spats_shape_seq import Spats
        from spats_shape_seq.pair import Pair

        json_base = { 'target' : self.config['target'], 'config' : { 'algorithm' : 'find_partial', 'debug' : True }, 'expect' : {}}

        spats_fp = Spats(cotrans = self.cotrans)
        spats_lookup = Spats(cotrans = self.cotrans)
        self._update_run_config(spats_fp.run)
        self._update_run_config(spats_lookup.run, json_base['config'])
        spats_fp.run.algorithm = 'find_partial'
        spats_lookup.run.algorithm = 'lookup'

        spats_fp.addTargets(self.config['target'])
        spats_lookup.addTargets(self.config['target'])

        count = 0
        match = 0
        with FastFastqParser(self.r1, self.r2) as parser:
            total = parser.appx_number_of_pairs()
            for batch in parser.iterator(5000):
                for item in batch:
                    pair_fp = Pair()
                    pair_lookup = Pair()
                    pair_fp.set_from_data(str(item[0]), item[1], item[2])
                    pair_lookup.set_from_data(str(item[0]), item[1], item[2])
                    try:
                        spats_fp.process_pair(pair_fp)
                        spats_lookup.process_pair(pair_lookup)
                    except:
                        print('Error after {}/{}'.format(match, count))
                        raise
                    if (pair_fp.has_site == pair_lookup.has_site):
                        if not pair_fp.has_site:
                            count += 1
                            continue
                        elif (pair_fp.target.name == pair_lookup.target.name and
                              pair_fp.end == pair_lookup.end and
                              pair_fp.site == pair_lookup.site and
                              pair_fp.mutations == pair_lookup.mutations):
                            count += 1
                            match += 1
                            continue
                    json_base["id"] = str(item[0])
                    json_base["R1"] = str(item[1])
                    json_base["R2"] = str(item[2])
                    print('After {}/{} matches; mismatched pair: {} != {}\n{}'.format(match, count, pair_fp, pair_lookup,
                                                                                      json.dumps(json_base, sort_keys = True,indent = 4, separators = (',', ': '))))
                    return
                print('{}/{}-{}...'.format(match, count, total))
        spats_fp.counters.total_pairs = count
        spats_lookup.counters.total_pairs = count
        print('All match {}/{}.'.format(match, count))
        print(spats_fp._report_counts())
        print(spats_lookup._report_counts())

    def _test_case_registry(self):
        import spats_shape_seq.tests.test_harness
        return spats_shape_seq.tests.test_harness.registry()

    def extract_case(self):
        """Extracts a test case from the registry.
        """

        self._skip_log = True
        if not self._command_args or len(self._command_args) < 2:
            raise Exception("extract requires a test case id and an output filename")
        case_id = self._command_args[0]
        test_case_file = self._command_args[1]

        reg = self._test_case_registry()
        case = reg.extract_case(case_id)
        if not case:
            raise Exception('Unknown case id: {}'.format(case_id))

        open(test_case_file, 'w').write(json.dumps(case.jsonDict(), sort_keys = True, indent = 4, separators = (',', ': ')))
        print("Extracted case '{}' to '{}'".format(case_id, test_case_file))

    def add_case(self):
        """Adds a test case from the registry.
        """

        self._skip_log = True
        if not self._command_args:
            raise Exception("add requires a test case file")

        test_case_file = self._command_args[0]
        test_case = json.loads(open(test_case_file, 'r').read())

        reg = self._test_case_registry()
        reg.add_case(test_case)
        print("Added case '{}' to set '{}' in the test case registry.".format(test_case['id'], test_case['set_name']))

    def show_test_case(self):
        """Shows the diagram and result for analysis of a unit test case.
        """

        self._skip_log = True
        if not self._command_args:
            raise Exception("show_test_case requires a test id")

        test_case_id = self._command_args[0]
        print(test_case_id)
        reg = self._test_case_registry()
        test_case = reg.extract_case(test_case_id)
        test_case.run_opts['debug'] = True
        print(json.dumps(test_case.jsonDict(), sort_keys = True, indent = 4, separators = (',', ': ')))
        self._show_case(test_case)


    def show(self):
        """Shows the diagram and result for analysis of a test case pair.
        """

        import spats_shape_seq.tests.test_harness
        self._skip_log = True

        if not self._command_args:
            raise Exception("show requires the path to a test case")

        test_case_file = self._command_args[0]
        test_case_dict = json.loads(open(test_case_file, 'r').read())
        test_case = spats_shape_seq.tests.test_harness.SpatsCase(test_case_dict)
        self._show_case(test_case)


    def _show_case(self, test_case):

        from spats_shape_seq import Spats
        from spats_shape_seq.diagram import diagram

        alg = test_case.run_opts.get('algorithm')
        algs = [ alg ] if alg else test_case.run_opts.get('algorithms', [ 'find_partial', 'lookup' ])
        for algorithm in algs:
            spats = Spats()
            spats.run.algorithm = algorithm

            for key, value in test_case.run_opts.items():
                if str(key) == 'algorithms':
                    continue
                if isinstance(value, unicode):
                    value = str(value)
                if not hasattr(spats.run, key):
                    raise Exception('Invalid run_opt: {}'.format(key))
                setattr(spats.run, key, value)

            for name, seq in test_case.targets.iteritems():
                spats.addTarget(name, seq)

            pair = test_case.pair()
            if len(algs) > 1:
                print('\n[[ ALGORITHM: {} ]]'.format(algorithm))
            spats.process_pair(pair)
            if test_case.comment:
                print('Comment: {}'.format(test_case.comment))
            print(diagram(pair, spats.run))

            if test_case.expect:
                # should mirror `_check_expect` in test_harness.py...
                expects = test_case.expect
                fail = False

                try:

                    if expects['site'] is None:
                        if pair.site is not None:
                            raise Exception("pair.site={} when expecting none.".format(pair.site))
                    else:
                        if pair.site is None:
                            raise Exception("pair.site is none when expecting {}.".format(expects['site']))
                        if pair.site != expects['site']:
                            raise Exception("pair.site={} != expect.site={}".format(pair.site, expects['site']))
                        if 'end' in expects and pair.end != expects['end']:
                            raise Exception("pair.end={} != expect.end={}".format(pair.end, expects['end']))
                        if 'muts' in expects:
                            if expects['muts'] is not None  and  len(expects['muts']) > 0:
                                if not sorted(expects['muts']) == (sorted(pair.mutations) if pair.mutations else pair.mutations):
                                    raise Exception("mismatching mutations:  expected={}, pair.mutations={}".format(expects['muts'], pair.mutations))
                            else:
                                if not (pair.mutations is None  or len(pair.mutations) == 0):
                                    raise Exception("unexpected mutations: {}".format(pair.mutations))
                        if 'r1_indels' in expects:
                            r1inds = objdict_to_dict(pair.r1.indels)
                            if expects['r1_indels']:
                                if expects['r1_indels'] != r1inds:
                                    raise Exception("mismatching R1 indels:  expected={}, pair.r1.indels={}".format(expects['r1_indels'], r1inds))
                            elif pair.r1.indels:
                                raise Exception("unexpected R1 indels:  pair.r1.indels={}".format(pair.r1.indels))
                        if 'r2_indels' in expects:
                            r2inds = objdict_to_dict(pair.r2.indels)
                            if expects['r2_indels']:
                                if expects['r2_indels'] != r2inds:
                                    raise Exception("mismatching R2 indels:  expected={}, pair.r2.indels={}".format(expects['r2_indels'], r2inds))
                            elif pair.r2.indels:
                                raise Exception("unexpected R2 indels:  pair.r2.indels={}".format(r2inds))
                        if 'counters' in expects:
                            for counter, value in expects['counters'].iteritems():
                                if getattr(spats.counters, str(counter)) != value:
                                    raise Exception("counter '{}' value off: expected={} != got={}".format(counter, value, getattr(spats.counters, counter)))
                        if 'pair.target' in expects:
                            tname = pair.target.name if pair.target else None
                            if tname != expects['pair.target']:
                                raise Exception("pair.target={} != expect.pair.target={}".format(tname, expects['pair.target']))

                except Exception as e:
                    print('FAIL: {}'.format(e))
                    sys.exit(1)

                print('PASS')

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
