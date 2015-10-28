#!/usr/bin/env python

#~ Copyright 2015 Wieger Wesselink.
#~ Distributed under the Boost Software License, Version 1.0.
#~ (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

import os
import os.path
import re
import sys
sys.path += [os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python'))]

from random_bes_generator import make_bes
from random_pbes_generator import make_pbes
from random_process_generator import make_process_specification, generator_map, make_action, make_delta, make_tau, \
    make_process_instance, make_sum, make_if_then, make_if_then_else, make_choice, make_seq
from testing import run_pbes_test_with_counter_example_minimization
from testcommand import YmlTest
from text_utility import write_text

MCRL2_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
MCRL2_INSTALL_DIR = os.path.join(MCRL2_ROOT, 'install', 'bin')

def ymlfile(file):
    return '{}/tests/specifications/{}.yml'.format(MCRL2_ROOT, file)

def mcrl2file(file):
    return os.path.join(MCRL2_ROOT, file)

class RandomTest(YmlTest):
    def __init__(self, name, ymlfile, settings = dict()):
        super(RandomTest, self).__init__(name, ymlfile, [], settings)

    # create input files for the random test, and add the filenames to self.inputfiles
    def create_inputfiles(self, runpath = '.'):
        raise NotImplemented

    # removes input files that are in the runpath directory
    def remove_inputfiles(self, runpath = '.'):
        for filename in self.inputfiles:
            if os.path.abspath(runpath) == os.path.abspath(os.path.dirname(filename)):
                os.remove(filename)

    def execute(self, runpath = '.'):
        self.create_inputfiles(runpath)
        super(RandomTest, self).execute(runpath)
        self.remove_inputfiles(runpath)

class ProcessTest(RandomTest):
    def __init__(self, name, ymlfile, settings = dict()):
        super(ProcessTest, self).__init__(name, ymlfile, settings)
        self.generator_map = generator_map
        self.actions = ['a', 'b', 'c', 'd']
        self.process_identifiers = ['P', 'Q', 'R']
        self.process_size = 13
        self.init = None
        self.generate_process_parameters = False

    def create_inputfiles(self, runpath = '.'):
        filename = '{0}.mcrl2'.format(self.name, self.settings)
        p = make_process_specification(self.generator_map, self.actions, self.process_identifiers, self.process_size, init = self.init, generate_process_parameters = self.generate_process_parameters)
        write_text(filename, str(p))
        self.inputfiles += [filename]

# generates random process with higher probability of tau transitions
class ProcessTauTest(ProcessTest):
    def __init__(self, name, testfile, settings = dict()):
        super(ProcessTauTest, self).__init__(name, testfile, settings)
        self.actions = ['a', 'b', 'c']
        self.init = 'hide({a}, allow({a, b, c}, P || Q || R))'
        self.generator_map = {
                               make_action: 8,
                               make_delta: 1,
                               make_tau: 4,
                               make_process_instance: 1,
                               make_sum: 0,
                               make_if_then: 0,
                               make_if_then_else: 0,
                               make_choice: 5,
                               make_seq: 5,
                             }

class AlphabetTest(ProcessTest):
    def __init__(self, name, settings = dict()):
        super(AlphabetTest, self).__init__(name, ymlfile('alphabet'), settings)

class LpsSuminstTest(ProcessTest):
    def __init__(self, name, settings = dict()):
        super(LpsSuminstTest, self).__init__(name, ymlfile('lpssuminst'), settings)

class LpsSumelmTest(ProcessTest):
    def __init__(self, name, settings = dict()):
        super(LpsSumelmTest, self).__init__(name, ymlfile('lpssumelm'), settings)

class LpsParelmTest(ProcessTest):
    def __init__(self, name, settings = dict()):
        super(LpsParelmTest, self).__init__(name, ymlfile('lpsparelm'), settings)
        self.generate_process_parameters = True

class LpsConfcheckTest(ProcessTauTest):
    def __init__(self, name, confluence_type, settings = dict()):
        self.option_map = { 'commutative' : 'C',
                            'commutative-disjoint' : 'c',
                            'disjoint' : 'd',
                            'triangular' : 'T',
                            'trivial' : 'Z'
                          }
        assert confluence_type in self.option_map
        super(LpsConfcheckTest, self).__init__(name, ymlfile('lpsconfcheck'), settings)
        self.set_command_line_options('t2', ['-x' + self.option_map[confluence_type]])

class LtscompareTest(ProcessTauTest):
    def __init__(self, name, equivalence_type, settings = dict()):
        assert equivalence_type in ['bisim', 'branching-bisim', 'dpbranching-bisim', 'weak-bisim', 'dpweak-bisim', 'sim', 'trace', 'weak-trace']
        super(LtscompareTest, self).__init__(name, ymlfile('ltscompare'), settings)
        self.set_command_line_options('t3', ['-e' + equivalence_type])
        self.set_command_line_options('t4', ['-e' + equivalence_type])

class BisimulationTest(ProcessTauTest):
    def __init__(self, name, equivalence_type, settings = dict()):
        assert equivalence_type in ['bisim', 'branching-bisim', 'weak-bisim']
        bisimulation_type = equivalence_type
        if bisimulation_type == 'bisim':
            bisimulation_type = 'strong-bisim'
        super(BisimulationTest, self).__init__(name, ymlfile('bisimulation'), settings)
        self.set_command_line_options('t3', ['-e' + equivalence_type])
        self.set_command_line_options('t4', ['-e' + equivalence_type])
        self.set_command_line_options('t7', ['-b' + bisimulation_type])

class LpsConstelmTest(ProcessTest):
    def __init__(self, name, settings = dict()):
        super(LpsConstelmTest, self).__init__(name, ymlfile('lpsconstelm'), settings)
        self.generate_process_parameters = True

class LpsBinaryTest(ProcessTest):
    def __init__(self, name, settings = dict()):
        super(LpsBinaryTest, self).__init__(name, ymlfile('lpsbinary'), settings)
        self.generate_process_parameters = True

class Lps2pbesTest(ProcessTest):
    def __init__(self, name, settings = dict()):
        super(Lps2pbesTest, self).__init__(name, ymlfile('lps2pbes'), settings)

    def create_inputfiles(self, runpath = '.'):
        super(Lps2pbesTest, self).create_inputfiles(runpath)
        self.inputfiles.append(mcrl2file('examples/modal-formulas/nodeadlock.mcf'))

class PbesTest(RandomTest):
    def __init__(self, name, ymlfile, settings = dict()):
        super(PbesTest, self).__init__(name, ymlfile, settings)
        self.equation_count = 4
        self.atom_count = 4
        self.propvar_count = 3
        self.use_quantifiers = True

    def execute(self, runpath = '.'):
        p = make_pbes(self.equation_count, self.atom_count, self.propvar_count, self.use_quantifiers)
        run_pbes_test_with_counter_example_minimization(self.name, self.ymlfile, p, self.settings)

# N.B. does not work yet due to unusable abstraction map
class PbesabsintheTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(PbesabsintheTest, self).__init__(name, ymlfile('pbesabsinthe'), settings)

class PbesabstractTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(PbesabstractTest, self).__init__(name, ymlfile('pbesabstract'), settings)

class PbesconstelmTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(PbesconstelmTest, self).__init__(name, ymlfile('pbesconstelm'), settings)

class PbesparelmTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(PbesparelmTest, self).__init__(name, ymlfile('pbesparelm'), settings)

class PbespareqelmTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(PbespareqelmTest, self).__init__(name, ymlfile('pbespareqelm'), settings)

class PbesrewrTest(PbesTest):
    def __init__(self, name, rewriter, settings = dict()):
        super(PbesrewrTest, self).__init__(name, ymlfile('pbesrewr'), settings)
        self.set_command_line_options('t2', ['-p' + rewriter])

class PbesinstTest(PbesTest):
    def __init__(self, name, options = ['-slazy'], settings = dict()):
        super(PbesinstTest, self).__init__(name, ymlfile('pbesinst'), settings)
        self.set_command_line_options('t2', options)

class PbessolveTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(PbessolveTest, self).__init__(name, ymlfile('pbessolve'), settings)

class PbesstategraphTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(PbesstategraphTest, self).__init__(name, ymlfile('pbesstategraph'), settings)

# N.B does not work due to unknown expressions (F_or)
class SymbolicExplorationTest(PbesTest):
    def __init__(self, name, settings = dict()):
        super(SymbolicExplorationTest, self).__init__(name, ymlfile('symbolic_exploration'), settings)

class BesTest(RandomTest):
    def __init__(self, name, ymlfile, settings = dict()):
        super(BesTest, self).__init__(name, ymlfile, settings)
        self.equation_count = 4
        self.term_size = 3

    def create_inputfiles(self, runpath = '.'):
        filename = '{0}.txt'.format(self.name, self.settings)
        p = make_bes(self.equation_count, self.term_size)
        write_text(filename, str(p))
        self.inputfiles += [filename]

class BessolveTest(BesTest):
    def __init__(self, name, settings = dict()):
        super(BessolveTest, self).__init__(name, ymlfile('bessolve'), settings)

available_tests = {
    'alphabet'                          : lambda name, settings: AlphabetTest(name, settings)                                ,
    'lpssuminst'                        : lambda name, settings: LpsSuminstTest(name, settings)                              ,
    'lpssumelm'                         : lambda name, settings: LpsSumelmTest(name, settings)                               ,
    'lpsparelm'                         : lambda name, settings: LpsParelmTest(name, settings)                               ,
    'lpsconfcheck_commutative'          : lambda name, settings: LpsConfcheckTest(name, 'commutative', settings)             ,
    'lpsconfcheck_commutative-disjoint' : lambda name, settings: LpsConfcheckTest(name, 'commutative-disjoint', settings)    ,
    'lpsconfcheck_disjoint'             : lambda name, settings: LpsConfcheckTest(name, 'disjoint', settings)                ,
    'lpsconfcheck_triangular'           : lambda name, settings: LpsConfcheckTest(name, 'triangular', settings)              ,
    'lpsconfcheck_trivial'              : lambda name, settings: LpsConfcheckTest(name, 'trivial', settings)                 ,
    'lpsconstelm'                       : lambda name, settings: LpsConstelmTest(name, settings)                             ,
    'lpsbinary'                         : lambda name, settings: LpsBinaryTest(name, settings)                               ,
    'lps2pbes'                          : lambda name, settings: Lps2pbesTest(name, settings)                                ,
    'ltscompare_bisim'                  : lambda name, settings: LtscompareTest(name, 'bisim', settings)                     ,
    'ltscompare_branching-bisim'        : lambda name, settings: LtscompareTest(name, 'branching-bisim', settings)           ,
    'ltscompare_dpbranching-bisim'      : lambda name, settings: LtscompareTest(name, 'dpbranching-bisim', settings)         ,
    'ltscompare_weak-bisim'             : lambda name, settings: LtscompareTest(name, 'weak-bisim', settings)                ,
    'ltscompare_dpweak-bisim'           : lambda name, settings: LtscompareTest(name, 'dpweak-bisim', settings)              ,
    'ltscompare_sim'                    : lambda name, settings: LtscompareTest(name, 'sim', settings)                       ,
    'ltscompare_trace'                  : lambda name, settings: LtscompareTest(name, 'trace', settings)                     ,
    'ltscompare_weak-trace'             : lambda name, settings: LtscompareTest(name, 'weak-trace', settings)                ,
    'bisimulation_bisim'                : lambda name, settings: BisimulationTest(name, 'bisim', settings)                   ,
    'bisimulation_branching-bisim'      : lambda name, settings: BisimulationTest(name, 'branching-bisim', settings)         ,
    'bisimulation_weak-bisim'           : lambda name, settings: BisimulationTest(name, 'weak-bisim', settings)              ,
    'pbesabstract'                      : lambda name, settings: PbesabstractTest(name, settings)                            ,
    'pbesconstelm'                      : lambda name, settings: PbesconstelmTest(name, settings)                            ,
    'pbesparelm'                        : lambda name, settings: PbesparelmTest(name, settings)                              ,
    'pbespareqelm'                      : lambda name, settings: PbespareqelmTest(name, settings)                            ,
    'pbesrewr_simplify'                 : lambda name, settings: PbesrewrTest(name, 'simplify', settings)                    ,
    'pbesrewr_pfnf'                     : lambda name, settings: PbesrewrTest(name, 'pfnf', settings)                        ,
    'pbesrewr_quantifier-all'           : lambda name, settings: PbesrewrTest(name, 'quantifier-all', settings)              ,
    'pbesrewr_quantifier-finite'        : lambda name, settings: PbesrewrTest(name, 'quantifier-finite', settings)           ,
    'pbesrewr_quantifier-one-point'     : lambda name, settings: PbesrewrTest(name, 'quantifier-one-point', settings)        ,
    'pbesinst_lazy'                     : lambda name, settings: PbesinstTest(name, ['-slazy'], settings)                    ,
    'pbesinst_alternative_lazy'         : lambda name, settings: PbesinstTest(name, ['-salternative-lazy'], settings)        ,
    'pbesinst_finite'                   : lambda name, settings: PbesinstTest(name, ['-sfinite', '-f*(*:Bool)'], settings)   ,
    'pbessolve'                         : lambda name, settings: PbessolveTest(name, settings)                               ,
    'pbesstategraph'                    : lambda name, settings: PbesstategraphTest(name, settings)                          ,
    'bessolve'                          : lambda name, settings: BessolveTest(name, settings)                                ,
}

def main(tests):
    import argparse
    cmdline_parser = argparse.ArgumentParser()
    cmdline_parser.add_argument('-t', '--toolpath', dest='toolpath', help='The path where the mCRL2 tools are installed')
    cmdline_parser.add_argument('-r', '--repetitions', dest='repetitions', metavar='N', default='10', help='Perform N repetitions of each test')
    cmdline_parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Display additional progress messages.')
    cmdline_parser.add_argument('-k', '--keep-files', dest='keep_files', action='store_true', help='Keep the files produced by the test')
    cmdline_parser.add_argument('-p', '--pattern', dest='pattern', metavar='P', default='.', action='store', help='Run the tests that match with pattern P')
    args = cmdline_parser.parse_args()
    toolpath = args.toolpath
    if not toolpath:
        toolpath = MCRL2_INSTALL_DIR
    settings = {'toolpath': toolpath, 'verbose': args.verbose, 'cleanup_files': not args.keep_files, 'allow-non-zero-return-values': True}
    I = range(int(args.repetitions))

    testdir = 'output'
    if not os.path.exists(testdir):
        os.mkdir(testdir)
    os.chdir(testdir)

    for name in sorted(tests):
        if re.search(args.pattern, name):
            try:
                for i in I:
                    test = tests[name]('{}_{}'.format(name, i), settings)
                    test.execute_in_sandbox()
            except Exception as e:
                print 'Test {} failed!'.format(test.name)
                print e

    # SymbolicExplorationTest('symbolic_exploration', settings)
    # PbesrewrTest('pbesrewr', 'bqnf-quantifier', settings)
    # PbesabstractTest('pbesabsinthe', settings)

if __name__ == '__main__':
    main(available_tests)