import argparse
import cmepy as cme
import cmepy.domain
import cmepy.model
import cmepy.recorder
import cmepy.fsp.solver
import cmepy.fsp.support_expander
import cmepy.domain
import cmepy.statistics
import collections
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import shutil
import subprocess
import scipy.integrate
import scipy.special
import sys
import tempfile
import warnings
from xml.dom import minidom
from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, Comment


import plot_utils


# Name of the StochKit XML file.
_STOCHKIT_FILENAME = 'stochkit.xml'
# Name of folders created by StochKit.
_STOCHKIT_STATS_OUTPUT_MEAN = 'stochkit_output/stats/means.txt'
_STOCHKIT_STATS_OUTPUT_STD = 'stochkit_output/stats/variances.txt'
_STOCHKIT_TRAJECTORIES_OUTPUT = 'stochkit_output/trajectories/trajectory%d.txt'
# Limits the number of explored states.
_MEANODE_STATE_LIMIT = 100000

# Remove extremely low probabilities.
_MEANODE_P_THRESHOLD = 1e-8


# Defines a complex.
class Complex(object):
    def __init__(self, species):
        # Build a dictionary of <species, counts>.
        self.counts = {}
        self.species = sorted(species)
        for s in species:
            if s not in self.counts:
                # We define the CRN species as normal strings.
                self.counts[s] = 0
            self.counts[s] += 1

    def GetSpecies(self):
        species = []
        for k, v in self.counts.iteritems():
            species.extend([k] * v)
        return species

    def __str__(self):
        if self.counts:
            return ' + '.join((('' if c == 1 else str(c)) + s) for s, c in self.counts.iteritems())
        return '(null)'

    def AddXML(self, node, node_type):
        # Type should be either 'Reactants' or 'Products'.
        assert node_type == 'Reactants' or node_type == 'Products'
        subnode = SubElement(node, node_type)
        for s, c in self.counts.iteritems():
            SubElement(subnode, 'SpeciesReference', {'id': s, 'stoichiometry': '%d' % c})
        if not self.counts:
            subnode.append(Comment('no ' + node_type.lower()))


# Defines a reaction.
class Reaction(object):
    def __init__(self, inputs, outputs, rate):
        self.inputs = inputs
        self.outputs = outputs
        self.rate = rate

    def __str__(self):
        return str(self.inputs) + ' -> ' + str(self.outputs) + ('  (rate: %g)' % self.rate)

    def AddXML(self, node, id, ):
        reaction_node = SubElement(node, 'Reaction')
        SubElement(reaction_node, 'Id').text = 'R' + str(id + 1)
        SubElement(reaction_node, 'Description').text = self.__str__()
        SubElement(reaction_node, 'Type').text = 'mass-action'
        SubElement(reaction_node, 'Rate').text = '%g' % self.rate
        self.inputs.AddXML(reaction_node, 'Reactants')
        self.outputs.AddXML(reaction_node, 'Products')

    def SetRate(self, rate):
        self.rate = rate


# Defines a CRN.
class CRN(object):
    def __init__(self, stochkit_binary='ssa'):
        self.species = set()  # Set of strings.
        self.species_descriptions = {}  # Species can have a description (defaults to 'Species X').
        self.species_population = {}  # Species have an initial population (defaults to 0).
        self.reactions = {}   # List of Reactions.
        self.stochkit_binary = stochkit_binary
        self.ordered_species = []  # Always keep a consistent species order.
        self.observable_states = {}  # dict of lists.

    def CreateFromJSONFile(self, filename):
        with open(filename) as fp:
            model_config = json.loads(_RemoveComments(fp.read()))
            assert 'reactions' in model_config, 'Could not find reactions'

        for reaction in model_config['reactions']:
            assert isinstance(reaction['reactants'], list), 'Bad model file'
            assert isinstance(reaction['products'], list), 'Bad model file'
            reactants = [] if 'reactants' not in reaction else reaction['reactants']
            products = [] if 'products' not in reaction else reaction['products']
            assert reactants or products, 'Reaction without any reactant or product found'
            assert 'rate' in reaction, 'Reaction without any rate found'
            assert 'name' in reaction, 'Reaction without any name found'
            rates = []
            if isinstance(reaction['rate'], list):
                assert isinstance(reaction['name'], list), 'When specifying 2 rates, 2 names must be given'
                rates = reaction['rate']
                names = reaction['name']
                if len(rates) == 1:
                    rates.append(None)
                    names.append(None)
                assert len(rates) == 2, 'Reaction with more than 0 or more than 2 rates found'
                assert len(names) == 2, 'Reaction with a different number of names that rates found'
            else:
                rates = [float(reaction['rate']), None]
                names = [reaction['name'], None]
            self.AddReaction(reactants, products, rates[0], names[0])
            if rates[1]:
                self.AddReaction(products, reactants, rates[1], names[1])
        if 'populations' in model_config:
            assert isinstance(model_config['populations'], dict), 'Bad model file'
            for species, count in model_config['populations'].iteritems():
                self.SetInitialPopulation(species, int(count))
        if 'observable_states' in model_config:
            assert isinstance(model_config['observable_states'], dict), 'Bad model file'
            for species, species_to_add in model_config['observable_states'].iteritems():
                self.AddObservableList(species, species_to_add)

    def AddReaction(self, inputs, outputs, rate, name):
        assert name not in self.reactions, 'Reaction with the same name already exists'
        self.species |= set(inputs)
        self.species |= set(outputs)
        self.ordered_species = sorted(self.species)
        self.reactions[name] = Reaction(Complex(inputs), Complex(outputs), rate)

    def SetSpeciesDescription(self, species, description):
        assert species in self.species, 'No reactions with species "%s". Description cannot be set.' % species
        self.species_descriptions[species] = description

    def SetInitialPopulation(self, species, population):
        assert species in self.species, 'No reactions with species "%s". Population cannot be set.' % species
        self.species_population[species] = population

    def AddObservableList(self, name, species):
        for s in species:
            assert s in self.species, 'No reactions with species "%s". Observable state cannot be added.' % species
        self.observable_states[name] = set(species)

    def SetReactionRate(self, name, rate):
        assert name in self.reactions, 'Reaction with name "%s" not found. Cannot change rate' % name
        self.reactions[name].SetRate(rate)

    def HasSpecies(self, name):
        return name in self.species

    def HasReaction(self, name):
        return name in self.reactions

    def GetReactionRate(self, name):
        return self.reactions[name].rate

    def NumBaseSpecies(self):
        return len(self.species_population)

    def GetBasePopulations(self):
        p = []
        for s in self.ordered_species:
            if s in self.species_population:
                p.append(self.species_population[s])
        return np.array(p)

    def SetBasePopulations(self, values):
        index = 0
        for s in self.ordered_species:
            if s in self.species_population:
                self.SetInitialPopulation(s, values[index])
                index += 1

    def GetGiniCoefficientOfRates(self):
        rates = np.array([r.rate for r in self.reactions.values()]).reshape(len(self.reactions), 1)
        return np.sum(np.abs(rates - rates.T)) / (2. * np.sum(rates) * rates.shape[0])

    def GetReactionRates(self, regexp=r'.*'):
        return np.array([v.rate for k, v in sorted(self.reactions.iteritems()) if re.match(regexp, k)])

    def GetReactions(self, regexp=r'.*'):
        return np.array([v for k, v in sorted(self.reactions.iteritems()) if re.match(regexp, k)])

    def __str__(self):
        output = 'Reactions:\n\n' + '\n'.join(str(r) for r in self.reactions.values())
        if self.species_population:
            output += ('\n\nInitial populations:\n\n' +
                       '\n'.join((s + ' = ' + str(c)) for s, c in self.species_population.iteritems()))
        if self.observable_states:
            output += ('\n\nObservable states:\n\n' +
                       '\n'.join((n + ' = ' + ' + '.join(s)) for n, s in self.observable_states.iteritems()))
        return output

    def _XMLString(self, description):
        # This function creates a StochKit compatible XML code.
        model_node = Element('Model')
        model_node.append(Comment('Generated by crn.py'))
        SubElement(model_node, 'Description').text = description
        # General info.
        SubElement(model_node, 'NumberOfReactions').text = str(len(self.reactions))
        SubElement(model_node, 'NumberOfSpecies').text = str(len(self.species))
        # Reactions.
        reactionslist_node = SubElement(model_node, 'ReactionsList')
        for i, reaction in enumerate(self.reactions.values()):
            reaction.AddXML(reactionslist_node, i)
        # Species.
        specieslist_node = SubElement(model_node, 'SpeciesList')
        for i, species in enumerate(self.ordered_species):
            species_node = SubElement(specieslist_node, 'Species')
            SubElement(species_node, 'Id').text = species
            SubElement(species_node, 'Description').text = (
                self.species_descriptions[species] if species in self.species_descriptions else
                ('Species %s' % species))
            SubElement(species_node, 'InitialPopulation').text = (
                '%d' % (self.species_population[species] if species in self.species_population else 0))
        # Small trick to make the XML pretty.
        raw_output = ElementTree.tostring(model_node, 'utf-8')
        return minidom.parseString(raw_output).toprettyxml(indent="  ")

    # Runs the complete simulation using StochKit and return 3 values:
    # (a) The list of species
    # (b) The timestamps at which the simulation is snapshot
    # (c) The species population at these timestamps as a 3D tensor with dimensions nruns x ndatapoints x nspecies.
    def Simulate(self, duration=10.0, nruns=1, ndatapoints=None, output_directory=None):
        # Create XML file (overwrite old file if any).
        # If output_directory is None, create a temporary directory.
        cleandir = False
        if output_directory is None:
            output_directory = tempfile.mkdtemp()
            cleandir = True
        xmlpath = os.path.join(output_directory, _STOCHKIT_FILENAME)
        with open(xmlpath, 'w') as fp:
            fp.write(self._XMLString('Autogenerated CRN'))
        # Keep one datapoint per time unit if ndatapoints is not set.
        if ndatapoints is None:
            ndatapoints = int(duration) + 1
        # Prepare command.
        commandline = [
            self.stochkit_binary,
            '-m', xmlpath,
            '-t', '%g' % duration,
            '-r', '%d' % nruns,
            '-i', '%d' % (ndatapoints - 1),
            '--keep-trajectories',
            '-f']
        print 'Executing:', ' '.join(commandline)
        if subprocess.call(commandline) != 0:
            raise RuntimeError('Error occur while running StochKit')
        # Gather trajectories back into numpy arrays.
        timestamps = np.empty((ndatapoints, 1))
        data = np.empty((nruns, ndatapoints, len(self.species)))
        for run in xrange(nruns):
            run_filename = os.path.join(output_directory, _STOCHKIT_TRAJECTORIES_OUTPUT % run)
            d = np.loadtxt(run_filename)
            data[run, :, :] = d[:, 1:]
            if run == 0:
                timestamps = d[:, 0]
        if cleandir:
            shutil.rmtree(output_directory)
        return self.ordered_species, timestamps, data

    def _CMEModel(self, description):
        # Utilities lambda generators.
        def _SpeciesCountsFunction(i):
            return lambda *x: x[i]

        def _PropensitiesFunction(reaction, species_counts, species_to_index):
            def action_mass(*x):
                r = reaction.rate
                for s, c in reaction.inputs.counts.iteritems():
                    r *= species_counts[species_to_index[s]](*x) ** float(c)
                return r
            return action_mass

        def _GetTransitions(reaction, species_to_index):
            transitions = [0]*len(species_to_index)
            for s, c in reaction.inputs.counts.iteritems():
                transitions[species_to_index[s]] -= c
            for s, c in reaction.outputs.counts.iteritems():
                transitions[species_to_index[s]] += c
            return transitions

        species_to_index = dict((s, i) for i, s in enumerate(self.ordered_species))
        species_counts = [_SpeciesCountsFunction(i) for i in xrange(len(self.ordered_species))]
        reactions = []
        propensities = []
        transitions = []
        for i, reaction in enumerate(self.reactions.values()):
            reactions.append(str(reaction))
            propensities.append(_PropensitiesFunction(reaction, species_counts, species_to_index))
            transitions.append(_GetTransitions(reaction, species_to_index))
        initial_state = [0]*len(species_to_index)
        for s, c in self.species_population.iteritems():
            initial_state[species_to_index[s]] = c
        initial_state = tuple(initial_state)
        return cme.model.create(name=description, species=self.ordered_species, species_counts=species_counts,
                                reactions=reactions, propensities=propensities, transitions=transitions,
                                initial_state=initial_state)

    # Runs the complete simulation using CMEPy and return 2 values:
    # (a) The list of species
    # (b) The cmepy.recorder object. You can use BuildDistribution and the Plotter utility on it.
    # See http://fcostin.github.io/cmepy/modules/recorder.html#module-recorder for more details on it.
    def CME(self, time_steps, distribution_precision=1e-2, verbose=False):
        model = self._CMEModel('Autogenerated CRN')
        initial_states = cme.domain.from_iter((model.initial_state, ))
        # Create expander for FSP expansion strategy.
        expander = cme.fsp.support_expander.SupportExpander(
            model.transitions, depth=1, epsilon=1.e-5)
        # Create FSP solver/
        fsp_solver = cme.fsp.solver.create(model, initial_states, expander)
        # Error of the solution at the final time to be bounded.
        num_steps = np.size(time_steps)
        max_error_per_step = distribution_precision / float(num_steps)
        # Create recorder to record species counts.
        recorder = cmepy.recorder.create((model.species, model.species_counts))
        print 'Executing: CME'
        for i, t in enumerate(time_steps):
            if i % max(len(time_steps) / 20, 1) == 0 and verbose:
                print '%d%%' % (i * 100 / len(time_steps))
            fsp_solver.step(t, max_error_per_step)
            # Record the solution at this timestep.
            recorder.write(t, fsp_solver.y[0])
        print 'done!'
        return self.ordered_species, recorder

    # Runs a complete simulation but only estimates the final distribution.
    # It assumes that the system is complex balanced for the distribution. The mean is correct.
    # Returns:
    # (a) The list of species.
    # (b) The timestamps at which the simulation is snapshot.
    # (c) The average species population at these timestamps as a matrix with dimensions ndatapoints x (nspecies + nspecies^2).
    # (d) The final distribution as a dictionary of all possible population states.
    def Average(self, time_steps):
        def f(x, R, P, nspecies):
            dx = np.zeros_like(x)
            for r, p in zip(R, P):
                dx += p(*x) * r
            return dx

        # Compute steady-state (takes about 30ms for 5 species).
        model = self._CMEModel('Autogenerated CRN')
        cmepy.model.validate_model(model)
        R = [np.array(p) for p in model.transitions]
        nspecies = len(self.ordered_species)
        initial_state = model.initial_state
        dx_dt = lambda x, t: f(x, R, model.propensities, nspecies)
        mean_x = scipy.integrate.odeint(dx_dt, initial_state, time_steps)

        # Find all possible states (takes roughly 500ms for 20000 states).
        start = tuple(model.initial_state)
        stack = collections.deque([start])
        explored_states = set([start])
        while stack:
            assert len(explored_states) < _MEANODE_STATE_LIMIT, 'CRN is too complex or has too many individuals. Try reducing the initial population'
            current_state = stack.pop()
            current_state_array = np.array(current_state)
            for r in R:
                next_state = tuple(current_state_array + r)
                if (all(s >= 0 for s in next_state) and
                    next_state not in explored_states):
                    explored_states.add(next_state)
                    stack.append(next_state)

        # Compute distribution as if the system was complex balanced (takes roughly 450ms for 20000 states).
        ln_steady_state = np.log(mean_x[-1, :])
        gammaln_state = {}  # Cache Gamma computation.
        ln_distribution = {}
        sum_p = 0.
        for current_state in explored_states:
            ln_p = 0.
            for ln_ci, xi in zip(ln_steady_state, current_state):
                if xi not in gammaln_state:
                    gammaln_state[xi] = scipy.special.gammaln(xi + 1)
                # The following is mathematically equivalent:
                # for i in xrange(1, int(xi) + 1):
                #     ln_p += np.log(ci) - np.log(float(i))
                # To:
                #   ln_p += xi * np.log(ci) - scipy.special.gammaln(xi + 1)
                ln_p += xi * ln_ci - gammaln_state[xi]
            ln_distribution[current_state] = ln_p
            sum_p += np.exp(ln_p)
        distribution = {}
        for current_state in ln_distribution:
            p = np.exp(ln_distribution[current_state] - np.log(sum_p))
            if p > _MEANODE_P_THRESHOLD:
                distribution[current_state] = p

        return self.ordered_species, time_steps, mean_x, distribution

    def SteadyState(self, time_steps):
        def f(x, R, P, nspecies):
            dx = np.zeros_like(x)
            for r, p in zip(R, P):
                dx += p(*x) * r
            return dx

        # Compute steady-state (takes about 30ms for 5 species).
        model = self._CMEModel('Autogenerated CRN')
        cmepy.model.validate_model(model)
        R = [np.array(p) for p in model.transitions]
        nspecies = len(self.ordered_species)
        initial_state = model.initial_state
        dx_dt = lambda x, t: f(x, R, model.propensities, nspecies)
        mean_x = scipy.integrate.odeint(dx_dt, initial_state, time_steps)
        return self.ordered_species, mean_x[-1, :]

    # Convert an output of run to the list of observable states.
    def ConvertToObservable(self, output_from_run):
        # Helper function for sum relevant conditions in each aggregate.
        def SpeciesCountsFunction(aggregate):
            return lambda *x: reduce(lambda v, w: v + x[w], aggregate, 0)
        assert self.observable_states, 'No observable state was specified'
        ordered_observable_states = list(self.observable_states.iteritems())
        species_to_observable_states = {}
        for s in self.species:
            for i, (_, species_set) in enumerate(ordered_observable_states):
                if s in species_set:
                    assert s not in species_to_observable_states, 'Observable state lists are not disjoint'
                    species_to_observable_states[s] = i
        aggregated_species = [k for k, _ in ordered_observable_states]

        if len(output_from_run) == 3:
            # StochKit data.
            species, timestamps, data = output_from_run
            aggregated_data = np.zeros((data.shape[0], data.shape[1], len(self.observable_states)))
            for i, s in enumerate(species):
                if s in species_to_observable_states:
                    aggregated_data[:, :, species_to_observable_states[s]] += data[:, :, i]
            return aggregated_species, timestamps, aggregated_data

        elif len(output_from_run) == 2:
            # CME data.
            species, recorder = output_from_run
            # List all observable species.
            aggregates = [[] for _ in xrange(len(self.observable_states))]
            for i, s in enumerate(species):
                if s in species_to_observable_states:
                    aggregates[species_to_observable_states[s]].append(i)
            new_recorder = cmepy.recorder.create((aggregated_species, [SpeciesCountsFunction(a) for a in aggregates]))
            species_tuple = tuple(species)
            measurement = recorder[species_tuple]
            time_steps = measurement.times
            for i, t in enumerate(time_steps):
                new_recorder.write(t, measurement.distributions[i])
            return aggregated_species, new_recorder

        elif len(output_from_run) == 4:
            # MeanODE data.
            species, timestamps, mean_x, distribution = output_from_run
            aggregated_mean_x = np.zeros((mean_x.shape[0], len(self.observable_states)))
            for i, s in enumerate(species):
                if s in species_to_observable_states:
                    aggregated_mean_x[:, species_to_observable_states[s]] += mean_x[:, i]
            aggregated_distribution = {}
            for k, v in distribution.iteritems():
                new_k = [0] * len(self.observable_states)
                for i, s in enumerate(species):
                    new_k[species_to_observable_states[s]] += k[i]
                new_k = tuple(new_k)
                if new_k not in aggregated_distribution:
                    aggregated_distribution[new_k] = 0.
                aggregated_distribution[new_k] += v
            return aggregated_species, timestamps, aggregated_mean_x, aggregated_distribution

        return None


# This function builds an dict. This dict can be
# used to compare two steady-state distributions.
def BuildDistribution(output_from_run, from_timestamp=0.):
    if len(output_from_run) == 3:
        # StochKit.
        counts = {}
        species, timestamps, trajectories = output_from_run
        trajectories = trajectories[:, timestamps >= from_timestamp, :]
        nruns = trajectories.shape[0]
        ntimestamps = trajectories.shape[1]
        total_counts = 0
        for i in xrange(nruns):
            for j in xrange(ntimestamps):
                indices = tuple(trajectories[i, j, :].astype(int).tolist())
                if indices not in counts:
                    counts[indices] = 0
                counts[indices] += 1
                total_counts += 1
        for k, v in counts.iteritems():
            counts[k] = float(counts[k]) / float(total_counts)
        return counts
    elif len(output_from_run) == 2:
        # CMEPy.
        species, recorder = output_from_run
        return recorder[tuple(species)].distributions[-1]
    elif len(output_from_run) == 4:
        return output_from_run[-1]
    raise ValueError('Unable to detect the type of data.')


def _CreateDirectory(directory, force=False):
    if force and os.path.exists(directory):
        print 'Deleting previous directory %s.' % directory
        shutil.rmtree(directory)
    print 'Preparing directory %s.' % directory
    try:
        os.makedirs(directory)
    except os.error:
        print 'Cannot create directory %s. Make sure it does not exist (or use --force).' % directory
        return False
    return True


def _RemoveComments(text):
    def Replacer(match):
        s = match.group(0)
        if s.startswith('/'):
            return ' '
        else:
            return s
    pattern = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE)
    return re.sub(pattern, Replacer, text)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs StochKit/CMEPy on a CRN')
    parser.add_argument('--model', metavar='PATH', action='store', required=True, help='Path where the model is stored')
    parser.add_argument('--directory', metavar='PATH', action='store', default='/tmp/run_system_output', help='Path where the data is stored')
    parser.add_argument('--stochkit_path', metavar='PATH', action='store', default='./StochKit2.0.11/ssa', help='Path where the stochkit binary is stored')
    parser.add_argument('--nruns', metavar='N', type=int, action='store', default=100, help='Number of simulation runs')
    parser.add_argument('--duration', metavar='SECONDS', type=float, action='store', default=20., help='Duration in seconds of each simulation run')
    parser.add_argument('--force', action='store_true', help='If set, the directory is overwritten')
    parser.add_argument('--show_plots', action='store_true', help='If set, plots are shown after the simulations')
    parser.add_argument('--cmepy', action='store_true', help='If set, uses CMEPy instead of StochKit')
    parser.add_argument('--meanode', action='store_true', help='If set, uses mean-field ODE instead of StochKit or CMEPy')
    parser.add_argument('--integration_dt', metavar='SECONDS', type=float, action='store', default=0.01, help='dt used between integration steps.')
    parser.add_argument('--show_observable_plots', action='store_true', help='If set, plots are shown after the simulations with observable data aggregated')
    args = parser.parse_args()

    if not (args.cmepy or args.meanode) and not _CreateDirectory(args.directory, args.force):
        sys.exit(-1)
    print ''

    crn = CRN(stochkit_binary=args.stochkit_path)
    crn.CreateFromJSONFile(args.model)

    print 'Running simulation on:\n----------------------\n'
    print crn, '\n'

    time_steps = np.concatenate((
        np.arange(0., args.duration, args.integration_dt),
        np.array([args.duration])  # Last timestep.
    ))
    if args.cmepy:
        try:
            output = crn.CME(time_steps)
        except ValueError:
            print 'Integration error detected. Try reducing --integration_dt'
            sys.exit(-1)
    elif args.meanode:
        output = crn.Average(time_steps)
    else:
        try:
            output = crn.Simulate(duration=args.duration, nruns=args.nruns, output_directory=args.directory)
        except OSError:
            print 'Cannot find StochKit (%s)' % args.stochkit_path
            sys.exit(-1)

    if args.show_plots:
        plotter = plot_utils.Plotter(output)
        # Plot average.
        plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions(from_timestamp=args.duration / 2.0)

    # Show plots of the observable data.
    if args.show_observable_plots:
        aggregated_output = crn.ConvertToObservable(output)
        plotter = plot_utils.Plotter(aggregated_output)
        # Plot average.
        plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions(from_timestamp=args.duration / 2.0)

    if args.show_plots or args.show_observable_plots:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            plt.show()
