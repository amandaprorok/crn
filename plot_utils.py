import matplotlib.pyplot as plt
import numpy as np


# Utility class that takes the results from a CRN::Run() and
# shows different types of plots.
class Plotter(object):
    _STOCHKIT = 'StochKit'
    _CMEPY = 'CMEPy'
    _MEANODE = 'MeanODE'

    def __init__(self, output_from_run, exclude_list=()):
        if len(output_from_run) == 3:
            self.type = Plotter._STOCHKIT
            self.species, self.timestamps, self.trajectories = output_from_run
            print 'Plotter: Assuming StochKit data.'
        elif len(output_from_run) == 2:
            self.type = Plotter._CMEPY
            self.species, self.recorder = output_from_run
            ts = self.recorder[self.species[0]].times
            self.timestamps = np.array(ts)
            print 'Plotter: Assuming CMEPy data.'
        elif len(output_from_run) == 4:
            self.type = Plotter._MEANODE
            self.species = output_from_run[0]
            self.timestamps = output_from_run[1]
            self.mean_trajectory = output_from_run[2]
            self.last_distribution = output_from_run[3]
            print 'Plotter: Assuming MeanODE data.'
        self.nspecies = len(self.species)
        self.exclude_list = exclude_list

    def AverageTrajectory(self):
        if self.type == Plotter._STOCHKIT:
            mean_trajectory = np.mean(self.trajectories, axis=0)  # First axis is runs.
            std_trajectory = np.std(self.trajectories, axis=0)
        elif self.type == Plotter._CMEPY:
            mean_trajectory = np.empty((self.timestamps.shape[0], self.nspecies))
            std_trajectory = np.empty((self.timestamps.shape[0], self.nspecies))
            for i, s in enumerate(self.species):
                measurement = self.recorder[s]
                mean_trajectory[:, i] = measurement.expected_value
                std_trajectory[:, i] = measurement.standard_deviation
        elif self.type == Plotter._MEANODE:
            mean_trajectory = self.mean_trajectory[:, :self.nspecies]
            std_trajectory = np.zeros((self.timestamps.shape[0], self.nspecies))
        colors = _GetColors(self.nspecies)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in xrange(self.nspecies):
            if self.species[i] in self.exclude_list:
                continue
            plt.fill_between(self.timestamps,
                             mean_trajectory[:, i] + std_trajectory[:, i],
                             mean_trajectory[:, i] - std_trajectory[:, i],
                             facecolor=colors[i], alpha=0.3)
        for i in xrange(self.nspecies):
            if self.species[i] in self.exclude_list:
                continue
            plt.plot(self.timestamps, mean_trajectory[:, i], linewidth=2, color=colors[i], label=self.species[i])
        plt.legend(loc='upper right', shadow=False)
        plt.xlabel('Time [s]')
        plt.ylabel('Population')
        return fig, ax

    def Distributions(self, from_timestamp=0.):
        indices = self.timestamps >= from_timestamp  # Not used for CMEPY.
        if self.type == Plotter._STOCHKIT:
            min_population = float('inf')
            max_population = -1.
            for i in xrange(self.nspecies):
                if self.species[i] in self.exclude_list:
                    continue
                min_population = min(min_population, np.min(self.trajectories[:, indices, i]))
                max_population = max(max_population, np.max(self.trajectories[:, indices, i]))
        elif self.type == Plotter._CMEPY:
            min_population = float('inf')
            max_population = -1.
            for s in self.species:
                if s in self.exclude_list:
                    continue
                for k in self.recorder[s].distributions[-1]:
                    max_population = max(k[0], max_population)
                    min_population = min(k[0], min_population)
        elif self.type == Plotter._MEANODE:
            min_population = float('inf')
            max_population = -1.
            for i in xrange(self.nspecies):
                if self.species[i] in self.exclude_list:
                    continue
                for k, v in self.last_distribution.iteritems():
                    if v > 1e-3:
                        max_population = max(k[i], max_population)
                        min_population = min(k[i], min_population)
        colors = _GetColors(self.nspecies)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in xrange(self.nspecies):
            if self.species[i] in self.exclude_list:
                continue
            if self.type == Plotter._STOCHKIT:
                populations = self.trajectories[:, indices, i].flatten().astype(int)
                plt.hist(populations, normed=True, bins=np.arange(min_population - 0.5, max_population + 0.5, 1.),
                         label=self.species[i], color=colors[i], alpha=0.3)
            elif self.type == Plotter._CMEPY:
                values = np.zeros((max_population - min_population + 1, ))
                for k, v in self.recorder[self.species[i]].distributions[-1].iteritems():
                    values[k[0] - min_population] = v
                plt.bar(np.arange(min_population - 0.5, max_population - 0.4, 1.), values, width=1.,
                        label=self.species[i], color=colors[i], alpha=0.3)
            elif self.type == Plotter._MEANODE:
                values = np.zeros((max_population - min_population + 1, ))
                for k, v in self.last_distribution.iteritems():
                    if v > 1e-3:
                        values[k[i] - min_population] += v
                plt.bar(np.arange(min_population - 0.5, max_population - 0.4, 1.), values, width=1.,
                        label=self.species[i], color=colors[i], alpha=0.3)
        plt.legend(loc='upper right', shadow=False)
        plt.xlabel('Population')
        plt.ylabel('Probability')
        return fig, ax

def _GetColors(n):
    cm = plt.get_cmap('gist_rainbow')
    return [cm(float(i) / float(n)) for i in range(n)]
