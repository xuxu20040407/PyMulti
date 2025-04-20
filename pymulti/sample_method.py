from calendar import c
import numpy as np


class Sampler:
    def __init__(self):
        pass

    import numpy as np

    def uniform_sampling(self, max_values, min_values, n_samples):
        """
        Generate uniform sampling grids for multiple dimensions.

        Parameters:
        max_values (list or array): Maximum values for each dimension.
        min_values (list or array): Minimum values for each dimension.
        n_samples (int or list or array): Number of samples for each dimension.
            If an integer, the same number of samples will be used for all dimensions.
            If a list or array, it must have the same length as max_values and min_values.

        Returns:
        sampled_grid (numpy.ndarray): A 2D array containing the sampled grid points.
        """
        if not isinstance(max_values, (list, np.ndarray)) or not isinstance(min_values, (list, np.ndarray)):
            raise ValueError("max_values and min_values must be lists or numpy arrays")
        if isinstance(n_samples, int):
            n_samples = [n_samples] * len(max_values)
        elif not isinstance(n_samples, (list, np.ndarray)) or len(n_samples) != len(max_values):
            raise ValueError(
                "n_samples must be an integer or a list/array with the same length as max_values and min_values")

        # 检查每个维度的max和min值是否相等，并据此生成采样值
        sampling_values = []
        for min_val, max_val, n in zip(min_values, max_values, n_samples):
            if min_val == max_val:
                # 如果最大值和最小值相等，只生成一个值
                sampling_values.append(np.array([min_val]))
            else:
                # 否则，生成从min_val到max_val的n个样本点
                sampling_values.append(np.linspace(min_val, max_val, n, endpoint=True))

        grids = np.meshgrid(*sampling_values, indexing='ij')
        flattened_grids = [grid.ravel() for grid in grids]
        sampled_grid = np.column_stack(flattened_grids)
        return sampled_grid

    def random_sampling(self, max_values, min_values, n_samples):
        """
        Generate random sampling grids for multiple dimensions.

        Parameters:
        max_values (list or array): Maximum values for each dimension.
        min_values (list or array): Minimum values for each dimension.
        n_samples (int or list or array): Number of samples for each dimension.
            If an integer, the same number of samples will be used for all dimensions.
            If a list or array, it must have the same length as max_values and min_values.

        Returns:
        sampled_grid (numpy.ndarray): A 2D array containing the sampled grid points.
        """
        if not isinstance(max_values, (list, np.ndarray)) or not isinstance(min_values, (list, np.ndarray)):
            raise ValueError("max_values and min_values must be lists or numpy arrays")
        if isinstance(n_samples, int):
            n_samples = [n_samples] * len(max_values)
        elif not isinstance(n_samples, (list, np.ndarray)) or len(n_samples) != len(max_values):
            raise ValueError(
                "n_samples must be an integer or a list/array with the same length as max_values and min_values")

        # Generate random sampling values for each dimension
        sampled_values = [np.random.uniform(min_val, max_val, size=n) for min_val, max_val, n in
                          zip(min_values, max_values, n_samples)]

        # Stack the sampled values
        sampled_grid = np.column_stack(sampled_values)

        return sampled_grid



class MC_Sampler:
    """
    A class for performing Markov Chain sampling to optimize a target function.

    Attributes:
    min_values (numpy.ndarray): The minimum values for each dimension.
    max_values (numpy.ndarray): The maximum values for each dimension.
    step_sizes (numpy.ndarray): The step sizes for each dimension, calculated as a ratio of the range.
    current_sample (numpy.ndarray): The current sample point in the Markov Chain.
    samples (list): A list of all sample points generated during the sampling process.
    sample_values (list): A list of function values corresponding to each sample point.
    accepted_samples (list): A list indicating whether each sample was accepted (1) or not (0).
    optimization_direction (str): The direction of optimization, either 'max' for maximization or 'min' for minimization.
    T (numpy.ndarray): A modified number for accepted probability. T=10 is easily accepted and T=5 is default.
    """

    def __init__(self, min_values, max_values, step_size_ratio=0.1, optimization_direction='max',T=5):
        self.min_values = np.array(min_values)
        self.max_values = np.array(max_values)
        self.step_sizes = (self.max_values - self.min_values) * step_size_ratio
        self.current_sample = None
        self.samples = []
        self.sample_values = []
        self.accepted_samples = []
        self.optimization_direction = optimization_direction  # 'max' for maximization, 'min' for minimization
        self.T=T

    def propose_sample(self, current_sample):
        """
        Generate a new sample point by perturbing the current sample with a normal distribution.

        Parameters:
        current_sample (numpy.ndarray): The current sample point.

        Returns:
        new_sample (numpy.ndarray): A new proposed sample point.
        """
        new_sample = current_sample + np.random.normal(0, self.step_sizes)
        new_sample = np.clip(new_sample, self.min_values, self.max_values)
        return new_sample

    def is_sample_saved(self, new_sample):
        """
        Check if a sample has already been generated in the sampling process.

        Parameters:
        new_sample (numpy.ndarray): The new sample point to check.

        Returns:
        bool: True if the sample has been saved, False otherwise.
        """
        return tuple(new_sample) in [tuple(sample) for sample in self.samples]

    def get_saved_value(self, new_sample):
        """
        Get the function value of a previously saved sample.

        Parameters:
        new_sample (numpy.ndarray): The sample point for which to retrieve the function value.

        Returns:
        float: The function value of the sample point.
        """
        for i, sample in enumerate(self.samples):
            if np.array_equal(sample, new_sample):
                return self.sample_values[i]

    def MC_sampling(self, target_function):
        if self.current_sample is None:
            self.current_sample = np.random.uniform(self.min_values, self.max_values)
            self.samples.append(self.current_sample)
            self.sample_values.append(target_function(self.current_sample))

        new_sample = self.propose_sample(self.current_sample)
        if self.is_sample_saved(new_sample):
            new_value = self.get_saved_value(new_sample)
        else:
            new_value = target_function(new_sample)
            self.samples.append(new_sample)
            self.sample_values.append(new_value)

        current_value = self.sample_values[-2]
        if self.optimization_direction == 'max':
            acceptance_probability = min(1, np.exp((new_value - current_value)/current_value*self.T))  # For maximization
        elif self.optimization_direction == 'min':
            acceptance_probability = min(1, np.exp((current_value - new_value)/current_value*self.T))  # For minimization
        print('acceptance_probability: ',acceptance_probability)
        if np.random.rand() < acceptance_probability:
            self.current_sample = new_sample
            self.accepted_samples.append(1) # Mark the sample as accepted
            print('accepted')
        else:
            self.accepted_samples.append(0) # Mark the sample as not accepted
            print('unaccepted')