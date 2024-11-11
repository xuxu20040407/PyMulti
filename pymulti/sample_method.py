import numpy as np

class Sampler:
    def __init__(self):
        pass

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

        sampling_values = [np.linspace(min_val, max_val, n, endpoint=True) for min_val, max_val, n in
                           zip(min_values, max_values, n_samples)]
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
