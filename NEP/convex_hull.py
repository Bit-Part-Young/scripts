"""
ConvexHull 类实现
The source code is copied from `tools` module of icet package.
"""

import itertools
from typing import List, Sized, Union

import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import ConvexHull as ConvexHullSciPy
from scipy.spatial import QhullError


class ConvexHull:
    """This class provides functionality for extracting the convex hull
    of the (free) energy of mixing. It is based on the `convex hull
    calculator in SciPy
    <http://docs.scipy.org/doc/scipy-dev/reference/\
generated/scipy.spatial.ConvexHull.html>`_.

    Parameters
    ----------
    concentrations
        Concentrations for each structure listed as ``[[c1, c2], [c1, c2],
        ...]``; for binaries, in which case there is only one independent
        concentration, the format ``[c1, c2, c3, ...]`` works as well.
    energies
        Energy (or energy of mixing) for each structure.

    Attributes
    ----------
    concentrations : np.ndarray
        Concentrations of the structures on the convex hull.
    energies : np.ndarray
        Energies of the structures on the convex hull.
    dimensions : int
        Number of independent concentrations needed to specify a point in
        concentration space (1 for binaries, 2 for ternaries and so on).
    structures : List[int]
        Indices of structures that constitute the convex hull (indices are
        defined by the order of their concentrations and energies are fed when
        initializing the :class:`ConvexHull` object).

    Examples
    --------
    A :class:`ConvexHull` object is easily initialized by providing lists of
    concentrations and energies::

        >>> data = {'concentration': [0,    0.2,  0.2,  0.3,  0.4,  0.5,  0.8,  1.0],
        ...         'mixing_energy': [0.1, -0.2, -0.1, -0.2,  0.2, -0.4, -0.2, -0.1]}
        >>> hull = ConvexHull(data['concentration'], data['mixing_energy'])

    Now one can for example access the points along the convex hull directly::

        >>> for c, e in zip(hull.concentrations, hull.energies):
        ...     print(c, e)
        0.0 0.1
        0.2 -0.2
        0.5 -0.4
        1.0 -0.1

    or plot the convex hull along with the original data using e.g., matplotlib::

        >>> import matplotlib.pyplot as plt
        >>> plt.scatter(data['concentration'], data['mixing_energy'], color='darkred')
        >>> plt.plot(hull.concentrations, hull.energies)
        >>> plt.show(block=False)

    It is also possible to extract structures at or close to the convex hull::

        >>> low_energy_structures = hull.extract_low_energy_structures(
        ...     data['concentration'], data['mixing_energy'],
        ...     energy_tolerance=0.005)

    A complete example can be found in the :ref:`basic tutorial
    <tutorial_enumerate_structures>`.
    """

    def __init__(
        self,
        concentrations: Union[List[float], List[List[float]]],
        energies: List[float],
    ) -> None:
        assert len(concentrations) == len(energies)
        # Prepare data in format suitable for SciPy-ConvexHull
        concentrations = np.array(concentrations)
        energies = np.array(energies)
        points = np.column_stack((concentrations, energies))
        self.dimensions = len(points[0]) - 1

        # Construct convex hull
        hull = ConvexHullSciPy(points)

        # Collect convex hull points in handy arrays
        concentrations = []
        energies = []
        for vertex in hull.vertices:
            if self.dimensions == 1:
                concentrations.append(points[vertex][0])
            else:
                concentrations.append(points[vertex][0:-1])
            energies.append(points[vertex][-1])
        concentrations = np.array(concentrations)
        energies = np.array(energies)

        structures = hull.vertices
        # If there is just one independent concentration, we'd better sort
        # according to it
        if self.dimensions == 1:
            ces = list(zip(*sorted(zip(concentrations, energies, structures))))
            self.concentrations = np.array(ces[0])
            self.energies = np.array(ces[1])
            self.structures = np.array(ces[2])
        else:
            self.concentrations = concentrations
            self.energies = energies
            self.structures = structures

        # Remove points that are above the "pure components plane"
        self._remove_points_above_tie_plane()

    def __str__(self) -> str:
        width = 40
        s = []
        s += ["{s:=^{n}}".format(s=" Convex Hull ", n=width)]
        s += [" {:24} : {}".format("dimensions", self.dimensions)]
        s += [" {:24} : {}".format("number of points", len(self.concentrations))]
        s += [
            " {:24} : {}".format("smallest concentration", np.min(self.concentrations))
        ]
        s += [
            " {:24} : {}".format("largest concentration", np.max(self.concentrations))
        ]
        s += ["".center(width, "=")]
        return "\n".join(s)

    def _repr_html_(self) -> str:
        """HTML representation. Used, e.g., in jupyter notebooks."""
        s = ["<h4>Convex Hull</h4>"]
        s += ['<table border="1" class="dataframe">']
        s += ["<tbody>"]
        s += [
            f'<tr><td style="text-align: left;">Dimensions</td><td>{self.dimensions}</td></tr>'
        ]
        s += [
            '<tr><td style="text-align: left;">Number of points</td>'
            f"<td>{len(self.concentrations)}</td></tr>"
        ]
        s += [
            '<tr><td style="text-align: left;">Smallest concentration</td>'
            f"<td>{np.min(self.concentrations)}</td></tr>"
        ]
        s += [
            '<tr><td style="text-align: left;">Largest concentration</td>'
            f"<td>{np.max(self.concentrations)}</td></tr>"
        ]
        s += ["</tbody>"]
        s += ["</table>"]
        return "".join(s)

    def _remove_points_above_tie_plane(self, tol: float = 1e-6) -> None:
        """
        Remove all points on the convex hull that correspond to maximum rather
        than minimum energy.

        Parameters
        ----------
        tol
            Tolerance for what energy constitutes a lower one.
        """

        # Identify the "complex concentration hull", i.e. the extremal
        # concentrations. In the simplest case, these should simply be the
        # pure components.
        if self.dimensions == 1:
            # Then the ConvexHullScipy function doesn't work, so we just pick
            # the indices of the lowest and highest concentrations.
            vertices = []
            vertices.append(np.argmin(self.concentrations))
            vertices.append(np.argmax(self.concentrations))
            vertices = np.array(vertices)
        else:
            concentration_hull = ConvexHullSciPy(self.concentrations)
            vertices = concentration_hull.vertices

        # Remove all points of the convex energy hull that have an energy that
        # is higher than what would be gotten with pure components at the same
        # concentration. These points are mathematically on the convex hull,
        # but in the physically uninteresting upper part, i.e. they maximize
        # rather than minimize energy.
        to_delete = []
        for i, concentration in enumerate(self.concentrations):
            # The points on the convex concentration hull should always be
            # included, so skip them.
            if i in vertices:
                continue

            # The energy obtained as a linear combination of concentrations on
            # the convex hull is the "z coordinate" of the position on a
            # (hyper)plane in the (number of independent concentrations +
            # 1)-dimensional (N-D) space. This plane is spanned by N points.
            # If there are more vertices on the convex hull, we need to loop
            # over all combinations of N vertices.
            for plane in itertools.combinations(
                vertices, min(len(vertices), self.dimensions + 1)
            ):
                # Calculate energy that would be gotten with pure components
                # with ascribed concentration.
                energy_pure = griddata(
                    self.concentrations[np.array(plane)],
                    self.energies[np.array(plane)],
                    concentration,
                    method="linear",
                )

                # Prepare to delete if the energy was lowered. `griddata` gives
                # NaN if the concentration is outside the triangle formed by
                # the three vertices. The result of the below comparison is
                # then False, which is what we want.
                if energy_pure < self.energies[i] - tol:
                    to_delete.append(i)
                    break

        # Finally remove all points
        self.concentrations = np.delete(self.concentrations, to_delete, 0)
        self.energies = np.delete(self.energies, to_delete, 0)
        self.structures = list(np.delete(self.structures, to_delete, 0))

    def get_energy_at_convex_hull(
        self, target_concentrations: Union[List[float], List[List[float]]]
    ) -> np.ndarray:
        """Returns the energy of the convex hull at specified concentrations.
        If any concentration is outside the allowed range, NaN is
        returned.

        Parameters
        ----------
        target_concentrations
            Concentrations at target points.

            If there is one independent concentration, a list of
            floats is sufficient. Otherwise, the concentrations ought
            to be provided as a list of lists, such as ``[[0.1, 0.2],
            [0.3, 0.1], ...]``.
        """
        if self.dimensions > 1 and isinstance(target_concentrations[0], Sized):
            assert len(target_concentrations[0]) == self.dimensions

        # Loop over all complexes of N+1 points to make sure that the lowest
        # energy plane is used in the end. This is needed in two dimensions
        # but in higher.
        hull_candidate_energies = []
        for plane in itertools.combinations(
            range(len(self.energies)), min(len(self.energies), self.dimensions + 1)
        ):
            try:
                plane_energies = griddata(
                    self.concentrations[list(plane)],
                    self.energies[list(plane)],
                    np.array(target_concentrations),
                    method="linear",
                )
            except QhullError:
                # If the points lie on a line, the convex hull will fail, but
                # we do not need to care about these "planes" anyway
                continue
            hull_candidate_energies.append(plane_energies)

        # Pick out the lowest energies found
        hull_energies = np.nanmin(hull_candidate_energies, axis=0)
        return hull_energies

    def extract_low_energy_structures(
        self,
        concentrations: Union[List[float], List[List[float]]],
        energies: List[float],
        energy_tolerance: float,
    ) -> List[int]:
        """Returns the indices of energies that lie within a certain
        tolerance of the convex hull.

        Parameters
        ----------
        concentrations
            Concentrations of candidate structures.

            If there is one independent concentration, a list of
            floats is sufficient. Otherwise, the concentrations must
            be provided as a list of lists, such as ``[[0.1, 0.2],
            [0.3, 0.1], ...]``.
        energies
            Energies of candidate structures.
        energy_tolerance
            Include structures with an energy that is at most this far
            from the convex hull.
        """
        # Convert to numpy arrays, can be necessary if, for example,
        # they are Pandas Series with "gaps"
        concentrations = np.array(concentrations)
        energies = np.array(energies)

        n_points = len(concentrations)
        if len(energies) != n_points:
            raise ValueError("concentrations and energies must have " "the same length")

        # Calculate energy at convex hull for specified concentrations
        hull_energies = self.get_energy_at_convex_hull(concentrations)

        # Extract those that are close enough
        close_to_hull = [
            i
            for i in range(n_points)
            if energies[i] <= hull_energies[i] + energy_tolerance
        ]

        return close_to_hull
