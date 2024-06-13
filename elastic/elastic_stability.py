"""
Born 弹性稳定性判据

reference: https://github.com/MineralsCloud/elastic-stability-conditions/blob/master/src/escpy/crystal.py
"""

import abc
from typing import List

import numpy as np


class CrystalSystem:
    def __init__(self, stiffness_matrix):
        self.stiffness_matrix = stiffness_matrix

    @property
    @abc.abstractmethod
    def ns_conditions_text(self) -> List[str]: ...

    @property
    @abc.abstractmethod
    def ns_conditions(self) -> List[bool]: ...

    def check_ns_conditions(self, outfile=None) -> None:
        flag = True

        print(
            "Born elastic stabilitycoditions are:\n",
            self.ns_conditions_text,
            "\n",
            file=outfile,
        )

        flag_list = []
        for criterion, criterion_text in zip(
            self.ns_conditions, self.ns_conditions_text
        ):
            if not criterion:  # If `criterion` evaluates to `False`.
                print(
                    "Criterion ${0}$ is not satisfied!".format(criterion_text),
                    file=outfile,
                )
                flag = False

            flag_list.append(flag)

        if all(flag_list):
            print("\nAll criteria are satisfied!", file=outfile)
        else:
            print("\nSome criteria are not satisfied!", file=outfile)

    @property
    def compliance_matrix(self):
        return np.linalg.inv(self.stiffness_matrix)


class CubicSystem(CrystalSystem):
    @property
    def ns_conditions_text(self):
        return ["C_{11} > | C_{12} |", "C_{11} + 2 C_{12} > 0", "C_{44} > 0"]

    @property
    def ns_conditions(self):
        c = self.stiffness_matrix
        c11, c12, c44 = c[0, 0], c[0, 1], c[3, 3]

        return [c11 > abs(c12), c11 + 2 * c12 > 0, c44 > 0]


class HexagonalSystem(CrystalSystem):
    @property
    def ns_conditions_text(self):
        return [
            "C_{11} > | C_{12} |",
            "2 C_{13}^2 < C_{33} (C_{11} + C_{12})",
            "C_{44} > 0",
            "C_{66} > 0",
        ]

    @property
    def ns_conditions(self):
        c = self.stiffness_matrix
        c11, c12, c13, c33, c44 = c[0, 0], c[0, 1], c[0, 2], c[2, 2], c[3, 3]
        c66 = (c11 - c12) / 2

        return [
            c11 > abs(c12),
            2 * c13**2 < c33 * (c11 + c12),
            c44 > 0,
            c66 > 0,
        ]


class TetragonalSystem(CrystalSystem):
    @property
    def ns_conditions_text(self):
        if self.stiffness_matrix[0, 5] == 0:  # Tetragonal (I) class
            return [
                "C_{11} > | C_{12} |",
                "2 C_{13}^2 < C_{33} (C_{11} + C_{12})",
                "C_{44} > 0",
                "C_{66} > 0",
            ]

        return [  # Tetragonal (II) class
            "C_{11} > | C_{12} |",
            "2 C_{13}^2 < C_{33} (C_{11} + C_{12})",
            "C_{44} > 0",
            "2 C_{16}^2 < C_{66} (C_{11} - C_{12})",
        ]

    @property
    def ns_conditions(self):
        c = self.stiffness_matrix
        c11, c12, c13, c16, c33, c44, c66 = (
            c[0, 0],
            c[0, 1],
            c[0, 2],
            c[0, 5],
            c[2, 2],
            c[3, 3],
            c[5, 5],
        )

        if c16 == 0:  # Tetragonal (I) class
            return [
                c11 > abs(c12),
                2 * c13**2 < c33 * (c11 + c12),
                c44 > 0,
                c66 > 0,
            ]

        return [  # Tetragonal (II) class
            c11 > abs(c12),
            2 * c13**2 < c33 * (c11 + c12),
            c44 > 0,
            2 * c16**2 < c66 * (c11 - c12),
        ]


class RhombohedralSystem(CrystalSystem):
    @property
    def ns_conditions_text(self):
        if self.stiffness_matrix[0, 4] == 0:  # Rhombohedral (I) class
            return [
                "C_{11} > | C_{12} |",
                "C_{44} > 0",
                "C_{13}^2 < 1/2 C_{33} (C_{11} + C_{12})",
                "C_{14}^2 < 1/2 C_{44} (C_{11} - C_{12})",
                "1/2 C_{44} (C_{11} - C_{12}) = C_{44} * C_{66}",
            ]

        return [  # Rhombohedral (II) class
            "C_{11} > | C_{12} |",
            "C_{44} > 0",
            "C_{13}^2 < 1/2 C_{33} (C_{11} + C_{12})",
            "C_{14}^2 + C_{15}^2 < 1/2 C_{44} (C_{11} - C_{12})",
            "1/2 C_{44} (C_{11} - C_{12}) = C_{44} * C_{66}",
        ]

    @property
    def ns_conditions(self):
        c = self.stiffness_matrix
        c11, c12, c13, c14, c15, c33, c44, c66 = (
            c[0, 0],
            c[0, 1],
            c[0, 2],
            c[0, 3],
            c[0, 4],
            c[2, 2],
            c[3, 3],
            c[5, 5],
        )

        if c15 == 0:  # Rhombohedral (I) class
            return [
                c11 > abs(c12),
                c44 > 0,
                c13**2 < 0.5 * c33 * (c11 + c12),
                c14**2 < 0.5 * c44 * (c11 - c12),
                0.5 * c44 * (c11 - c12) == c44 * c66,
            ]

        return [  # Rhombohedral (II) class
            c11 > abs(c12),
            c44 > 0,
            c13**2 < 0.5 * c33 * (c11 + c12),
            c14**2 + c15**2 < 0.5 * c44 * (c11 - c12),
            0.5 * c44 * (c11 - c12) == c44 * c66,
        ]


class OrthorhombicSystem(CrystalSystem):
    @property
    def ns_conditions_text(self):
        return [
            "C_{11} > 0",
            "C_{11} C_{22} > C_{12}^2",
            "C_{11} C_{22} C_{33} + 2 C_{12} C_{13} C_{23} > C_{11} C_{23}^2 + C_{22} C_{13}^2 + C_{33} C_{12}^2",
            "C_{44} > 0",
            "C_{55} > 0",
            "C_{66} > 0",
        ]

    @property
    def ns_conditions(self):
        c = self.stiffness_matrix
        c11, c22, c33, c44, c55, c66 = np.diag(c)
        c12, c13, c23 = c[0, 1], c[0, 2], c[1, 2]

        return [
            c11 > 0,
            c11 * c22 > c12**2,
            c11 * c22 * c33 + 2 * c12 * c13 * c23
            > c11 * c23**2 + c22 * c13**2 + c33 * c12**2,
            c44 > 0,
            c55 > 0,
            c66 > 0,
        ]


class MonoclinicSystem(CrystalSystem):
    @property
    def ns_conditions_text(self):
        return [
            "C_{ii} > 0",
            "C_{11} + C_{22} + C_{33} + 2 * (C_{12} + C_{13} + C_{23})",
            "C_{33} * C_{55} - C_{35}^2 > 0",
            "C_{44} * C_{66} - C_{46}^2 > 0",
            "C_{22} + C_{33} - 2 * C_{23}  > 0",
            "C_{22} * (C_{33} * C_{55} - C_{35}^2) + 2 * C_{23} * C_{25} * C_{35} - c23^2 * c55 - c25^2 * c33> 0",
            "The last criterion",
        ]

    @property
    def ns_conditions(self):
        c = self.stiffness_matrix
        c11, c22, c33, c44, c55, c66 = np.diag(c)
        c12, c13, c15, c23, c25, c35, c46 = (
            c[0, 1],
            c[0, 2],
            c[0, 4],
            c[1, 2],
            c[1, 4],
            c[2, 4],
            c[3, 5],
        )
        g = (
            c11 * c22 * c33
            - c11 * c23 * c23
            - c22 * c13 * c13
            - c33 * c12 * c12
            + 2 * c12 * c13 * c23
        )

        return [
            all(_ > 0 for _ in (c11, c22, c33, c44, c55, c66)),
            c11 + c22 + c33 + 2 * (c12 + c13 + c23) > 0,
            c33 * c55 - c35**2 > 0,
            c44 * c66 - c46**2 > 0,
            c22 + c33 - 2 * c23 > 0,
            c22 * (c33 * c55 - c35**2)
            + 2 * c23 * c25 * c35
            - c23**2 * c55
            - c25**2 * c33
            > 0,
            2
            * (
                c15 * c25 * (c33 * c12 - c13 * c23)
                + c15 * c35 * (c22 * c13 - c12 * c23)
                + c25 * c35 * (c11 * c23 - c12 * c13)
            )
            - (
                c15 * c15 * (c22 * c33 - c23 ^ 2)
                + c25 * c25 * (c11 * c33 - c13 ^ 2)
                + c35 * c35 * (c11 * c22 - c12 ^ 2)
            )
            + c55 * g
            > 0,
        ]


class TriclinicSystem(CrystalSystem):
    @property
    def ns_conditions_text(self):
        return None

    @property
    def ns_conditions(self):
        return None


if __name__ == "__main__":

    # 六方晶系
    hexagonal_matrix = np.array(
        [
            [-297.85, 126.9, 104.5, 0.0, 0.0, 0.0],
            [126.9, -297.85, 104.5, 0.0, 0.0, 0.0],
            [104.5, 104.5, 286.9, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 59.225, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 59.225, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 85.475],
        ]
    )

    hexagonal_system = HexagonalSystem(hexagonal_matrix)
    hexagonal_system.check_ns_conditions()

    # 立方晶系
    cubic_matrix = np.array(
        [
            [1050.640, 126.640, 126.640, 0.000, 0.000, 0.000],
            [126.640, 1050.640, 126.640, 0.000, 0.000, 0.000],
            [126.640, 126.640, 1050.640, 0.000, 0.000, 0.000],
            [0.000, 0.000, 0.000, 559.861, 0.000, 0.000],
            [0.000, 0.000, 0.000, 0.000, 559.861, 0.000],
            [0.000, 0.000, 0.000, 0.000, 0.000, 559.861],
        ]
    )

    cubic_system = CubicSystem(cubic_matrix)
    cubic_system.check_ns_conditions()
