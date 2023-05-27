import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
y = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
lam = np.array(
    [
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
        [
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
            705,
            720,
            740,
            760,
            800,
            780,
            800,
            820,
            840,
            900,
        ],
    ]
)

lam_shift = np.array(
    [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
)


def gen_shift(x, y, lam, a, b, c):
    subpixel_length = (x[-1] - x[0]) / (len(x) - 1)
    for i in range(len(x)):
        X = x[i]
        shift = int((a * X**2 + b * X + c) / subpixel_length)
        for j in range(len(y)):
            if j + shift in y:
                lam[i][j] = lam[i][j + shift]
                lam_shift[i][j] = shift
    shifts = [lam_shift[i][0] for i in range(len(lam_shift))]
    plt.plot(x, shifts)
    plt.show()


gen_shift(x, y, lam, -1 / 14, 19 / 14, 0)
