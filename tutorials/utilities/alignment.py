import numpy as np
import numpy.typing as npt


def fit_template_points_to_observed(
    template_points: npt.NDArray,
    observed_points: npt.NDArray,
):
    # https://web.stanford.edu/class/cs273/refs/umeyama.pdf
    # https://doi.org/10.1109/34.88573
    p = template_points.T
    q = observed_points.T

    # centroids
    cen_P = np.mean(p, axis=1).reshape(-1, 1)
    cen_Q = np.mean(q, axis=1).reshape(-1, 1)

    # centered vectors
    X = p - cen_P
    Y = q - cen_Q

    # svd
    U, sigma, Vt = np.linalg.svd(X @ Y.T)

    # reflection correction
    d = np.identity(Vt.T.shape[1])
    d[-1, -1] = np.linalg.det(Vt.T @ U.T)

    # rotation and translation
    R = Vt.T @ d @ U.T
    t = cen_Q - R @ cen_P

    # convert to TRS matrix4x4
    translation = np.identity(4)
    translation[:3, 3] = t.reshape(-1)
    rotation = np.identity(4)
    rotation[:3, :3] = R

    return translation @ rotation
