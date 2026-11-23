"""Check the finite-frequency coarea identity in a charged Cartan branch."""

import numpy as np


def null_transverse(gradient):
    _, _, right = np.linalg.svd(gradient.conj()[None, :])
    return right.conj().T[:, 1:]


def curl_matrix(spatial, temporal):
    curl = np.zeros((6, 3), dtype=complex)
    curl[:3, :2] = temporal * spatial[0]
    curl[:3, 2] = -spatial[1]
    row = 3
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        curl[row, :2] = spatial[1][i] * spatial[0][j] - spatial[1][j] * spatial[0][i]
        row += 1
    return curl


def relative(left, right):
    return abs(left - right) / max(abs(left), abs(right))


def verify(n, k, angle, temporal):
    p0 = np.expm1(2j * np.pi * k / n)
    charge = np.expm1(1j * (2 * np.pi * k + angle) / n)
    carry = np.exp(0.5j * angle / n) * np.sinc(angle / (2 * np.pi * n))
    ordinary = carry[:, None] * null_transverse(p0)
    charged = null_transverse(charge)
    transverse = charged.conj().T @ ordinary
    combined = np.column_stack([ordinary, charge])
    coarea = abs(np.linalg.det(combined))
    coarea_factored = abs(np.linalg.det(transverse)) * np.linalg.norm(charge)

    curl = curl_matrix((ordinary, charge), np.expm1(1j * temporal))
    gaussian = np.sqrt(np.linalg.det(curl.conj().T @ curl).real)
    gaussian_factored = coarea_factored * (np.linalg.norm(charge) ** 2 + abs(np.expm1(1j * temporal)) ** 2)

    basis, _ = np.linalg.qr(curl)
    projector = basis @ basis.conj().T
    curve = np.eye(6) + 0.25 * np.sin(temporal) * np.diag(np.arange(1, 7) / 6)
    direct = np.linalg.det(curl.conj().T @ curve @ curl) / gaussian**2
    whitened = np.linalg.det(np.eye(6) + projector @ (curve - np.eye(6)) @ projector)
    return [relative(coarea, coarea_factored), relative(gaussian, gaussian_factored), relative(direct, whitened)]


def main():
    for n in (7, 23, 89, 257):
        errors = np.zeros(3)
        for k in [np.array(x) for x in ((1, 0, 0), (0, -1, 2), (-1, 1, 0))]:
            for angle in [np.array(x) for x in ((0, 0, 0), (3.7, -2.3, 0.63), (-3.7, 3.5, 3.2))]:
                for temporal in (0.0, 0.03, np.pi / 2, np.pi):
                    errors = np.maximum(errors, verify(n, k, angle, temporal))
        print(f"n={n:3d} coarea={errors[0]:.3e} gaussian={errors[1]:.3e} projection={errors[2]:.3e}")
        assert max(errors) < 3e-10


if __name__ == "__main__":
    main()
