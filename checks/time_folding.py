from fractions import Fraction
from random import Random


def station(sites, edges, index, radius):
    period = len(sites)
    forward = [
        (sites[(index + offset) % period], edges[(index + offset) % period])
        for offset in range(-radius, radius)
    ]
    backward = [
        (sites[(index - offset) % period], -edges[(index - offset - 1) % period])
        for offset in range(-radius, radius)
    ]

    def directed(window):
        point_term = sum((offset + 2) * (site + 3 * edge) ** 2
                         for offset, (site, edge) in enumerate(window))
        link_term = window[0][0] * window[-1][1] ** 3
        return point_term + link_term

    return directed(forward) + directed(backward)


def fold(sites, edges, start, length, period):
    pause_start = period // 8
    pause_end = period - 2 * length - pause_start
    steps = [0] * pause_start + [1] * length + [0] * pause_end + [-1] * length
    read_index = start
    folded_sites = []
    folded_edges = []
    for step in steps:
        folded_sites.append(sites[read_index % len(sites)])
        match step:
            case 1:
                folded_edges.append(edges[read_index % len(sites)])
            case -1:
                folded_edges.append(-edges[(read_index - 1) % len(sites)])
            case 0:
                folded_edges.append(0)
        read_index += step
    assert read_index == start
    return folded_sites, folded_edges


def verify(period, long_period, radius, seed):
    random = Random(seed)
    sites = [random.randrange(-80, 80) for _ in range(long_period)]
    edges = [random.randrange(-80, 80) for _ in range(long_period)]
    long_average = Fraction(
        sum(station(sites, edges, index, radius) for index in range(long_period)),
        long_period,
    )
    frozen_average = Fraction(
        sum(station([site] * period, [0] * period, 0, radius) for site in sites),
        long_period,
    )
    folded_averages = []
    seam_totals = []
    lengths = [period // 8, period // 4]
    for length in lengths:
        sum_all = 0
        for start in range(long_period):
            folded_sites, folded_edges = fold(sites, edges, start, length, period)
            sum_all += sum(station(folded_sites, folded_edges, index, radius)
                           for index in range(period))
        average = Fraction(sum_all, long_period)
        folded_averages.append(average)
        seam_totals.append(average - 2 * (length - 2 * radius - 1) * long_average
                           - (period - 2 * length - 4 * radius - 2) * frozen_average)
    reconstructed = frozen_average + (folded_averages[1] - folded_averages[0]) / (
        2 * (lengths[1] - lengths[0])
    )
    assert seam_totals[0] == seam_totals[1]
    assert reconstructed == long_average
    print(f"n={period} q={long_period} R={radius} both residuals exactly zero")


if __name__ == "__main__":
    verify(129, 3 * 129, 1, 31)
    verify(193, 2 * 193, 2, 33)
    verify(257, 2 * 257 + 11, 3, 53)
