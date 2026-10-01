"""Diffusion on a grid of voxels in one, two and three dimensions.

The operator is checked against exact discrete eigenvectors, the 3-D point
release, and its own conservation. It is checked for invariance under the
symmetries of the box, and for agreement between dimensions. The multigrid
solver is checked for grid-independent convergence and against a direct
solve. docs/validation.md, "Transport in one, two and three dimensions".
"""

import math

import numpy as np
import pytest

from marse.core.framestore import FrameStore as CoreFrameStore
from marse.ecosystem.framestore import FrameStore as EcosystemFrameStore
from marse.spatial.column import ColumnSystem
from marse.spatial.grid import Grid
from marse.spatial.multigrid import ImplicitSystem, inverse
from marse.spatial.transport import Diffusion
from marse.spatial.vtk import read_vti, write_pvd, write_vti
from marse.validation.analytical import cosine_mode_rate, point_source_diffusion_3d

SHAPES = {1: (16,), 2: (12, 16), 3: (12, 10, 16)}


def cosine_mode(grid: Grid, modes: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    """A product of cosines that the boundary conditions admit, and its wavenumbers."""
    factors, wavenumbers = [], []
    for axis in range(grid.dimensions):
        length = grid.size_um[axis]
        if axis == grid.height_axis:  # no flux at the substratum, zero at the top face
            k = (2 * modes[axis] + 1) * np.pi / (2 * length)
        else:  # periodic
            k = 2 * np.pi * modes[axis] / length
        wavenumbers.append(k)
        shape = [1] * grid.dimensions
        shape[axis] = -1
        factors.append(np.cos(k * grid.centers_um(axis)).reshape(shape))
    return np.prod(np.broadcast_arrays(*factors), axis=0), np.array(wavenumbers)


@pytest.mark.parametrize("dims", [1, 2, 3])
def test_a_cosine_mode_is_an_exact_eigenvector(dims):
    grid = Grid(SHAPES[dims], 2.0)
    diffusivity = np.array([3.0e6, 1.0e6, 0.0])
    mode, k = cosine_mode(grid, (2, 1, 3)[-dims:])
    c = np.stack([mode, 2 * mode, mode])
    rate = Diffusion(grid, diffusivity, np.zeros(3)).rate(c)
    for j, d in enumerate(diffusivity):
        expected = cosine_mode_rate(k, d, grid.voxel_um) * c[j]
        np.testing.assert_allclose(rate[j], expected, rtol=0, atol=1e-12 * 3.0e6)


@pytest.mark.numerical
def test_the_discrete_rate_converges_to_the_continuum_at_second_order():
    k = np.array([2 * np.pi / 40.0, 2 * np.pi / 40.0, 3 * np.pi / 80.0])
    exact = cosine_mode_rate(k, 1.0)
    errors = [abs(cosine_mode_rate(k, 1.0, h) - exact) for h in (1.0, 0.5, 0.25, 0.125)]
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(orders > 1.99), orders


@pytest.mark.numerical
def test_the_3d_operator_is_second_order_on_a_point_release():
    # Truncation error: the operator applied to the exact spreading release,
    # against its exact time derivative. The release is 10 um wide in a 128 um
    # box, so the boundaries are 1e-9 of the peak away.
    d, t, amount = 1.0e6, 5.0e-5, 1.0
    errors = []
    for n in (32, 64, 128):
        h = 128.0 / n
        grid = Grid((n, n, n), h)
        centre = 64.0
        axes = [grid.centers_um(a) - centre for a in range(3)]
        r2 = axes[0][:, None, None] ** 2 + axes[1][None, :, None] ** 2 + axes[2][None, None, :] ** 2
        c = point_source_diffusion_3d(np.sqrt(r2), t, amount, d)
        exact = c * (r2 / (4 * d * t**2) - 1.5 / t)
        rate = Diffusion(grid, np.array([d]), np.zeros(1)).rate(c[None])[0]
        errors.append(np.max(np.abs(rate - exact)) / np.max(np.abs(exact)))
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(orders > 1.9), orders


@pytest.mark.parametrize("dims", [1, 2, 3])
def test_diffusion_moves_matter_only_through_the_top_face(dims):
    grid = Grid(SHAPES[dims], 1.5)
    rng = np.random.default_rng(dims)
    c = rng.uniform(0.0, 3.0, size=(4, *grid.shape))
    diffusion = Diffusion(grid, np.array([2.0e6, 5.0e5, 1.0, 0.0]), np.array([1.0, 0.0, 2.0, 5.0]))
    inside = diffusion.amounts(diffusion.rate(c))
    through_top = diffusion.import_rate(c)
    np.testing.assert_allclose(
        inside, through_top, rtol=1e-13, atol=1e-13 * np.abs(through_top).max()
    )
    assert inside[3] == 0.0  # a component without diffusivity does not move


def test_a_laterally_uniform_field_diffuses_exactly_as_one_column():
    rng = np.random.default_rng(3)
    column = rng.uniform(0.0, 1.0, size=(2, 16))
    diffusivity, bulk = np.array([1.0e6, 3.0e5]), np.array([0.2, 1.0])
    one = Diffusion(Grid((16,), 2.0), diffusivity, bulk).rate(column)
    box = np.broadcast_to(column[:, None, None, :], (2, 5, 4, 16)).copy()
    three = Diffusion(Grid((5, 4, 16), 2.0), diffusivity, bulk).rate(box)
    for i, j in np.ndindex(5, 4):
        np.testing.assert_array_equal(three[:, i, j, :], one)


def test_swapping_or_reflecting_lateral_axes_moves_the_rate_with_them():
    rng = np.random.default_rng(4)
    c = rng.uniform(0.0, 1.0, size=(2, 6, 6, 8))
    diffusion = Diffusion(Grid((6, 6, 8), 2.0), np.array([1.0e6, 2.0e5]), np.array([0.5, 0.0]))
    rate = diffusion.rate(c)
    swapped = diffusion.rate(np.swapaxes(c, 1, 2))
    np.testing.assert_allclose(swapped, np.swapaxes(rate, 1, 2), rtol=1e-14, atol=1e-9)
    reflected = diffusion.rate(c[:, ::-1])
    np.testing.assert_allclose(reflected, rate[:, ::-1], rtol=1e-14, atol=1e-9)


def test_the_explicit_limit_that_implicit_steps_avoid():
    # Oxygen at 37 C (2624 um^2/s) on 2 um voxels in 3-D: a quarter of a millisecond.
    diffusion = Diffusion(Grid((4, 4, 4), 2.0), np.array([2624.0 * 3600]), np.zeros(1))
    assert diffusion.explicit_step_limit_h() * 3600 == pytest.approx(2.54e-4, rel=1e-2)


def random_system(shape, components=4, a=0.03, seed=1):
    rng = np.random.default_rng(seed)
    diffusivity = np.array([4.0e6, 1.0e6, 3.0e5, 0.0][:components])
    blocks = rng.normal(0.0, 50.0, size=(components, components, *shape))
    blocks[np.arange(components), np.arange(components)] = (
        -np.abs(blocks[np.arange(components), np.arange(components)]) * 20
    )
    return ImplicitSystem(shape, (2.0,) * len(shape), diffusivity, a, blocks), rng


@pytest.mark.parametrize("shape", [(8, 8, 8), (16, 16, 16), (32, 32, 16), (256,), (32, 32)])
def test_a_multigrid_cycle_reduces_the_error_by_a_factor_independent_of_the_grid(shape):
    system, rng = random_system(shape)
    assert len(system.levels) >= 2  # the cycle, not only the exact coarsest solve, is tested
    b = rng.normal(0.0, 1.0, size=(4, *shape))
    x = np.zeros_like(b)
    for _ in range(4):
        residual = b - system.apply(x)
        if np.linalg.norm(residual) < 1e-11 * np.linalg.norm(b):
            break  # converged to rounding; further ratios would compare noise
        x = x + system.precondition(residual)
        factor = np.linalg.norm(b - system.apply(x)) / np.linalg.norm(residual)
        assert factor < 0.15, (shape, factor)


def test_the_linear_solver_agrees_with_a_direct_solve():
    shape = (4, 6, 8)
    system, rng = random_system(shape, components=3)
    n = 3 * 4 * 6 * 8
    matrix = np.column_stack(
        [system.apply(np.eye(n)[k].reshape(3, *shape)).ravel() for k in range(n)]
    )
    b = rng.normal(0.0, 1.0, size=(3, *shape))
    direct = np.linalg.solve(matrix, b.ravel()).reshape(b.shape)
    solved, iterations = system.solve(b, scale=np.ones_like(b), tolerance=1e-12)
    np.testing.assert_allclose(solved, direct, rtol=1e-9, atol=1e-11)
    assert iterations < 20


def test_the_linear_solver_is_deterministic():
    system, rng = random_system((16, 16, 8))
    b = rng.normal(0.0, 1.0, size=(4, 16, 16, 8))
    first, _ = system.solve(b, scale=np.ones_like(b), tolerance=1e-10)
    second, _ = system.solve(b, scale=np.ones_like(b), tolerance=1e-10)
    np.testing.assert_array_equal(first, second)


@pytest.mark.numerical
def test_the_coarsest_level_is_inverted_as_lapack_would():
    rng = np.random.default_rng(7)
    for n in (1, 5, 33, 64, 200):  # one panel, several, and a partial last one
        matrix = rng.normal(size=(n, n)) + 0.5 * n * np.diag(rng.choice([-1.0, 1.0], n))
        expected = np.linalg.inv(matrix)
        np.testing.assert_allclose(
            inverse(matrix), expected, rtol=1e-12, atol=1e-12 * np.abs(expected).max()
        )
    swapped = np.array([[0.0, 2.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 3.0]])  # needs a row swap
    np.testing.assert_allclose(inverse(swapped) @ swapped, np.eye(3), atol=1e-15)
    with pytest.raises(np.linalg.LinAlgError):
        inverse(np.array([[1.0, 2.0], [2.0, 4.0]]))


def test_the_coarsest_matrix_is_the_operator_applied_to_every_unknown():
    shape, components = (4, 2, 3), 3
    system, _ = random_system(shape, components=components)
    level = system.levels[-1]
    assert level.dense_inverse is not None  # 72 unknowns: one level, solved exactly
    n = components * math.prod(shape)
    one_by_one = np.column_stack(
        [system.apply(np.eye(n)[k].reshape(components, *shape)).ravel() for k in range(n)]
    )
    np.testing.assert_array_equal(system._dense(level), one_by_one)
    np.testing.assert_allclose(level.dense_inverse @ one_by_one, np.eye(n), atol=1e-12)


@pytest.mark.parametrize("dims", [1, 2, 3])
def test_vtk_frames_read_back(tmp_path, dims):
    grid = Grid(SHAPES[dims], 2.5)
    rng = np.random.default_rng(dims)
    fields = {
        "oxygen_mol_per_m3": rng.uniform(0, 0.2, grid.shape),
        "lactate": rng.uniform(0, 9, grid.shape),
    }
    write_vti(tmp_path / "exact.vti", grid, fields, dtype="<f8")
    read, counts, spacing = read_vti(tmp_path / "exact.vti")
    assert spacing == 2.5
    assert counts[-1] == grid.shape[-1]
    for name, values in fields.items():
        np.testing.assert_array_equal(read[name].reshape(grid.shape), values)
    write_vti(tmp_path / "small.vti", grid, fields)
    small, _, _ = read_vti(tmp_path / "small.vti")
    np.testing.assert_allclose(small["lactate"].reshape(grid.shape), fields["lactate"], rtol=1e-7)


def test_a_pvd_file_lists_the_frames_with_their_times(tmp_path):
    path = write_pvd(tmp_path / "run.pvd", [(0.0, "frame_0000.vti"), (0.5, "frame_0001.vti")])
    text = path.read_text("ascii")
    assert text.index('timestep="0.0"') < text.index('timestep="0.5"')
    assert 'file="frame_0001.vti"' in text


def test_both_engines_record_frames_through_one_store():
    assert EcosystemFrameStore is CoreFrameStore


@pytest.mark.parametrize(
    ("shape", "voxel", "message"),
    [
        ((8, 8, 8, 8), 1.0, "one, two or three axes"),
        ((8, 1), 1.0, "at least two voxels"),
        ((8,), 0.0, "positive"),
    ],
)
def test_impossible_grids_are_refused(shape, voxel, message):
    with pytest.raises(ValueError, match=message):
        Grid(shape, voxel)


def test_the_footprint_counts_a_missing_lateral_axis_as_one_voxel_deep():
    assert Grid((10,), 2.0).footprint_um2 == 4.0
    assert Grid((5, 10), 2.0).footprint_um2 == 20.0
    assert Grid((5, 3, 10), 2.0).footprint_um2 == 60.0


# --- a closed top, and the direct solve of a column ---------------------------------


@pytest.mark.parametrize("dims", [1, 2, 3])
def test_a_closed_top_keeps_everything_in_the_box(dims):
    grid = Grid(SHAPES[dims], 1.5)
    rng = np.random.default_rng(dims)
    c = rng.uniform(0.0, 3.0, size=(3, *grid.shape))
    diffusion = Diffusion(grid, np.array([2.0e6, 5.0e5, 0.0]), np.ones(3), closed_top=True)
    assert np.all(diffusion.import_rate(c) == 0.0)
    rate = diffusion.rate(c)
    size = diffusion.amounts(np.abs(rate)).max()  # what the sum could lose to rounding
    np.testing.assert_allclose(diffusion.amounts(rate), 0.0, atol=1e-13 * size)
    # The diagonal is the rate's own derivative, closed face and all.
    unit = np.zeros((3, *grid.shape))
    corner = (0,) * (dims - 1) + (-1,)
    unit[(0, *corner)] = 1.0
    assert (
        diffusion.rate(unit, homogeneous=True)[(0, *corner)] == diffusion.diagonal()[(0, *corner)]
    )


def _column_blocks(n, components, seed):
    rng = np.random.default_rng(seed)
    blocks = rng.normal(0.0, 50.0, size=(components, components, n))
    diagonal = np.arange(components)
    blocks[diagonal, diagonal] = -np.abs(blocks[diagonal, diagonal]) * 20
    return blocks, rng


@pytest.mark.parametrize("closed_top", [False, True])
def test_a_column_is_solved_exactly_and_agrees_with_multigrid(closed_top):
    n, components, a = 64, 4, 0.03
    diffusivity = np.array([4.0e6, 1.0e6, 3.0e5, 0.0])
    blocks, rng = _column_blocks(n, components, 2)
    column = ColumnSystem((n,), (2.0,), diffusivity, a, blocks, closed_top=closed_top)
    grid = ImplicitSystem((n,), (2.0,), diffusivity, a, blocks, closed_top=closed_top)
    x = rng.normal(size=(components, n))
    np.testing.assert_array_equal(column.apply(x), grid.apply(x))
    b = rng.normal(size=(components, n))
    solved, iterations = column.solve(b, scale=np.ones_like(b))
    assert iterations == 1
    np.testing.assert_allclose(column.apply(solved), b, rtol=0, atol=1e-12)
    by_multigrid, _ = grid.solve(b, scale=np.ones_like(b), tolerance=1e-13)
    np.testing.assert_allclose(solved, by_multigrid, rtol=1e-9, atol=1e-11)


def test_a_column_is_solved_the_same_way_every_time():
    blocks, rng = _column_blocks(32, 3, 5)
    column = ColumnSystem((32,), (2.0,), np.array([1.0e6, 2.0e5, 0.0]), 0.1, blocks)
    b = rng.normal(size=(3, 32))
    np.testing.assert_array_equal(column.precondition(b), column.precondition(b))
