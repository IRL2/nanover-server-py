import numpy as np
import pytest
from hypothesis import given
from nanover.jupyter.utilities import TransformsUtility
from nanover.testing.servers import make_app_server
from nanover.testing.strategies import transformation
from nanover.utilities.transforms import Transform


@pytest.fixture(scope="module")
def reusable_transforms_utility():
    with make_app_server() as app_server:
        yield TransformsUtility.from_client(app_server)


@given(a=transformation(), b=transformation(), c=transformation())
def test_fetch_transform_root(reusable_transforms_utility, a, b, c):
    transforms = reusable_transforms_utility
    transforms.update_transform("a", transform=Transform.from_local_to_parent_matrix(a))
    transforms.update_transform(
        "b", transform=Transform.from_local_to_parent_matrix(b), parent="a"
    )
    transforms.update_transform(
        "c", transform=Transform.from_local_to_parent_matrix(c), parent="b"
    )
    abc = transforms.fetch_transform_root("c")

    assert np.allclose(abc.local_to_parent_matrix, a @ b @ c, atol=1e-7)
