try:
    import pytest

    pytest.register_assert_rewrite("nanover.testing.asserts")
except ImportError:
    pass

from .asserts import (
    assert_equal_soon as assert_equal_soon,
    assert_in_soon as assert_in_soon,
    assert_not_in_soon as assert_not_in_soon,
)
