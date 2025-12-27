"""
MaxFuse: Cross-modal single-cell data integration with spatial awareness

This package provides tools for integrating single-cell datasets from different
modalities with no overlapping features and/or under low signal-to-noise ratio regimes.
"""

from . import core
from . import mario
from .core.model import Fusor
from .mario.match import Mario

__version__ = "0.1.0"
__all__ = ["Fusor", "Mario", "core", "mario"]
