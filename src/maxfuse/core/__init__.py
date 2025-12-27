"""
Core MaxFuse functionality for cross-modal integration.
"""

from .model import Fusor
from . import match_utils, graph, utils, spatial_utils, metrics

__all__ = ["Fusor", "match_utils", "graph", "utils", "spatial_utils", "metrics"]
