"""
MARIO: Matching with Approximate Regional Importance with Overlap

Statistical matchability testing and bipartite matching for cross-modal integration.
"""

from .match import Mario
from . import match_utils, cluster, embed, utils

__all__ = ["Mario", "match_utils", "cluster", "embed", "utils"]
