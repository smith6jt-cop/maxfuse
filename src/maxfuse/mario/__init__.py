"""
MARIO: Matching with Approximate Regional Importance with Overlap

Statistical matchability testing and bipartite matching for cross-modal integration.
"""

from .match import Mario, pipelined_mario
from . import match_utils, cluster, embed, utils

__all__ = ["Mario", "pipelined_mario", "match_utils", "cluster", "embed", "utils"]
