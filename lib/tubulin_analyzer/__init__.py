"""
Tubulin Analyzer - A library for analyzing tubulin protein structures

This library provides tools for:
- Parsing PDB/mmCIF structures
- Analyzing protein-protein interfaces
- Geometric analysis and alignment
- Grid generation from 3D structures
- Visualization utilities
"""

from lib.tubulin_analyzer.structural_analyzer import SpatialGridGenerator
from lib.tubulin_analyzer.interface_analyzer import InterfaceAnalyzer
from lib.tubulin_analyzer.geometric_analyzer import GeometricAnalyzer
from lib.tubulin_analyzer.visualization import VisualizationUtils
from lib.tubulin_analyzer.models import SubunitData, GridData, DebugData

__version__ = "0.1.0"
__author__ = "Your Name"

__all__ = [
    "SpatialGridGenerator",
    "InterfaceAnalyzer", 
    "GeometricAnalyzer",
    "VisualizationUtils",
    "SubunitData",
    "GridData",
    "DebugData",
]