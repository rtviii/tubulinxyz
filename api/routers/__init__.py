# api/routers/__init__.py
from api.routers.router_structures import router_structures
from api.routers.router_polymers import router_polymers
from api.routers.router_ligands import router_ligands
from api.routers.router_msa import router_msa
from api.routers.router_grid import router_grid

__all__ = [
    "router_structures",
    "router_polymers",
    "router_ligands",
    "router_msa",
    "router_grid",
]