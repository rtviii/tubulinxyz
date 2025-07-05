#!/usr/bin/env python3
"""
Start the FastAPI server for tubulin analyzer
Run this from the project root directory
"""

import uvicorn
from pathlib import Path

if __name__ == "__main__":
    # Ensure debug output directory exists
    Path("debug_output").mkdir(exist_ok=True)
    
    # Start the FastAPI server
    uvicorn.run(
        "api.main:app", 
        host="127.0.0.1", 
        port=8000, 
        reload=True,
        reload_dirs=["tubulin_analyzer", "api"]
    )