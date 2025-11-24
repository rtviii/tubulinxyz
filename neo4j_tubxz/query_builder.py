# neo4j_tubxz/query_builder.py
from typing import List, Dict, Any, Optional

class CypherQueryBuilder:
    """Helper to build complex Cypher queries with filters"""
    
    def __init__(self):
        self.parts: List[str] = []
        self.params: Dict[str, Any] = {}
        
    def match(self, pattern: str):
        self.parts.append(f"MATCH {pattern}")
        return self
    
    def where(self, *conditions: str):
        """Add WHERE clause with multiple conditions"""
        valid = [c for c in conditions if c]
        if valid:
            if "WHERE" not in " ".join(self.parts):
                self.parts.append("WHERE " + " AND ".join(valid))
            else:
                self.parts.append("AND " + " AND ".join(valid))
        return self
    
    def with_clause(self, *items: str):
        self.parts.append("WITH " + ", ".join(items))
        return self
    
    def order_by(self, *fields: str):
        self.parts.append("ORDER BY " + ", ".join(fields))
        return self
    
    def return_clause(self, *items: str):
        self.parts.append("RETURN " + ", ".join(items))
        return self
    
    def add_param(self, key: str, value: Any):
        self.params[key] = value
        return self
    
    def build(self) -> tuple[str, Dict[str, Any]]:
        return "\n".join(self.parts), self.params