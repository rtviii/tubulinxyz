from typing import List
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from api.musle_alignment import Annotation, TubulinAlignmentMapper

class AlignmentRequest(BaseModel):
    sequence: str
    sequence_id: str
    annotations: List[dict] = []

class AlignmentResponse(BaseModel):
    aligned_sequence: str
    mapped_annotations: List[dict]
    statistics: dict
    mapping: List[int]

app    = FastAPI()
mapper = TubulinAlignmentMapper("data/master_alpha.msa", "muscle3.8.1")

@app.post("/api/align-sequence", response_model=AlignmentResponse)
async def align_sequence(request: AlignmentRequest):
    try:
        aligned_seq, mapping = mapper.align_sequence(
            request.sequence_id, 
            request.sequence
        )
        annotations = [
            Annotation(
                start=ann['start'],
                end=ann['end'],
                label=ann['label'],
                metadata=ann.get('metadata', {})
            ) for ann in request.annotations
        ]
        
        mapped_annotations = mapper.map_annotations(annotations, mapping)
        stats              = mapper.get_alignment_statistics(aligned_seq, mapping)
        
        return AlignmentResponse(
            aligned_sequence=aligned_seq,
            mapped_annotations=[{
                'start': ann.start,
                'end': ann.end, 
                'label': ann.label,
                'metadata': ann.metadata
            } for ann in mapped_annotations],
            statistics=stats,
            mapping=mapping
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))