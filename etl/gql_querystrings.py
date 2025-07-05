EntryInfoString = """
{
  entry(entry_id: "$RCSB_ID") {
  rcsb_id
  rcsb_accession_info{
    deposit_date
  }
    struct_keywords {
      pdbx_keywords
      text
    }
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_external_references {
      link
      type
      id
    }
    exptl {
      method
    }
    citation {
      rcsb_authors
      year
      title
      pdbx_database_id_DOI
    }
    struct_keywords {
      pdbx_keywords
      text
    }
  }
}
"""

AssemblyIdentificationString = """
{
  entry(entry_id: "$RCSB_ID") {
    assemblies {
      rcsb_id
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers {
          comp_id
          auth_asym_id
          rcsb_id
          auth_seq_id
          entity_id
        }
      }
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          asym_id
          auth_asym_id
          entry_id
          entity_id
        }
      }
    }
  }
}
"""

PolymerEntitiesString = """{
  entry(entry_id: "$RCSB_ID") {
    
    rcsb_id
    assemblies {
      polymer_entity_instances {
        rcsb_id
      }
    }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
        entry_id
        entity_id
      }
      pfams {
        rcsb_pfam_accession
        rcsb_pfam_comment
        rcsb_pfam_description
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        scientific_name
      }
      rcsb_entity_host_organism {
        ncbi_taxonomy_id
        scientific_name
      }
      uniprots {
        rcsb_id
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
      rcsb_polymer_entity_annotation {
        annotation_id
        assignment_version
        description
        name
        provenance_source
        type
      }
    }
  }
}
"""

NonpolymerEntitiesString = """{
  entry(entry_id: "$RCSB_ID") {
    nonpolymer_entities {
      pdbx_entity_nonpoly {
        comp_id
        name
        entity_id
      }
      rcsb_nonpolymer_entity {
        details
        formula_weight
        pdbx_description
        pdbx_number_of_molecules
      }
      nonpolymer_comp {
        chem_comp {
          id
          name
          three_letter_code
        }
        rcsb_chem_comp_target {
          interaction_type
          name
          provenance_source
          reference_database_accession_code
          reference_database_name
        }
        drugbank {
          drugbank_container_identifiers {
            drugbank_id
          }
          drugbank_info {
            cas_number
            description
            indication
            mechanism_of_action
            name
            pharmacology
          }
        }
      }
    }
  }
}
"""

LigandsChemInfo = """{
  chem_comps(comp_ids: $COMP_IDS) {
    chem_comp{
      id
    }
    rcsb_chem_comp_descriptor {
      SMILES
      InChI
      InChIKey
      SMILES_stereo
    }
  }
}
"""

# https://data.rcsb.org/graphql/index.html
#///////////////////////////////////////////////////////////////////////////////////////////////

#!DEPRECATED - Use each subcomponent query separately.
WholeStructureTemplate = """{
  entry(entry_id: "$RCSB_ID") {
    assemblies{
        rcsb_id 
       	nonpolymer_entity_instances{
      
          rcsb_nonpolymer_entity_instance_container_identifiers{
            auth_asym_id
            auth_seq_id
            entity_id
          }
        }
        polymer_entity_instances{
           rcsb_polymer_entity_instance_container_identifiers {        
            auth_asym_id
            entity_id
          }
        }
      }


      



    rcsb_id
      rcsb_accession_info{
        deposit_date
      }
      struct_keywords {
        pdbx_keywords
        text
      }
      rcsb_entry_info {
        resolution_combined
      }
      rcsb_external_references {
        link
        type
        id
      }
      exptl {
        method
      }
      citation {
        rcsb_authors
        year
        title
        pdbx_database_id_DOI
      }
      struct_keywords {
        pdbx_keywords
        text
      }




      


    polymer_entities {
      entry{
        rcsb_id
      }
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
        entry_id
        entity_id
      }
      pfams {
        rcsb_pfam_accession
        rcsb_pfam_comment
        rcsb_pfam_description
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        scientific_name
      }

      rcsb_entity_host_organism {
        ncbi_taxonomy_id
        scientific_name
      }

      uniprots {
        rcsb_id
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
    }


    nonpolymer_entities {
      pdbx_entity_nonpoly {
        comp_id
        name
        entity_id
      }
      rcsb_nonpolymer_entity {
        details
        formula_weight
        pdbx_description
        pdbx_number_of_molecules
      }
      nonpolymer_comp {
        chem_comp {
          id
          name
          three_letter_code
        }
        rcsb_chem_comp_target {
          interaction_type
          name
          provenance_source
          reference_database_accession_code
          reference_database_name
        }

        drugbank {
          drugbank_container_identifiers {
            drugbank_id
          }
          drugbank_info {
            cas_number
            description
            indication
            mechanism_of_action
            name
            pharmacology
          }
        }
      }
      
    }


  }
}"""