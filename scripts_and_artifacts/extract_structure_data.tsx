// scripts_and_artifacts/extract_structure_data.tsx

import { CIF } from 'molstar/lib/mol-io/reader/cif';
import { trajectoryFromMmCIF } from 'molstar/lib/mol-model-formats/structure/mmcif';
import { Structure, StructureElement, StructureProperties, StructureSelection, QueryContext, Unit } from 'molstar/lib/mol-model/structure';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { compile } from 'molstar/lib/mol-script/runtime/query/compiler';

import * as fs from 'fs/promises';
import * as path from 'path';

// ============================================================
// Output Types
// ============================================================

interface ObservedResidue {
    auth_seq_id : number;
    label_seq_id: number;
    comp_id     : string;
    one_letter  : string;
}

interface ObservedSequenceData {
    auth_asym_id: string;
    entity_id   : string;
    residues    : ObservedResidue[];
}

interface NeighborhoodResidue {
    auth_asym_id  : string;
    observed_index: number;
    comp_id       : string;
}

interface SimplifiedLigandNeighborhood {
    ligand_comp_id       : string;
    ligand_auth_asym_id  : string;
    ligand_auth_seq_id   : number;
    neighborhood_residues: NeighborhoodResidue[];
}

interface ExtractionResult {
    rcsb_id             : string;
    sequences           : ObservedSequenceData[];
    ligand_neighborhoods: SimplifiedLigandNeighborhood[];
}

// ============================================================
// Residue Conversion
// ============================================================

const STANDARD_AA: Record<string, string> = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
};

const MODIFIED_AA: Record<string, string> = {
    'MSE': 'M', 'CSX': 'C', 'CSO': 'C', 'SEP': 'S', 'TPO': 'T',
    'PTR': 'Y', 'KCX': 'K', 'LLP': 'K', 'PCA': 'E', 'HIC': 'H',
    'MLY': 'K', 'M3L': 'K', 'ALY': 'K', 'CGU': 'E', 'CME': 'C',
};

function toOneLetter(compId: string): string | null {
    const upper = compId.toUpperCase();
    return STANDARD_AA[upper] ?? MODIFIED_AA[upper] ?? null;
}

function isProteinResidue(compId: string): boolean {
    return toOneLetter(compId) !== null;
}

// ============================================================
// Skip Lists
// ============================================================

const SKIP_LIGANDS = new Set([
    'HOH', 'DOD', 'WAT', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE',
    'MN', 'CO', 'NI', 'CU', 'SO4', 'PO4', 'NO3', 'CO3', 'GOL', 'EDO',
    'PEG', 'PGE', 'ACT', 'FMT', 'MES', 'TRS', 'HEP', 'EPE', 'CIT',
    'TAR', 'DMS', 'DMF', 'BME', 'DTT'
]);

// ============================================================
// Structure Loading
// ============================================================

async function loadStructureFromCif(cifContent: string): Promise<Structure> {
    const parsed = await CIF.parse(cifContent).run();

    if (parsed.isError) {
        throw new Error(`CIF parsing failed: ${parsed.message}`);
    }

    const block      = parsed.result.blocks[0];
    const trajectory = await trajectoryFromMmCIF(block).run();
    return Structure.ofModel(trajectory.representative);
}

// ============================================================
// Sequence Extraction
// ============================================================

function extractObservedSequences(structure: Structure): ObservedSequenceData[] {
    const chainDataMap = new Map<string, ObservedSequenceData>();

    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;

        const loc = StructureElement.Location.create(structure, unit, unit.elements[0]);
        const seenResidues = new Set<string>();

        for (let i = 0; i < unit.elements.length; i++) {
            loc.element = unit.elements[i];

            const authAsymId = StructureProperties.chain.auth_asym_id(loc);
            const entityId   = StructureProperties.entity.id(loc);
            const authSeqId  = StructureProperties.residue.auth_seq_id(loc);
            const labelSeqId = StructureProperties.residue.label_seq_id(loc);
            const compId     = StructureProperties.atom.auth_comp_id(loc);
            const resKey     = `${authAsymId}:${authSeqId}`;

            if (seenResidues.has(resKey)) continue;

            const oneLetter = toOneLetter(compId);
            if (!oneLetter) continue;

            seenResidues.add(resKey);

            if (!chainDataMap.has(authAsymId)) {
                chainDataMap.set(authAsymId, {
                    auth_asym_id: authAsymId,
                    entity_id: String(entityId),
                    residues: []
                });
            }

            chainDataMap.get(authAsymId)!.residues.push({
                auth_seq_id: authSeqId,
                label_seq_id: labelSeqId,
                comp_id: compId,
                one_letter: oneLetter
            });
        }
    }

    for (const chainData of chainDataMap.values()) {
        chainData.residues.sort((a, b) => a.auth_seq_id - b.auth_seq_id);
    }

    return Array.from(chainDataMap.values());
}

// ============================================================
// Ligand Neighborhood Extraction
// ============================================================

interface LigandInstance {
    compId: string;
    authAsymId: string;
    authSeqId: number;
}

function findLigandInstances(structure: Structure): LigandInstance[] {
    const ligands: LigandInstance[] = [];
    const seen = new Set<string>();

    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;

        const loc = StructureElement.Location.create(structure, unit, unit.elements[0]);

        for (let i = 0; i < unit.elements.length; i++) {
            loc.element = unit.elements[i];

            const compId = StructureProperties.atom.auth_comp_id(loc);

            if (isProteinResidue(compId)) continue;
            if (SKIP_LIGANDS.has(compId.toUpperCase())) continue;

            const authAsymId = StructureProperties.chain.auth_asym_id(loc);
            const authSeqId = StructureProperties.residue.auth_seq_id(loc);

            const key = `${compId}_${authAsymId}_${authSeqId}`;
            if (seen.has(key)) continue;
            seen.add(key);

            ligands.push({ compId, authAsymId, authSeqId });
        }
    }

    return ligands;
}

function extractLigandNeighborhood(
    structure: Structure,
    ligand: LigandInstance
): SimplifiedLigandNeighborhood | null {
    try {
        const ligandQuery = MS.struct.generator.atomGroups({
            'residue-test': MS.core.logic.and([
                MS.core.rel.eq([MS.ammp('auth_seq_id'), ligand.authSeqId]),
                MS.core.rel.eq([MS.ammp('auth_comp_id'), ligand.compId])
            ]),
            'chain-test': MS.core.rel.eq([MS.ammp('auth_asym_id'), ligand.authAsymId])
        });

        const surroundingsExpr = MS.struct.modifier.includeSurroundings({
            0: ligandQuery,
            radius: 5,
            'as-whole-residues': true
        });

        const surroundingsQuery = compile(surroundingsExpr);
        const surroundingsSelection = surroundingsQuery(new QueryContext(structure));
        const neighborhoodStructure = StructureSelection.unionStructure(surroundingsSelection);

        if (neighborhoodStructure.elementCount === 0) {
            return null;
        }

        const neighborhoodResidues: NeighborhoodResidue[] = [];
        const seenRes = new Set<string>();

        for (const unit of neighborhoodStructure.units) {
            if (!Unit.isAtomic(unit)) continue;

            const loc = StructureElement.Location.create(neighborhoodStructure, unit, unit.elements[0]);

            for (let i = 0; i < unit.elements.length; i++) {
                loc.element = unit.elements[i];

                const authAsymId = StructureProperties.chain.auth_asym_id(loc);
                const authSeqId = StructureProperties.residue.auth_seq_id(loc);
                const compId = StructureProperties.atom.auth_comp_id(loc);

                if (authAsymId === ligand.authAsymId &&
                    authSeqId === ligand.authSeqId &&
                    compId === ligand.compId) {
                    continue;
                }

                if (!isProteinResidue(compId)) continue;

                const key = `${authAsymId}:${authSeqId}`;
                if (seenRes.has(key)) continue;
                seenRes.add(key);

                neighborhoodResidues.push({
                    auth_asym_id: authAsymId,
                    observed_index: authSeqId,
                    comp_id: compId
                });
            }
        }

        neighborhoodResidues.sort((a, b) => {
            const chainCmp = a.auth_asym_id.localeCompare(b.auth_asym_id);
            if (chainCmp !== 0) return chainCmp;
            return a.observed_index - b.observed_index;
        });

        return {
            ligand_comp_id: ligand.compId,
            ligand_auth_asym_id: ligand.authAsymId,
            ligand_auth_seq_id: ligand.authSeqId,
            neighborhood_residues: neighborhoodResidues
        };

    } catch (e) {
        console.error(`Error extracting neighborhood for ${ligand.compId}_${ligand.authAsymId}: ${e}`);
        return null;
    }
}

// ============================================================
// Main
// ============================================================

async function runExtraction(cifPath: string, rcsbId: string, outputPath: string) {
    try {
        console.error(`Loading: ${cifPath}`);
        const fileContent = await fs.readFile(cifPath, 'utf8');

        const structure = await loadStructureFromCif(fileContent);
        console.error(`  Loaded structure with ${structure.units.length} units`);

        console.error(`Extracting sequences...`);
        const sequences = extractObservedSequences(structure);
        console.error(`  Found ${sequences.length} polymer chains`);

        console.error(`Extracting ligand neighborhoods...`);
        const ligandInstances = findLigandInstances(structure);
        console.error(`  Found ${ligandInstances.length} ligand instances`);

        const ligandNeighborhoods: SimplifiedLigandNeighborhood[] = [];
        for (const lig of ligandInstances) {
            const neighborhood = extractLigandNeighborhood(structure, lig);
            if (neighborhood && neighborhood.neighborhood_residues.length > 0) {
                ligandNeighborhoods.push(neighborhood);
                console.error(`    ${lig.compId}_${lig.authAsymId}: ${neighborhood.neighborhood_residues.length} nearby residues`);
            }
        }

        const result: ExtractionResult = {
            rcsb_id: rcsbId.toUpperCase(),
            sequences,
            ligand_neighborhoods: ligandNeighborhoods
        };

        await fs.mkdir(path.dirname(outputPath), { recursive: true });
        await fs.writeFile(outputPath, JSON.stringify(result, null, 2));
        console.error(`Success! Written to ${outputPath}`);

    } catch (e) {
        console.error("Error:", e);
        process.exit(1);
    }
}

const args = process.argv.slice(2);
if (args.length < 3) {
    console.error("Usage: npx tsx extract_structure_data.tsx <cif_path> <rcsb_id> <output_json>");
    process.exit(1);
}

runExtraction(args[0], args[1], args[2]);