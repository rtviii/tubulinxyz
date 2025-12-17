import { HeadlessPluginContext } from 'molstar/lib/mol-plugin/headless-plugin-context';
import { DefaultPluginSpec } from 'molstar/lib/mol-plugin/spec';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { StructureElement, StructureProperties, StructureSelection, QueryContext, Structure, Unit } from 'molstar/lib/mol-model/structure';
import { compile } from 'molstar/lib/mol-script/runtime/query/compiler';
import { InteractionsProvider } from 'molstar/lib/mol-model-props/computed/interactions';
import { interactionTypeLabel } from 'molstar/lib/mol-model-props/computed/interactions/common';
import { Features } from 'molstar/lib/mol-model-props/computed/interactions/features';
import { SyncRuntimeContext } from 'molstar/lib/mol-task/execution/synchronous';
import { AssetManager } from 'molstar/lib/mol-util/assets';
import { OrderedSet } from 'molstar/lib/mol-data/int';

import gl from 'gl';
import pngjs from 'pngjs';
import jpegjs from 'jpeg-js';

import * as fs from 'fs/promises';
import * as path from 'path';

// --- Compact Output Types ---

// [auth_asym_id, auth_seq_id, auth_comp_id]
type ResidueTuple = [string, number, string];

// [auth_asym_id, auth_seq_id, auth_comp_id, atom_id, isLigand]
type ParticipantTuple = [string, number, string, string, boolean];

interface Interaction {
    type: string;
    participants: [ParticipantTuple, ParticipantTuple];
}

interface BindingSiteResult {
    ligand: ResidueTuple;
    interactions: Interaction[];
    neighborhood: ResidueTuple[];
}

// --- Internal Types (for processing) ---

interface ResidueInfo {
    auth_asym_id: string;
    auth_seq_id: number;
    auth_comp_id: string;
}

interface AtomInfo extends ResidueInfo {
    atom_id: string;
}

// --- ALL INTERACTIONS ENABLED ---
const allInteractionsParams = {
    providers: {
        'ionic': { name: 'on', params: {} },
        'cation-pi': { name: 'on', params: {} },
        'pi-stacking': { name: 'on', params: {} },
        'hydrogen-bonds': { name: 'on', params: {} },
        'halogen-bonds': { name: 'on', params: {} },
        'hydrophobic': { name: 'on', params: {} },
        'metal-coordination': { name: 'on', params: {} },
        'weak-hydrogen-bonds': { name: 'on', params: {} },
    }
};

// --- Helper Functions ---

function getResidueInfo(loc: StructureElement.Location): ResidueInfo {
    return {
        auth_asym_id: StructureProperties.chain.auth_asym_id(loc),
        auth_seq_id: StructureProperties.residue.auth_seq_id(loc),
        auth_comp_id: StructureProperties.atom.auth_comp_id(loc)
    };
}

function toResidueTuple(info: ResidueInfo): ResidueTuple {
    return [info.auth_asym_id, info.auth_seq_id, info.auth_comp_id];
}

function toParticipantTuple(atom: AtomInfo, isLigand: boolean): ParticipantTuple {
    return [atom.auth_asym_id, atom.auth_seq_id, atom.auth_comp_id, atom.atom_id, isLigand];
}

function getAtomInfoFromFeature(
    structure: Structure,
    unit: Unit,
    featureIndex: number,
    features: Features
): AtomInfo {
    const atomIndex = features.members[features.offsets[featureIndex]];
    const elementIndex = unit.elements[atomIndex];
    const loc = StructureElement.Location.create(structure, unit, elementIndex);

    return {
        auth_asym_id: StructureProperties.chain.auth_asym_id(loc),
        auth_seq_id: StructureProperties.residue.auth_seq_id(loc),
        auth_comp_id: StructureProperties.atom.auth_comp_id(loc),
        atom_id: StructureProperties.atom.label_atom_id(loc)
    };
}

function findLigandInStructure(structure: Structure, ligandInfo: ResidueInfo): { unit: Unit, indices: OrderedSet } | null {
    for (const unit of structure.units) {
        const atomIndices: number[] = [];
        const loc = StructureElement.Location.create(structure, unit, unit.elements[0]);

        for (let i = 0; i < unit.elements.length; i++) {
            loc.element = unit.elements[i];
            if (StructureProperties.chain.auth_asym_id(loc) === ligandInfo.auth_asym_id &&
                StructureProperties.residue.auth_seq_id(loc) === ligandInfo.auth_seq_id &&
                StructureProperties.atom.auth_comp_id(loc) === ligandInfo.auth_comp_id) {
                atomIndices.push(i);
            }
        }

        if (atomIndices.length > 0) {
            return {
                unit,
                indices: OrderedSet.ofSortedArray(atomIndices.sort((a, b) => a - b))
            };
        }
    }
    return null;
}

function extractNeighborhoodResidues(structure: Structure, ligandInfo: ResidueInfo): ResidueTuple[] {
    const seen = new Set<string>();
    const residues: ResidueTuple[] = [];

    for (const unit of structure.units) {
        const loc = StructureElement.Location.create(structure, unit, unit.elements[0]);

        for (let i = 0; i < unit.elements.length; i++) {
            loc.element = unit.elements[i];
            const info = getResidueInfo(loc);

            if (info.auth_asym_id === ligandInfo.auth_asym_id &&
                info.auth_seq_id === ligandInfo.auth_seq_id &&
                info.auth_comp_id === ligandInfo.auth_comp_id) continue;

            const key = `${info.auth_asym_id}:${info.auth_seq_id}:${info.auth_comp_id}`;
            if (!seen.has(key)) {
                seen.add(key);
                residues.push(toResidueTuple(info));
            }
        }
    }

    return residues;
}

function extractInteractionsFromStructure(
    structure: Structure,
    ligandUnitId: number,
    ligandIndices: OrderedSet
): Interaction[] {
    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) {
        console.log('  WARNING: No interactions object found');
        return [];
    }

    const results: Interaction[] = [];
    const seen = new Set<string>();
    const { contacts: interUnitContacts, unitsFeatures } = interactions;

    console.log(`  Inter-unit contacts: ${interUnitContacts.edgeCount}`);

    const isAtomInLigand = (unitId: number, featureIndex: number, features: Features): boolean => {
        if (unitId !== ligandUnitId) return false;
        const atomIndex = features.members[features.offsets[featureIndex]];
        return OrderedSet.has(ligandIndices, atomIndex);
    };

    for (const bond of interUnitContacts.edges) {
        const unitA = structure.unitMap.get(bond.unitA);
        const unitB = structure.unitMap.get(bond.unitB);
        if (!unitA || !unitB) continue;

        const featuresA = unitsFeatures.get(bond.unitA);
        const featuresB = unitsFeatures.get(bond.unitB);
        if (!featuresA || !featuresB) continue;

        const isALigand = isAtomInLigand(bond.unitA, bond.indexA, featuresA);
        const isBLigand = isAtomInLigand(bond.unitB, bond.indexB, featuresB);

        if ((isALigand && !isBLigand) || (!isALigand && isBLigand)) {
            const atomA = getAtomInfoFromFeature(structure, unitA, bond.indexA, featuresA);
            const atomB = getAtomInfoFromFeature(structure, unitB, bond.indexB, featuresB);

            const keyParts = [
                `${atomA.auth_asym_id}:${atomA.auth_seq_id}:${atomA.atom_id}`,
                `${atomB.auth_asym_id}:${atomB.auth_seq_id}:${atomB.atom_id}`
            ].sort();
            const key = `${bond.props.type}:${keyParts[0]}:${keyParts[1]}`;

            if (seen.has(key)) continue;
            seen.add(key);

            results.push({
                type: interactionTypeLabel(bond.props.type),
                participants: [
                    toParticipantTuple(atomA, isALigand),
                    toParticipantTuple(atomB, isBLigand)
                ]
            });
        }
    }

    return results;
}

// --- Main Processor ---

async function runExtraction(pdbPath: string, ligandCompId: string, authAsymId: string, outputPath: string) {
    const externalModules = { gl, pngjs, 'jpeg-js': jpegjs };
    const spec = DefaultPluginSpec();

    const plugin = new HeadlessPluginContext(
        externalModules,
        spec,
        { width: 800, height: 600 }
    );
    await plugin.init();

    try {
        console.log(`Reading file: ${pdbPath}`);

        const fileContent = await fs.readFile(pdbPath, 'utf8');

        const data = await plugin.builders.data.rawData({
            data: fileContent,
            label: 'CIF Data'
        });

        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structureNode = await plugin.builders.structure.createStructure(model);

        const structure = structureNode.data;
        if (!structure) throw new Error("Could not retrieve structure data");

        const ligandQuery = MS.struct.generator.atomGroups({
            'residue-test': MS.core.rel.eq([MS.ammp('auth_comp_id'), ligandCompId]),
            'chain-test': MS.core.rel.eq([MS.ammp('auth_asym_id'), authAsymId]),
            'group-by': MS.ammp('residueKey')
        });

        const ligandSelection = compile(ligandQuery)(new QueryContext(structure));
        const ligandLoci = StructureSelection.toLociWithSourceUnits(ligandSelection);

        if (ligandLoci.elements.length === 0) {
            console.warn(`No ligand found with comp_id=${ligandCompId} and auth_asym_id=${authAsymId}`);
        }

        const findings: BindingSiteResult[] = [];

        for (const element of ligandLoci.elements) {
            const unit = element.unit;

            const residueIndicesMap = new Map<number, number[]>();
            const location = StructureElement.Location.create(structure, unit, unit.elements[0]);

            OrderedSet.forEach(element.indices, (index) => {
                location.element = unit.elements[index];
                const residueIndex = StructureProperties.residue.key(location);
                if (!residueIndicesMap.has(residueIndex)) residueIndicesMap.set(residueIndex, []);
                residueIndicesMap.get(residueIndex)!.push(index);
            });

            for (const [_resKey, atomIndicesArr] of residueIndicesMap.entries()) {
                location.element = unit.elements[atomIndicesArr[0]];
                const currentLigandInfo = getResidueInfo(location);

                console.log(`Processing: ${currentLigandInfo.auth_comp_id} | Chain: ${currentLigandInfo.auth_asym_id} | Res: ${currentLigandInfo.auth_seq_id}`);

                const specificLigandQuery = MS.struct.generator.atomGroups({
                    'residue-test': MS.core.logic.and([
                        MS.core.rel.eq([MS.ammp('auth_seq_id'), currentLigandInfo.auth_seq_id]),
                        MS.core.rel.eq([MS.ammp('auth_comp_id'), currentLigandInfo.auth_comp_id])
                    ]),
                    'chain-test': MS.core.rel.eq([MS.ammp('auth_asym_id'), currentLigandInfo.auth_asym_id])
                });

                const surroundingsExpression = MS.struct.modifier.includeSurroundings({
                    0: specificLigandQuery,
                    radius: 5,
                    'as-whole-residues': true
                });

                const surroundingsQuery = compile(surroundingsExpression);
                const surroundingsSelection = surroundingsQuery(new QueryContext(structure));
                const neighborhoodStructure = StructureSelection.unionStructure(surroundingsSelection);

                const ligandInNeighborhood = findLigandInStructure(neighborhoodStructure, currentLigandInfo);
                if (!ligandInNeighborhood) {
                    console.warn(`  Could not find ligand in neighborhood structure`);
                    continue;
                }

                console.log(`  Found ligand unit ${ligandInNeighborhood.unit.id} with ${OrderedSet.size(ligandInNeighborhood.indices)} atoms`);

                const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };
                const allInteractionsParams = {
                    providers: {
                        'ionic': { name: 'on', params: { distanceMax: 5.0 } },
                        'cation-pi': { name: 'on', params: { distanceMax: 6.0 } },
                        'pi-stacking': { name: 'on', params: { distanceMax: 5.5, offsetMax: 2.0, angleDevMax: 30 } },
                        'hydrogen-bonds': {
                            name: 'on', params: {
                                distanceMax: 3.5,
                                angleMax: 45,
                                water: false,
                                sulfurDistanceMax: 4.1
                            }
                        },
                        'halogen-bonds': { name: 'on', params: { distanceMax: 4.0, angleMax: 30 } },
                        'hydrophobic': { name: 'on', params: { distanceMax: 4.0 } },
                        'metal-coordination': { name: 'on', params: { distanceMax: 2.5 } },
                        'weak-hydrogen-bonds': {
                            name: 'on', params: {
                                distanceMax: 3.5,
                                angleMax: 45
                            }
                        },
                    }
                };
                await InteractionsProvider.attach(ctx, neighborhoodStructure, allInteractionsParams, true);

                const interactions = extractInteractionsFromStructure(
                    neighborhoodStructure,
                    ligandInNeighborhood.unit.id,
                    ligandInNeighborhood.indices
                );
                const neighborhood = extractNeighborhoodResidues(neighborhoodStructure, currentLigandInfo);

                console.log(`  Found ${interactions.length} ligand-polymer interactions`);

                findings.push({
                    ligand: toResidueTuple(currentLigandInfo),
                    interactions,
                    neighborhood
                });
            }
        }

        await fs.mkdir(path.dirname(outputPath), { recursive: true });
        await fs.writeFile(outputPath, JSON.stringify(findings, null, 2));
        console.log(`Success! Data written to ${outputPath}`);

    } catch (e) {
        console.error("Error processing structure:", e);
        process.exit(1);
    } finally {
        plugin.dispose();
    }
}

const args = process.argv.slice(2);
if (args.length < 4) {
    console.log("Usage: npx tsx extract-interactions.tsx <path_to_cif> <ligand_comp_id> <auth_asym_id> <output_json>");
    process.exit(1);
}

runExtraction(args[0], args[1], args[2], args[3]);