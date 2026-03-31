// scripts_and_artifacts/render_thumbnail.tsx
//
// Headless Molstar renderer for structure thumbnails.
// Called by the ETL pipeline via molstar_bridge.py.
//
// Usage:
//   tsx render_thumbnail.tsx <cif_path> <rcsb_id> <profile_json_path> <output_path>
//
// Produces a transparent PNG at 480x320 with ghost-mode colors
// (beige alpha, gray-blue beta, 55% transparency) and vibrant MAPs/ligands.

// Set up minimal DOM environment before any Molstar imports.
// Molstar's Canvas3D.create calls syncCanvasBackground(canvas, ...) which
// accesses canvas.style. In the headless path, HeadlessScreenshotHelper
// doesn't provide a canvas to Canvas3D.create, so it's undefined.
// We patch Canvas3DContext to inject a mock canvas element.
import { JSDOM } from 'jsdom';
const dom = new JSDOM('<!DOCTYPE html><html><body><canvas id="mock"></canvas></body></html>', { pretendToBeVisual: true });
(globalThis as any).window = dom.window;
(globalThis as any).document = dom.window.document;
(globalThis as any).self = dom.window;
(globalThis as any).HTMLElement = dom.window.HTMLElement;
(globalThis as any).HTMLCanvasElement = dom.window.HTMLCanvasElement;
(globalThis as any).Node = dom.window.Node;
(globalThis as any).getComputedStyle = dom.window.getComputedStyle;
if (!(globalThis as any).navigator) (globalThis as any).navigator = dom.window.navigator;
if (!(globalThis as any).DOMParser) (globalThis as any).DOMParser = dom.window.DOMParser;
if (!(globalThis as any).XMLSerializer) (globalThis as any).XMLSerializer = dom.window.XMLSerializer;
if (!(globalThis as any).Element) (globalThis as any).Element = dom.window.Element;
if (!(globalThis as any).DocumentFragment) (globalThis as any).DocumentFragment = dom.window.DocumentFragment;
if (!(globalThis as any).requestAnimationFrame) (globalThis as any).requestAnimationFrame = (cb: Function) => setTimeout(cb, 0);
if (!(globalThis as any).cancelAnimationFrame) (globalThis as any).cancelAnimationFrame = (id: number) => clearTimeout(id);

import gl from 'gl';
import pngjs from 'pngjs';
import jpegjs from 'jpeg-js';

// Monkey-patch Canvas3D.create to inject a mock canvas element when none is provided.
// HeadlessScreenshotHelper in molstar 5.6.x doesn't pass canvas to Canvas3D.create,
// but syncCanvasBackground accesses canvas.style which crashes without this.
import { Canvas3D } from 'molstar/lib/mol-canvas3d/canvas3d';
const _origCreate = Canvas3D.create;
Canvas3D.create = function(ctx: any, ...rest: any[]) {
    if (!ctx.canvas) {
        ctx.canvas = document.createElement('canvas');
        document.body.appendChild(ctx.canvas);
    }
    return (_origCreate as any).call(this, ctx, ...rest);
};

import { HeadlessPluginContext } from 'molstar/lib/mol-plugin/headless-plugin-context';
import { DefaultPluginSpec } from 'molstar/lib/mol-plugin/spec';
import { StructureRepresentationPresetProvider } from 'molstar/lib/mol-plugin-state/builder/structure/representation-preset';
import { ParamDefinition as PD } from 'molstar/lib/mol-util/param-definition';
import { Structure } from 'molstar/lib/mol-model/structure';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { StateObjectRef } from 'molstar/lib/mol-state';
import { PluginStateObject } from 'molstar/lib/mol-plugin-state/objects';
import { Color } from 'molstar/lib/mol-util/color';
import { structureLayingTransform, changeCameraRotation } from 'molstar/lib/mol-plugin-state/manager/focus-camera/orient-axes';
import { setStructureTransparency } from 'molstar/lib/mol-plugin-state/helpers/structure-transparency';
import * as fs from 'fs';
import * as path from 'path';
import { setFSModule } from 'molstar/lib/mol-util/data-source';
setFSModule(fs);

// ============================================================
// Configuration
// ============================================================

const IMAGE_WIDTH = 640;
const IMAGE_HEIGHT = 420;

// Ghost colors for alpha/beta tubulin (muted, catalogue-friendly)
const GHOST_COLORS: Record<string, Color> = {
    tubulin_alpha: Color(0xD4C4A8),  // warm beige/taupe
    tubulin_beta:  Color(0xB8C4D0),  // cool grayish-blue
};

const GHOST_TRANSPARENCY = 0.55;

// Vibrant colors for non-ghost chains
const TUBULIN_COLORS: Record<string, Color> = {
    tubulin_alpha:   Color(0x4F92E8),
    tubulin_beta:    Color(0xE07850),
    tubulin_gamma:   Color(0xA07CC0),
    tubulin_delta:   Color(0x5EAB70),
    tubulin_epsilon: Color(0xD4C060),
};

const MAP_COLORS: Record<string, Color> = {
    map_eb_family:            Color(0x00CED1),
    map_camsap1:              Color(0x20B2AA),
    map_camsap2:              Color(0x48D1CC),
    map_camsap3:              Color(0x40E0D0),
    map_kinesin13:            Color(0xFF1493),
    map_katanin_p60:          Color(0xFF69B4),
    map_spastin:              Color(0xDA70D6),
    map_tau:                  Color(0xFF8C00),
    map_map2:                 Color(0xE67E22),
    map_doublecortin:         Color(0xD35400),
    map_gcp2_3:               Color(0x1F618D),
    map_gcp4:                 Color(0x2874A6),
    map_gcp5_6:               Color(0x2E86C1),
    map_vash_detyrosinase:    Color(0x27AE60),
    map_atat1:                Color(0x2ECC71),
    map_ttll_glutamylase_long: Color(0xA9DFBF),
};

const DEFAULT_COLOR = Color(0xBDC3C7);

// Ligands to skip (crystallography artifacts)
const LIGAND_IGNORE_IDS = new Set([
    'EDO', 'GOL', 'MPD', 'PEG', 'PG4', 'PG0', 'PGE', '1PE', 'P6G',
    'DMS', 'SO4', 'PO4', 'ACT', 'MES', 'CIT', 'BME',
    'HOH', 'DOD', 'WAT',
]);

// Ligand colors (subset of most common)
const LIGAND_COLORS: Record<string, Color> = {
    GTP: Color(0x2979FF), GDP: Color(0xFF9100), GCP: Color(0x82B1FF),
    GSP: Color(0x448AFF), ANP: Color(0xB3E5FC), ACP: Color(0xFFE082),
    TXL: Color(0x00C853), VLB: Color(0xF50057), CN2: Color(0xFF6D00),
    COL: Color(0xFF9E80), MG:  Color(0x00E5FF), ZN:  Color(0xEEFF41),
};

// Postprocessing (outlines + ambient occlusion, no shadows)
const POSTPROCESSING = {
    outline: {
        name: 'on' as const,
        params: {
            scale: 1,
            color: Color(0x555555),
            threshold: 0.12,
            includeTransparent: true,
        },
    },
    occlusion: {
        name: 'on' as const,
        params: {
            multiScale: { name: 'off' as const, params: {} },
            radius: 5,
            bias: 0.8,
            blurKernelSize: 15,
            blurDepthBias: 0.5,
            samples: 32,
            resolutionScale: 1,
            color: Color(0x000000),
        },
    },
    shadow: { name: 'off' as const, params: {} },
};

// ============================================================
// Helpers
// ============================================================

type Classification = Record<string, string>;

function classificationFromProfile(profile: any): Classification {
    const classification: Classification = {};
    if (!profile?.polypeptides) return classification;
    for (const poly of profile.polypeptides) {
        const entity = profile.entities?.[poly.entity_id];
        if (entity?.family) {
            classification[poly.auth_asym_id] = entity.family;
        }
    }
    return classification;
}

function isAlphaBeta(family: string | undefined): boolean {
    return family === 'tubulin_alpha' || family === 'tubulin_beta';
}

function getChainColor(family: string | undefined, ghost: boolean): Color {
    if (!family) return DEFAULT_COLOR;
    if (ghost && isAlphaBeta(family)) {
        return GHOST_COLORS[family] ?? DEFAULT_COLOR;
    }
    return TUBULIN_COLORS[family] ?? MAP_COLORS[family] ?? DEFAULT_COLOR;
}

function getLigandColor(compId: string): Color {
    return LIGAND_COLORS[compId] ?? Color(0x94A3B8);
}

interface LigandInstance {
    compId: string;
    auth_asym_id: string;
    auth_seq_id: number;
    uniqueKey: string;
}

function getLigandInstances(structure: Structure): LigandInstance[] {
    const instances: LigandInstance[] = [];
    const seen = new Set<string>();
    const { _rowCount, auth_asym_id } = structure.model.atomicHierarchy.chains;

    for (let cI = 0; cI < _rowCount; cI++) {
        const eI = structure.model.atomicHierarchy.index.getEntityFromChain(cI);
        if (structure.model.entities.data.type.value(eI) !== 'non-polymer') continue;

        const chainId = auth_asym_id.value(cI);
        const residueOffsets = structure.model.atomicHierarchy.chainAtomSegments.offsets;
        const residueIndex = structure.model.atomicHierarchy.residueAtomSegments;

        const atomStart = residueOffsets[cI];
        const atomEnd = residueOffsets[cI + 1];

        for (let aI = atomStart; aI < atomEnd; aI++) {
            const rI = residueIndex.index[aI];
            const compId = structure.model.atomicHierarchy.atoms.label_comp_id.value(aI);
            const seqId = structure.model.atomicHierarchy.residues.auth_seq_id.value(rI);
            const key = `${compId}_${chainId}_${seqId}`;
            if (seen.has(key)) continue;
            seen.add(key);
            instances.push({ compId, auth_asym_id: chainId, auth_seq_id: seqId, uniqueKey: key });
        }
    }
    return instances;
}

// ============================================================
// Thumbnail Preset - ghost mode built in
// ============================================================

const ThumbnailPreset = StructureRepresentationPresetProvider({
    id: 'thumbnail-ghost-preset',
    display: { name: 'Thumbnail Ghost', group: 'TubulinXYZ' },
    params: () => ({
        ...StructureRepresentationPresetProvider.CommonParams,
        pdbId: PD.Text(''),
        tubulinClassification: PD.Value<Classification>({}, { isHidden: true }),
    }),

    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const structure = structureCell.obj!.data as Structure;
        const { update } = StructureRepresentationPresetProvider.reprBuilder(plugin, params);

        const polymerRefs: { ref: string; chainId: string; family?: string }[] = [];

        // Polymer chains
        const { auth_asym_id } = structure.model.atomicHierarchy.chains;
        const chainCount = structure.model.atomicHierarchy.chains._rowCount;
        const seenChains = new Set<string>();

        for (let cI = 0; cI < chainCount; cI++) {
            const chainId = auth_asym_id.value(cI);
            if (seenChains.has(chainId)) continue;
            seenChains.add(chainId);

            const eI = structure.model.atomicHierarchy.index.getEntityFromChain(cI);
            if (structure.model.entities.data.type.value(eI) !== 'polymer') continue;

            const family = params.tubulinClassification[chainId];
            const color = getChainColor(family, true);

            const component = await plugin.builders.structure.tryCreateComponentFromExpression(
                structureCell,
                MS.struct.generator.atomGroups({
                    'chain-test': MS.core.rel.eq([
                        MS.struct.atomProperty.macromolecular.auth_asym_id(),
                        chainId,
                    ]),
                }),
                `${params.pdbId}_${chainId}`,
                { label: `${family || 'Polymer'} (${chainId})` },
            );

            if (component) {
                await plugin.builders.structure.representation.addRepresentation(component, {
                    type: 'cartoon',
                    color: 'uniform',
                    colorParams: { value: color },
                });
                polymerRefs.push({ ref: component.ref, chainId, family });
            }
        }

        // Ligands
        const ligandInstances = getLigandInstances(structure);
        for (const instance of ligandInstances) {
            if (LIGAND_IGNORE_IDS.has(instance.compId)) continue;

            const component = await plugin.builders.structure.tryCreateComponentFromExpression(
                structureCell,
                MS.struct.generator.atomGroups({
                    'residue-test': MS.core.logic.and([
                        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), instance.compId]),
                        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), instance.auth_asym_id]),
                        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), instance.auth_seq_id]),
                    ]),
                }),
                `${params.pdbId}_${instance.uniqueKey}`,
                { label: `Ligand ${instance.compId}` },
            );

            if (component) {
                await plugin.builders.structure.representation.addRepresentation(component, {
                    type: 'ball-and-stick',
                    color: 'uniform',
                    colorParams: { value: getLigandColor(instance.compId) },
                    typeParams: { sizeFactor: 0.3 },
                });
            }
        }

        await update.commit({ revertOnError: true });

        // Apply ghost transparency to alpha/beta chains
        const hierarchy = plugin.managers.structure.hierarchy.current;
        if (hierarchy.structures.length > 0) {
            for (const pr of polymerRefs) {
                if (!isAlphaBeta(pr.family)) continue;

                const chainComponents = hierarchy.structures[0].components.filter(
                    c => c.cell.transform.ref === pr.ref,
                );
                if (chainComponents.length === 0) continue;

                try {
                    await setStructureTransparency(
                        plugin,
                        chainComponents,
                        GHOST_TRANSPARENCY,
                        async (structure) => Structure.toStructureElementLoci(structure),
                    );
                } catch (_) {
                    // Non-fatal: thumbnail still renders, just without transparency
                }
            }
        }

        return {};
    },
});

// ============================================================
// Main
// ============================================================

async function main() {
    const [cifPath, rcsbId, profilePath, outputPath] = process.argv.slice(2);

    if (!cifPath || !rcsbId || !profilePath || !outputPath) {
        console.error('Usage: tsx render_thumbnail.tsx <cif_path> <rcsb_id> <profile_json_path> <output_path>');
        process.exit(1);
    }

    // Read profile and build classification
    const profileRaw = JSON.parse(fs.readFileSync(profilePath, 'utf-8'));
    const classification = classificationFromProfile(profileRaw);

    // Create headless plugin
    const plugin = new HeadlessPluginContext(
        { gl, pngjs, 'jpeg-js': jpegjs },
        DefaultPluginSpec(),
        { width: IMAGE_WIDTH, height: IMAGE_HEIGHT },
        {
            imagePass: {
                transparentBackground: true,
                cameraHelper: { axes: { name: 'off', params: {} } },
                multiSample: { mode: 'on', sampleLevel: 4 },
            },
        },
    );

    await plugin.init();
    plugin.builders.structure.representation.registerPreset(ThumbnailPreset);

    try {
        // Load structure from local CIF file
        const data = await plugin.builders.data.download({
            url: `file://${path.resolve(cifPath)}`,
            isBinary: false,
        });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);

        // Apply thumbnail preset
        await plugin.builders.structure.representation.applyPreset(
            structure,
            'thumbnail-ghost-preset',
            { pdbId: rcsbId, tubulinClassification: classification },
        );

        // Set ignoreLight on all representations
        const reprUpdate = plugin.state.data.build();
        const representations = plugin.state.data.selectQ(q =>
            q.ofType(PluginStateObject.Molecule.Structure.Representation3D),
        );
        for (const repr of representations) {
            reprUpdate.to(repr).update(old => {
                if (old.type?.params) old.type.params.ignoreLight = true;
            });
        }
        await reprUpdate.commit();

        plugin.managers.structure.component.setOptions({
            ...plugin.managers.structure.component.state.options,
            ignoreLight: true,
        });

        // Orient camera
        const structCells = plugin.state.data.selectQ(q =>
            q.ofType(PluginStateObject.Molecule.Structure),
        );
        const structures = structCells
            .filter(cell => cell.obj && !cell.obj.data.parent)
            .map(cell => cell.obj!.data);

        if (structures.length > 0 && plugin.canvas3d) {
            const { rotation } = structureLayingTransform(structures);
            const currentSnapshot = plugin.canvas3d.camera.getSnapshot();
            const newSnapshot = changeCameraRotation(currentSnapshot, rotation);
            plugin.canvas3d.camera.setState(newSnapshot);
        }

        // Save
        fs.mkdirSync(path.dirname(outputPath), { recursive: true });
        await plugin.saveImage(
            outputPath,
            { width: IMAGE_WIDTH, height: IMAGE_HEIGHT },
            POSTPROCESSING,
            'png',
            90,
        );

        console.log(`OK ${outputPath}`);
    } catch (err) {
        console.error(`FAIL ${rcsbId}: ${err}`);
        process.exit(1);
    } finally {
        plugin.dispose();
    }
}

main();
