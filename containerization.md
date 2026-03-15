Ok so i have this tubulinxyz application that i'd like to containerize and hand off to the people that are going to be hosting it for the foreseaable future. The whole app is basically a fastapi + neo4j backend and a nextjs frontend. 

It's going to be an institutional network/VM first while they check that there are no security holes in my system, and then if all is well they will open some ports to the world for other researchers to be able to interact with it.

I'm not sure what's the best way for us to do this and to check that everything works, but basically i want to have a self-sufficient self-healing system eventually that regenerates its own data, but for now i just want to send the basic dockerfiles to the guy so we can start the process of building/deploying the application on their infra.

In particular, i've done this kinda thing once with my other application (ribosome.xyz), but in that case only the backend and the neo4j db were containerized..Then the admin has provided me with a VM on which i would run an nginx server on "bare metal" that would talk to my backend and the DB. All the while the frontend itself -- which is also an nextjs application -- is running uncontainerized as a `pm2` process. That's the best i've come up with but it works quite well so far...


The nextjs app communicates with the fastapi via a codegen'ed rtk-query api that is generated from the fastapi's openapi schema. The FASTAPI in turn hits the neo4j database for whatever the frontend needs via a set of pretty general endpoints.

Let me show you both projects' layouts and we can go from there..

```
(venv) ᢹ saeta.rtviii[ dev/tubulinxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta|npet|*.mdx|*.ts.map|*.d.ts|nightingale|NPET2'
.
├── Carsten_report.md
├── cli.py
├── containerization.md
├── data
│   ├── alpha_tubulin
│   │   ├── alpha_modifications.json
│   │   └── alpha_mutations.json
│   ├── beta_tubulin
│   ├── delta_tubulin
│   ├── epsilon_tubulin
│   ├── gamma_tubulin
│   ├── hmms
│   │   ├── classification_cache
│   │   ├── maps
│   │   └── tubulin
│   ├── map_sequences
│   ├── maxim_data
│   │   ├── 7sj7_tails.pdb
│   │   ├── 7sj7_with_metadata.cif
│   │   ├── add_computed_res_annotation.py
│   │   ├── curved_hum_wt_GDP_6S8L.pdb
│   │   ├── fold_htuba1a_model_0.cif
│   │   ├── fold_htuba1a.zip
│   │   ├── mini_mt_patch_hum_wt_GDP_7SJ7.pdb
│   │   └── straight_hum_wt_GDP_7SJ7.pdb
│   └── sequences
│       ├── maps
│       └── tubulin
├── dirtocontext.py
├── docs
├── indexing_report.md
├── lib
│   ├── etl
│   │   ├── assets.py
│   │   ├── augmentation.py
│   │   ├── classification.py
│   │   ├── collector.py
│   │   ├── constants.py
│   │   ├── libtax.py
│   │   ├── molstar_bridge.py
│   │   └── sequence_alignment.py
│   ├── ingestion_logs
│   └── types.py
├── muscle3.8.1
├── neo4j_schema.json
├── neo4j_tubxz
│   ├── db_lib_builder.py
│   ├── db_lib_reader.py
│   ├── models.py
│   ├── node_binding_site.py
│   ├── node_ligand.py
│   ├── node_modification.py
│   ├── node_phylogeny.py
│   ├── node_polymer.py
│   ├── node_structure.py
│   ├── node_variant.py
│   ├── query_builder.py
│   ├── structure_query_builder.py
│   └── test_one_struct.py
├── notes
│   ├── 0_general_context.md
│   ├── 1_ligand_census.md
│   ├── 2_hmm_classifier_eval_report.txt
│   ├── 3_tubulin_classes_seq_clustering.md
│   ├── 4_MAP_report.md
│   ├── doc_hmm_building.md
│   ├── doc_IndexingModel.md
│   ├── indexing_pipeline_summary.md
│   ├── journal.pone.0295279-5.pdf
│   ├── lig_seq_id_ingestion.md
│   ├── MINUTES_dec19.md
│   ├── MINUTES_jan_26th.md
│   ├── MINUTES_november24th.md
│   ├── MINUTES_october22nd.md
│   ├── q_index_processing_pipleine.md
│   └── sequence_annotations_meetup.md
├── package.json
├── plan.md
├── README.md
├── requirements.txt
├── scripts_and_artifacts
│   ├── archive
│   │   ├── analyze_ligands.py
│   │   ├── tubulin_ligand_census.py
│   │   └── visualize_ligands.py
│   ├── extract_ixs.tsx
│   ├── extract_structure_data.tsx
│   ├── hmm_building
│   │   ├── analyze_cluster_composition.py
│   │   ├── eval_classifier.py
│   │   ├── fetch_mipmaps.py
│   │   ├── fetch_tubulin.py
│   │   ├── process_mipmaps.py
│   │   ├── process_tubulin.py
│   │   └── validate_pdb_structs.py
│   └── morisette_stuff
│       ├── morisette_alpha_beta_gamma_uniprot.md
│       ├── morisette_alpha.py
│       ├── mset_consensus.py
│       └── mset_parser.py
├── taxdump.tar.gz
└── verify_ingestion.py

26 directories, 77 files
(venv) ᢹ saeta.rtviii[ dev/tubulinxyz ]  pwd                                                                                                                     [main]
/Users/rtviii/dev/tubulinxyz
```

frontend:
```
(venv) ᢹ saeta.rtviii[ dev/fend_tubulinxyz ]  pwd                                                                                                                [main]
/Users/rtviii/dev/fend_tubulinxyz
(venv) ᢹ saeta.rtviii[ dev/fend_tubulinxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta|npet|*.mdx|*.ts.map|*.d.ts|nightingale|NPET2'
.
├── components.json
├── eslint.config.mjs
├── minutes_Feb27th.md
├── next.config.js
├── nonpolymer_chemical_ids.json
├── notes
│   └── ngl_doc.md
├── openapi-config.ts
├── package.json
├── PLAN_residue_detail_panel.md
├── plan.md
├── postcss.config.js
├── postcss.config.mjs
├── present_ligands.md
├── public
│   ├── __old_output
│   ├── file.svg
│   ├── globe.svg
│   ├── landing
│   └── output
├── q_chains_sync_ui.md
├── q_chains_sync.md
├── q_modes copy.md
├── q_modes.md
├── q_modes2.md
├── q_modes3_fix_colors.md
├── q_polymer_filtering_dialog.md
├── README.md
├── records.json
├── REPORT_next_refactors.md
├── scripts
│   ├── bulk-render.ts
│   ├── extract_ixs.tsx
│   └── extract_strucutre_data.tsx
├── src
│   ├── app
│   │   ├── debug
│   │   │   ├── LigandDebug.tsx
│   │   │   └── page.tsx
│   │   ├── globals.css
│   │   ├── landing
│   │   │   ├── LandingViewer.tsx
│   │   │   └── TubulinLandingViewer.tsx
│   │   ├── layout.tsx
│   │   ├── page.tsx
│   │   └── structures
│   │       ├── [rcsb_id]
│   │       │   └── page.tsx
│   │       ├── page.tsx
│   │       ├── structure_filters.tsx
│   │       └── StructureFiltersPanel.tsx
│   ├── components
│   │   ├── annotations
│   │   │   ├── LigandsPanel.tsx
│   │   │   └── VariantsPanel.tsx
│   │   ├── explorer
│   │   │   ├── CanonicalSiteSearch.tsx
│   │   │   ├── ExplorerPanel.tsx
│   │   │   ├── heatmapColors.ts
│   │   │   ├── questions
│   │   │   │   ├── useCanonicalBindingSite.ts
│   │   │   │   ├── useInterfaceContacts.ts
│   │   │   │   └── useNucleotideHighlight.ts
│   │   │   └── types.ts
│   │   ├── molstar
│   │   │   ├── coloring
│   │   │   │   ├── ColorschemeManager.ts
│   │   │   │   ├── schemes
│   │   │   │   │   ├── interactionScheme.ts
│   │   │   │   │   └── mutationScheme.ts
│   │   │   │   └── types.ts
│   │   │   ├── colors
│   │   │   │   ├── palette.ts
│   │   │   │   ├── preset_structure.tsx
│   │   │   │   ├── preset-helpers.ts
│   │   │   │   ├── schemes
│   │   │   │   │   └── ligandHitTheme.ts
│   │   │   │   └── tubulin-color-theme.tsx
│   │   │   ├── core
│   │   │   │   ├── MolstarViewer.ts
│   │   │   │   ├── queries.ts
│   │   │   │   └── types.ts
│   │   │   ├── labels
│   │   │   │   └── LabelManager.ts
│   │   │   ├── mstar.css
│   │   │   ├── overlay
│   │   │   │   └── ResidueInfoOverlay.tsx
│   │   │   ├── rendering
│   │   │   │   └── postprocessing-config.ts
│   │   │   ├── services
│   │   │   │   ├── MolstarInstance.ts
│   │   │   │   └── MolstarInstanceManager.tsx
│   │   │   ├── spec.tsx
│   │   │   └── state
│   │   │       ├── molstarInstancesSlice.ts
│   │   │       └── selectors.ts
│   │   ├── monomer
│   │   │   ├── AlignmentDialog.tsx
│   │   │   ├── AlignStructureForm.tsx
│   │   │   ├── ChainRow.tsx
│   │   │   ├── MonomerMSAPanel.tsx
│   │   │   ├── MonomerSidebar.tsx
│   │   │   └── PolymerBrowser.tsx
│   │   ├── msa
│   │   │   ├── AnnotationPanel.tsx
│   │   │   ├── index.ts
│   │   │   ├── MSALabels.tsx
│   │   │   ├── MSAToolbar.tsx
│   │   │   ├── ResizableMSAContainer.tsx
│   │   │   └── types.ts
│   │   ├── structure
│   │   │   └── StructureSidebar.tsx
│   │   └── ui
│   │       ├── badge.tsx
│   │       ├── button.tsx
│   │       ├── card.tsx
│   │       ├── checkbox.tsx
│   │       ├── collapsible.tsx
│   │       ├── FloatingNav.tsx
│   │       ├── input.tsx
│   │       ├── resizable.tsx
│   │       ├── separator.tsx
│   │       └── skeleton.tsx
│   ├── config.ts
│   ├── hooks
│   │   ├── useAnnotationVisibility.ts
│   │   ├── useChainAlignment.ts
│   │   ├── useChainFocusSync.ts
│   │   ├── useMultiChainAnnotations.tsx
│   │   ├── useNightingaleComponents.ts
│   │   ├── usePolymerSearch.ts
│   │   ├── useStructureHoverSync.ts
│   │   └── useViewerSync.ts
│   ├── lib
│   │   ├── chain_key.ts
│   │   ├── formatters.ts
│   │   ├── profile_utils.ts
│   │   ├── types
│   │   │   └── annotations.ts
│   │   ├── useDebounce.ts
│   │   ├── utils.ts
│   │   └── utils.tsx
│   ├── services
│   │   └── profile_service.ts
│   ├── store
│   │   ├── emptyApi.ts
│   │   ├── slices
│   │   │   ├── annotationsSlice.ts
│   │   │   ├── chainFocusSlice.ts
│   │   │   ├── colorRulesSelector.ts
│   │   │   ├── sequence_registry.ts
│   │   │   └── slice_structures.ts
│   │   ├── store.ts
│   │   └── tubxz_api.ts
│   └── types
├── tailwind.config.js
├── TODOs.md
├── tsconfig.json
└── yarn.lock
```

One thing to note is that my frontend houses a `nightingale` "fork" -- that is i cloned their repo, tore out the github folder and made some changes to the source code that were necessary for my app and its working ok overall but is giving me hella headaches ever since then on every build. Let's be careful with that. I can show you both package.json files if needed...

Please let me know if you'd like to see any files in particular, ask me questions about the system or the particulars of the architecture. We even have a deepwiki for each of the parts of the application so we can get precise references to the codeblocks.


So far i have this going:
```
ᢹ saeta.rtviii[ dev/tubxz_deployment ]  tree -L 6
.
├── deploy
    .env
    .env.example
│   ├── docker-compose.yml
│   └── nginx
│       ├── certs
│       └── nginx.conf
├── notes.md
├── setup.sh
├── tubulinxyz
│   └── Dockerfile
└── tubulinxyz_fend
    └── Dockerfile

6 directories, 6 files

```
