
Let me show you my project (it's a react/next.js app):
```
ᢹ saeta.rtviii[ dev/fend_tubulinxyz ]  ls -la                                                                                                                    [main]
total 4160
-rw-r--r--@   1 rtviii  staff  636044 Aug  5 23:35 α tubulin interactions – Tubulin Mutations(Sheet1).csv
-rw-r--r--@   1 rtviii  staff  282306 Aug  5 23:35 α tubulin modifications – Tubulin Mutations(Sheet1).csv
-rw-r--r--@   1 rtviii  staff  655624 Aug  5 23:34 α tubulin mutations – Tubulin Mutations(Sheet1).csv
drwxr-xr-x@  25 rtviii  staff     800 Aug  6 00:12 .
drwxr-xr-x@  45 rtviii  staff    1440 Aug  5 23:41 ..
-rw-r--r--@   1 rtviii  staff    6148 Jul 22 00:16 .DS_Store
drwxr-xr-x   14 rtviii  staff     448 Aug  6 10:16 .git
-rw-r--r--    1 rtviii  staff     480 Jun 17 21:05 .gitignore
drwxr-xr-x   11 rtviii  staff     352 Aug 28 17:13 .next
-rw-r--r--    1 rtviii  staff     448 Jun 27 07:25 components.json
drwxr-xr-x    6 rtviii  staff     192 Jun 22 20:25 data
-rw-r--r--    1 rtviii  staff     393 Jun 17 21:05 eslint.config.mjs
-rw-r--r--    1 rtviii  staff     211 Jun 28 18:17 next-env.d.ts
-rw-r--r--    1 rtviii  staff     500 Jun 28 18:57 next.config.js
drwxr-xr-x  509 rtviii  staff   16288 Jul 25 16:12 node_modules
-rw-r--r--    1 rtviii  staff  292342 Jun 28 18:20 package-lock.json
-rw-r--r--    1 rtviii  staff     804 Jul 26 17:17 package.json
-rw-r--r--    1 rtviii  staff     102 Jun 27 06:53 postcss.config.js
-rw-r--r--    1 rtviii  staff      81 Jun 17 21:05 postcss.config.mjs
drwxr-xr-x    7 rtviii  staff     224 Jun 17 21:05 public
-rw-r--r--    1 rtviii  staff    1450 Jun 17 21:05 README.md
drwxr-xr-x    8 rtviii  staff     256 Aug  3 13:34 src
-rw-r--r--    1 rtviii  staff    1849 Jun 28 18:19 tailwind.config.js
-rw-r--r--    1 rtviii  staff     710 Jun 28 18:51 tsconfig.json
-rw-r--r--    1 rtviii  staff  200029 Jul 26 17:17 yarn.lock
ᢹ saeta.rtviii[ dev/fend_tubulinxyz ]  tree src                                                                                                                  [main]
src
├── app
│   ├── data.ts
│   ├── globals.css
│   ├── layout.tsx
│   ├── page.tsx
│   ├── research_panel_tubdb_processor.tsx
│   └── research_panel.tsx
├── components
│   ├── chain_panel.tsx
│   ├── entities_panel.tsx
│   ├── molstar
│   │   ├── colors
│   │   │   └── colorscheme.ts
│   │   ├── molstar_controller.tsx
│   │   ├── molstar_preset_computed_residues.tsx
│   │   ├── molstar_preset.tsx
│   │   ├── molstar_service.tsx
│   │   ├── molstar_spec.tsx
│   │   ├── molstar_viewer.tsx
│   │   ├── mstar.css
│   │   ├── preset-helpers.ts
│   │   └── tubulin-color-theme.tsx
│   ├── nonpolymer_panel.tsx
│   ├── protofilament_grid.tsx
│   ├── sequence_viewer.tsx
│   └── seqviz.styles.css
├── hooks
│   └── useMolstarSync.ts
├── lib
│   ├── utils.ts
│   └── utils.tsx
├── services
│   ├── gql_parser.ts
│   ├── protofilament_grid_parser.tsx
│   └── rcsb_graphql_service.ts
└── store
    ├── slices
    │   ├── interaction_slice.ts
    │   ├── molstar_refs.ts
    │   ├── nonpolymer_states.ts
    │   ├── polymer_states.ts
    │   ├── sequence_structure_sync.ts
    │   ├── sequence_viewer.ts
    │   └── tubulin_structures.ts
    └── store.ts

10 directories, 36 files
```


It's very quickly evolving so please let me kwnow if you want to see any file in particular and i'll attach it to the context.

Currently it's just the main page where we have a single structure of tubulin being uploaded to the molstar viewer where i worked out the main interactions between things: ligand visualization neighborhood visualization, sequence display, chain differentiation and selection, upload etc.  Most of these are facets of the molstar viewer and controller and are most of the interactions that should be there for any visualized tubulin strucutre. Ultimately, in retail terms -- it's an "item" page.


However we now want to "scale" it out in the following ways: 
- basically have a landing page (with description of the project, etc) 
- selection of wether to inspect a curved or straight tubulin structure 
- then a "catalogue" page where the items/available tubulin structures will be listed as cards with some description and filtes (we will just mock up the filters UI for now, no time to implement it currently). Each card would then if clicked link to the item page (the current main page).


    

